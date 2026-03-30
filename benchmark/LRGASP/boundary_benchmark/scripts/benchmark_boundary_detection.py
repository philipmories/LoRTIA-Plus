#!/usr/bin/env python3
import argparse
import gzip
import sys
import re
from pathlib import Path
from collections import defaultdict
from bisect import bisect_left
import heapq

ANNOTATORS_DEFAULT = ["BAMBU", "IsoQuant", "FLAIR", "LoRTIA", "NAGATA"]

def open_maybe_gzip(path: Path, mode="rt"):
    p = str(path)
    if p.endswith(".gz"):
        return gzip.open(p, mode, encoding="utf-8", errors="replace")
    return open(p, mode, encoding="utf-8", errors="replace")

def split_ws(line: str):
    # TAB/space + CRLF tolerant
    return re.split(r"\s+", line.strip())

def norm_contig(contig: str, chr_mode: str) -> str:
    if chr_mode == "keep":
        return contig
    if chr_mode == "add":
        return contig if contig.startswith("chr") else "chr" + contig
    if chr_mode == "strip":
        return contig[3:] if contig.startswith("chr") else contig
    raise ValueError("bad chr_mode")

def infer_annotator_from_path(p: Path, annotators):
    s = p.name
    for a in annotators:
        if re.search(rf"(^|[-_.]){re.escape(a)}($|[-_.])", s, flags=re.IGNORECASE):
            return a
    # fallback: try parent folder names
    for part in p.parts[::-1]:
        for a in annotators:
            if part.lower() == a.lower():
                return a
    return "."

def parse_manifest_4col(manifest: Path, verbose: bool):
    """
    Columns (TAB): Chemistry  Cell-line  GTF  Reference
    Header allowed (first row).
    """
    rows = []
    with open_maybe_gzip(manifest, "rt") as fh:
        for ln, line in enumerate(fh, start=1):
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = split_ws(s)
            if len(parts) < 4:
                if verbose:
                    print(f"[WARN] manifest line {ln}: need 4 cols, skipping", file=sys.stderr)
                continue

            # header detection
            if ln == 1 and parts[0].lower().startswith("chem") and parts[2].lower() in {"gtf", "path"}:
                continue

            chem, cell, gtf_p, ref_p = parts[0], parts[1], parts[2], parts[3]
            rows.append((chem, cell, Path(gtf_p), Path(ref_p)))

    if not rows:
        raise SystemExit("ERROR: manifest parsed 0 rows.")
    return rows

def load_ref_intervals(ref_bed: Path, mode: str, window: int, ref_coord: str,
                       ref_point_mode: str, chr_mode: str):
    """
    Returns dict: (chrom, strand) -> list of (lo, hi) inclusive (1-based)
    TSS: ref is cluster BED6. bed0 => [start+1, end], bed1 => [start, end], then ±window
    TES: ref is point BED6. as_is => pos=start, bed0_to_1based => pos=start+1, then ±window
    """
    d = defaultdict(list)
    with open_maybe_gzip(ref_bed, "rt") as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            cols = split_ws(line)
            if len(cols) < 6:
                continue

            chrom = norm_contig(cols[0], chr_mode)
            try:
                a = int(cols[1]); b = int(cols[2])
            except ValueError:
                continue
            strand = cols[5].strip()
            if strand not in {"+", "-"}:
                continue

            if mode == "TSS":
                if ref_coord == "bed0":
                    lo = (a + 1) - window
                    hi = b + window
                else:  # bed1
                    lo = a - window
                    hi = b + window
            else:  # TES
                if ref_point_mode == "as_is":
                    pos = a
                else:  # bed0_to_1based
                    pos = a + 1
                lo = pos - window
                hi = pos + window

            if lo < 1:
                lo = 1
            if hi < lo:
                continue

            d[(chrom, strand)].append((lo, hi))

    out = {}
    for k, intervals in d.items():
        out[k] = sorted(set(intervals), key=lambda x: (x[0], x[1]))
    return out

def load_pred_points_from_gtf(gtf_path: Path, mode: str, chr_mode: str):
    """
    Reads transcript features only.
    TSS: + start, - end
    TES: + end,   - start
    Returns dict: (chrom, strand) -> sorted unique positions (1-based)
    """
    pos_sets = defaultdict(set)
    with open_maybe_gzip(gtf_path, "rt") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            chrom, _, feature, start_s, end_s, _, strand, _, _ = cols
            if feature != "transcript":
                continue
            if strand not in {"+", "-"}:
                continue
            try:
                start = int(start_s); end = int(end_s)
            except ValueError:
                continue

            chrom = norm_contig(chrom, chr_mode)
            if mode == "TSS":
                pos = start if strand == "+" else end
            else:
                pos = end if strand == "+" else start

            pos_sets[(chrom, strand)].add(pos)

    return {k: sorted(v) for k, v in pos_sets.items()}

def max_match_points_to_intervals(points, intervals):
    """
    Maximum 1-to-1 matching between point set and interval set.
    Sweep points; keep candidate intervals in a min-heap by interval end.
    """
    if not points or not intervals:
        return 0
    tp = 0
    j = 0
    heap = []
    n = len(intervals)

    for p in points:
        while j < n and intervals[j][0] <= p:
            heapq.heappush(heap, intervals[j][1])
            j += 1
        while heap and heap[0] < p:
            heapq.heappop(heap)
        if heap:
            heapq.heappop(heap)
            tp += 1
    return tp

def compute_metrics(tp, n_pred, n_ref):
    fp = n_pred - tp
    fn = n_ref - tp
    prec = tp / (tp + fp) if (tp + fp) else 0.0
    rec = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = (2 * prec * rec / (prec + rec)) if (prec + rec) else 0.0
    return fp, fn, prec, rec, f1

def main():
    ap = argparse.ArgumentParser(description="Benchmark boundary calls (TSS/TES) using 4-col manifest (Chemistry, Cell-line, GTF, Reference).")
    ap.add_argument("-m", "--manifest", required=True, help="TSV: Chemistry  Cell-line  GTF  Reference")
    ap.add_argument("--mode", choices=["TSS", "TES"], required=True)
    ap.add_argument("--window", type=int, default=0, help="Tolerance window in nt (default 0).")
    ap.add_argument("--ref-coord", choices=["bed0", "bed1"], default="bed0",
                    help="TSS only: interpret ref BED columns as bed0 (default) or bed1.")
    ap.add_argument("--ref-point-mode", choices=["as_is", "bed0_to_1based"], default="as_is",
                    help="TES only: interpret point ref start=end as-is (default) or bed0->1based.")
    ap.add_argument("--chr-mode", choices=["keep", "add", "strip"], default="keep",
                    help="Normalize contigs: keep/add/strip 'chr' between ref and pred.")
    ap.add_argument("--annotators", nargs="*", default=ANNOTATORS_DEFAULT, help="Used only to infer annotator name from GTF filename.")
    ap.add_argument("-o", "--output", default="-", help="Output TSV (default stdout).")
    ap.add_argument("--strict", action="store_true", help="Fail if any GTF/Reference file is missing.")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    if args.window < 0:
        raise SystemExit("ERROR: --window must be >= 0")

    manifest = Path(args.manifest).resolve()
    if not manifest.exists():
        raise SystemExit(f"ERROR: manifest not found: {manifest}")

    rows = parse_manifest_4col(manifest, args.verbose)

    # cache refs by exact ref path (same file repeated 5x annotátoronként)
    ref_cache = {}

    out_fh = sys.stdout if args.output == "-" else open(args.output, "w", encoding="utf-8")
    try:
        out_fh.write("\t".join([
            "mode","window","chem","cell","annotator",
            "tp","fp","fn","precision","recall","f1",
            "n_pred","n_ref","gtf","ref"
        ]) + "\n")

        for chem, cell, gtf_p, ref_p in rows:
            gtf_p = Path(str(gtf_p).strip())
            ref_p = Path(str(ref_p).strip())
            annot = infer_annotator_from_path(gtf_p, args.annotators)

            if not gtf_p.exists():
                msg = f"[MISSING GTF] {gtf_p}"
                if args.strict:
                    raise SystemExit(msg)
                if args.verbose:
                    print(msg, file=sys.stderr)
                continue
            if not ref_p.exists():
                msg = f"[MISSING REF] {ref_p}"
                if args.strict:
                    raise SystemExit(msg)
                if args.verbose:
                    print(msg, file=sys.stderr)
                continue

            ref_key = (ref_p.resolve(), args.mode, args.window, args.ref_coord, args.ref_point_mode, args.chr_mode)
            if ref_key not in ref_cache:
                ref_cache[ref_key] = load_ref_intervals(
                    ref_p, args.mode, args.window, args.ref_coord, args.ref_point_mode, args.chr_mode
                )
                if args.verbose:
                    n_ref = sum(len(v) for v in ref_cache[ref_key].values())
                    print(f"[INFO] loaded ref: {ref_p} n_ref={n_ref}", file=sys.stderr)

            ref_intervals = ref_cache[ref_key]
            pred_points = load_pred_points_from_gtf(gtf_p, args.mode, args.chr_mode)

            n_ref = sum(len(v) for v in ref_intervals.values())
            n_pred = sum(len(v) for v in pred_points.values())

            tp = 0
            for key, points in pred_points.items():
                intervals = ref_intervals.get(key, [])
                tp += max_match_points_to_intervals(points, intervals)

            fp, fn, prec, rec, f1 = compute_metrics(tp, n_pred, n_ref)

            out_fh.write("\t".join([
                args.mode, str(args.window), chem, cell, annot,
                str(tp), str(fp), str(fn),
                f"{prec:.6f}", f"{rec:.6f}", f"{f1:.6f}",
                str(n_pred), str(n_ref),
                str(gtf_p), str(ref_p)
            ]) + "\n")

            if args.verbose:
                print(f"[OK] {chem}×{cell} {annot} TP={tp} FP={fp} FN={fn} F1={f1:.4f}", file=sys.stderr)

    finally:
        if out_fh is not sys.stdout:
            out_fh.close()

if __name__ == "__main__":
    main()