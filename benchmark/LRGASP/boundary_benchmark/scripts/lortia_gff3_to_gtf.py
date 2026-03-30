#!/usr/bin/env python3
import argparse
import gzip
import re
import sys
from pathlib import Path
from collections import defaultdict

def open_maybe_gzip(path: Path, mode="rt"):
    p = str(path)
    if p.endswith(".gz"):
        return gzip.open(p, mode, encoding="utf-8", errors="replace")
    return open(p, mode, encoding="utf-8", errors="replace")

def strip_ext(name: str) -> str:
    n = name
    if n.lower().endswith(".gff3.gz"):
        return n[:-8]
    if n.lower().endswith(".gff3"):
        return n[:-5]
    if n.lower().endswith(".gz"):
        return n[:-3]
    return Path(n).stem

def parse_attrs_gff3(attr_field: str) -> dict:
    """
    Handles:
      - standard GFF3: key=value;key2=value2
      - also tolerates single-token IDs (like 'chr10_1') -> {'ID': 'chr10_1'}
    """
    if not attr_field or attr_field == ".":
        return {}

    s = attr_field.strip()

    # LoRTIA TES/TSS case: single token (no ;, no =, no space) => treat as ID
    if (";" not in s) and ("=" not in s) and (" " not in s):
        return {"ID": s}

    d = {}
    for part in s.split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
        elif " " in part:
            k, v = part.split(" ", 1)
        else:
            k, v = part, ""
        d[k.strip()] = v.strip().strip('"')
    return d

def gtf_attr_str(d: dict) -> str:
    # Stable-ish order helps diffs
    keys = ["gene_id", "transcript_id", "exon_number", "tag"]
    out = []
    used = set()
    for k in keys:
        if k in d:
            out.append(f'{k} "{d[k]}";')
            used.add(k)
    for k in sorted(d.keys()):
        if k in used:
            continue
        out.append(f'{k} "{d[k]}";')
    return " ".join(out)

def detect_default_tag(in_path: Path):
    u = in_path.name.upper()
    if "TSS" in u:
        return "TSS"
    if "TES" in u:
        return "TES"
    return None

def gather_inputs(root: Path):
    if root.is_file():
        return [root]
    files = list(root.rglob("*.gff3")) + list(root.rglob("*.gff3.gz"))
    return sorted(files)

def convert_one_gff3_to_gtf(in_path: Path, out_path: Path, overwrite: bool, verbose: bool):
    if out_path.exists() and not overwrite:
        if verbose:
            print(f"[SKIP] exists: {out_path}", file=sys.stderr)
        return

    # First pass: collect records
    genes = set()
    transcripts = {}            # tid -> dict(seqid, source, start, end, strand, gid)
    exons_by_tid = defaultdict(list)  # tid -> list of (start,end)

    point_records = []  # list of dict for tss/tes-like records: (seqid,source,ftype,start,end,strand,id,score)

    default_tag = detect_default_tag(in_path)

    with open_maybe_gzip(in_path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            seqid, source, ftype, start_s, end_s, score, strand, phase, attrs = cols
            try:
                start = int(start_s); end = int(end_s)
            except ValueError:
                continue

            a = parse_attrs_gff3(attrs)
            f = ftype.lower()

            if f == "gene":
                gid = a.get("ID") or a.get("gene_id") or a.get("Name")
                if gid:
                    genes.add(gid)
                continue

            if f in ("mrna", "transcript"):
                tid = a.get("transcript_id") or a.get("ID") or a.get("Name")
                if not tid:
                    tid = f"tx_{seqid}_{start}_{end}_{strand}"
                parent = a.get("Parent")
                gid = a.get("gene_id") or a.get("gene") or (parent if parent in genes else None) or f"gene_{tid}"
                transcripts[tid] = dict(
                    seqid=seqid,
                    source=source or "LoRTIA",
                    start=start,
                    end=end,
                    strand=strand if strand in ("+","-") else ".",
                    gid=gid
                )
                continue

            if f == "exon":
                parent = a.get("Parent") or a.get("transcript_id")
                if not parent:
                    continue
                # Parent can be comma-separated
                for tid in [p.strip() for p in parent.split(",") if p.strip()]:
                    exons_by_tid[tid].append((start, end))
                continue

            # LoRTIA boundary outputs: 'tss' / 'tes' (and often no key=value attrs)
            if f in ("tss", "tes"):
                rid = a.get("ID") or a.get("Name") or f"{seqid}_{start}_{strand}"
                point_records.append(dict(
                    seqid=seqid, source=source or "LoRTIA", ftype=f, start=start, end=end,
                    strand=strand if strand in ("+","-") else ".", rid=rid, score=score
                ))
                continue

    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Decide mode
    has_isoform = len(transcripts) > 0
    has_points = len(point_records) > 0

    # If isoform GFF3: write transcript+exon (Bambu-like: no explicit 'gene' feature lines)
    with open(out_path, "w", encoding="utf-8") as out:

        if has_isoform:
            # deterministic order
            items = sorted(
                transcripts.items(),
                key=lambda kv: (kv[1]["seqid"], kv[1]["start"], kv[1]["end"], kv[1]["strand"], kv[0])
            )

            for tid, info in items:
                seqid = info["seqid"]
                source = info["source"]
                t_start = info["start"]
                t_end = info["end"]
                strand = info["strand"]
                gid = info["gid"]

                exs = exons_by_tid.get(tid, [])
                if not exs:
                    # monoexon fallback
                    exs = [(t_start, t_end)]

                # exon ordering 5'->3'
                if strand == "+":
                    exs_sorted = sorted(exs, key=lambda x: (x[0], x[1]))
                else:
                    exs_sorted = sorted(exs, key=lambda x: (-x[1], -x[0]))

                tr_attrs = {"gene_id": gid, "transcript_id": tid}
                if default_tag:
                    tr_attrs["tag"] = default_tag

                out.write("\t".join([
                    seqid, source, "transcript",
                    str(t_start), str(t_end),
                    ".", strand, ".", gtf_attr_str(tr_attrs)
                ]) + "\n")

                for exon_number, (xs, xe) in enumerate(exs_sorted, start=1):
                    ex_attrs = {
                        "gene_id": gid,
                        "transcript_id": tid,
                        "exon_number": str(exon_number),
                    }
                    out.write("\t".join([
                        seqid, source, "exon",
                        str(xs), str(xe),
                        ".", strand, ".", gtf_attr_str(ex_attrs)
                    ]) + "\n")

        # If boundary points: each record -> transcript line (start=end already in your TES/TSS files)
        if has_points and not has_isoform:
            # keep original order (or sort if you prefer)
            for r in point_records:
                seqid = r["seqid"]
                source = r["source"]
                strand = r["strand"]
                pos_start = r["start"]
                pos_end = r["end"]
                rid = r["rid"]

                # bambu-like ids
                gid = f"gene_{rid}"
                tid = f"tx_{rid}"

                tag = "TSS" if r["ftype"] == "tss" else "TES"
                tr_attrs = {"gene_id": gid, "transcript_id": tid, "tag": tag}

                out.write("\t".join([
                    seqid, source, "transcript",
                    str(pos_start), str(pos_end),
                    ".", strand, ".", gtf_attr_str(tr_attrs)
                ]) + "\n")

        if (not has_isoform) and (not has_points):
            # nothing recognized
            if verbose:
                print(f"[WARN] No convertible features found in: {in_path}", file=sys.stderr)

def main():
    ap = argparse.ArgumentParser(
        description="Convert LoRTIA GFF3 (isoform mRNA/exon OR boundary tss/tes) to Bambu-like GTF. Recursively processes folders."
    )
    ap.add_argument("-i", "--input", required=True, help="Input .gff3/.gff3.gz file OR directory (recursive).")
    ap.add_argument("-o", "--outdir", default=None, help="Optional output root. If omitted: write .gtf next to each .gff3 (inplace).")
    ap.add_argument("--overwrite", action="store_true", help="Overwrite existing .gtf files.")
    ap.add_argument("-v", "--verbose", action="store_true", help="Verbose logging.")
    args = ap.parse_args()

    in_path = Path(args.input).resolve()
    inputs = gather_inputs(in_path)
    if not inputs:
        raise SystemExit(f"Nincs feldolgozható GFF3: {in_path}")

    base = in_path if in_path.is_dir() else in_path.parent
    out_root = Path(args.outdir).resolve() if args.outdir else None

    for gff in inputs:
        if out_root is None:
            # inplace: same dir
            out_path = gff.parent / f"{strip_ext(gff.name)}.gtf"
        else:
            rel = gff.relative_to(base)
            out_path = out_root / rel.parent / f"{strip_ext(rel.name)}.gtf"

        if args.verbose:
            print(f"[INFO] {gff} -> {out_path}", file=sys.stderr)

        convert_one_gff3_to_gtf(gff, out_path, args.overwrite, args.verbose)

if __name__ == "__main__":
    main()