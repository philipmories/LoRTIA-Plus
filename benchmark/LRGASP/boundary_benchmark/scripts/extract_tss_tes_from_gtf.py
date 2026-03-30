#!/usr/bin/env python3
import argparse
import gzip
import sys
from pathlib import Path

def open_maybe_gzip(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")

def gtf_stem(path: Path) -> str:
    name = path.name
    if name.endswith(".gtf.gz"):
        return name[:-7]
    if name.endswith(".gtf"):
        return name[:-4]
    return path.stem

def safe_relpath(p: Path, root: Path) -> Path:
    try:
        return p.relative_to(root)
    except ValueError:
        return Path(p.name)

def append_tag(attrs: str, tag_value: str) -> str:
    a = attrs.rstrip()
    if not a.endswith(";"):
        a += ";"
    a += f' tag "{tag_value}";'
    return a

def process_one_gtf(gtf_path: Path, in_root: Path, out_root: Path,
                    keep_structure: bool, dedup: bool, verbose: bool):
    rel = safe_relpath(gtf_path, in_root)
    subdir = rel.parent if keep_structure else Path(".")
    out_dir = out_root / subdir
    out_dir.mkdir(parents=True, exist_ok=True)

    stem = gtf_stem(gtf_path)
    tss_out = out_dir / f"{stem}.TSS.gtf"
    tes_out = out_dir / f"{stem}.TES.gtf"

    seen_tss = set()  # (chrom, strand, tss_pos)
    seen_tes = set()  # (chrom, strand, tes_pos)

    kept_tss = kept_tes = 0
    dup_tss = dup_tes = 0
    skipped = 0
    total_tx = 0

    with open_maybe_gzip(gtf_path) as fh, \
         open(tss_out, "w", encoding="utf-8") as tss_fh, \
         open(tes_out, "w", encoding="utf-8") as tes_fh:

        for line in fh:
            if not line:
                continue

            if line.startswith("#"):
                tss_fh.write(line)
                tes_fh.write(line)
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue

            chrom, source, feature, start_s, end_s, score, strand, frame, attrs = cols
            if feature != "transcript":
                continue
            if strand not in {"+", "-"}:
                skipped += 1
                continue

            try:
                start = int(start_s)
                end = int(end_s)
            except ValueError:
                skipped += 1
                continue

            total_tx += 1

            # TSS/TES pozíció (1-based) strand szerint
            if strand == "+":
                tss_pos = start
                tes_pos = end
            else:
                tss_pos = end
                tes_pos = start

            # dedup kulcsok
            tss_key = (chrom, strand, tss_pos)
            tes_key = (chrom, strand, tes_pos)

            # TSS sor: start=end=tss_pos (GTF 1-based)
            if (not dedup) or (tss_key not in seen_tss):
                if dedup:
                    seen_tss.add(tss_key)
                tss_cols = cols[:]
                tss_cols[3] = str(tss_pos)
                tss_cols[4] = str(tss_pos)
                tss_cols[8] = append_tag(attrs, "TSS")
                tss_fh.write("\t".join(tss_cols) + "\n")
                kept_tss += 1
            else:
                dup_tss += 1

            # TES sor: start=end=tes_pos
            if (not dedup) or (tes_key not in seen_tes):
                if dedup:
                    seen_tes.add(tes_key)
                tes_cols = cols[:]
                tes_cols[3] = str(tes_pos)
                tes_cols[4] = str(tes_pos)
                tes_cols[8] = append_tag(attrs, "TES")
                tes_fh.write("\t".join(tes_cols) + "\n")
                kept_tes += 1
            else:
                dup_tes += 1

    if verbose:
        print(
            f"[OK] {gtf_path}\n"
            f"  transcripts_total={total_tx} skipped={skipped}\n"
            f"  -> {tss_out} kept={kept_tss}" + (f" dups_removed={dup_tss}" if dedup else "") + "\n"
            f"  -> {tes_out} kept={kept_tes}" + (f" dups_removed={dup_tes}" if dedup else ""),
            file=sys.stderr
        )

def main():
    ap = argparse.ArgumentParser(
        description="Recursively extract TSS/TES point coordinates from transcript features in GTF, write per-file .TSS.gtf/.TES.gtf (STRICTLY GTF)."
    )
    ap.add_argument("-i", "--input", required=True, help="Input folder with .gtf/.gtf.gz files (recursively).")
    ap.add_argument("-o", "--output", required=True, help="Output folder for *.TSS.gtf and *.TES.gtf files.")
    ap.add_argument("--no-structure", action="store_true",
                    help="Do NOT preserve input subfolder structure; dump all outputs into output root (risk of name collisions).")
    ap.add_argument("--no-dedup", action="store_true", help="Disable deduplication (write all points).")
    ap.add_argument("-v", "--verbose", action="store_true", help="Verbose logging to stderr.")
    args = ap.parse_args()

    in_root = Path(args.input).resolve()
    out_root = Path(args.output).resolve()
    out_root.mkdir(parents=True, exist_ok=True)

    if not in_root.exists() or not in_root.is_dir():
        print(f"ERROR: input path is not a directory: {in_root}", file=sys.stderr)
        sys.exit(2)

    gtfs = sorted(list(in_root.rglob("*.gtf")) + list(in_root.rglob("*.gtf.gz")))
    if not gtfs:
        print(f"WARNING: no .gtf or .gtf.gz files found under: {in_root}", file=sys.stderr)
        sys.exit(0)

    keep_structure = not args.no_structure
    dedup = not args.no_dedup

    for gtf_path in gtfs:
        process_one_gtf(gtf_path, in_root, out_root, keep_structure, dedup, args.verbose)

if __name__ == "__main__":
    main()