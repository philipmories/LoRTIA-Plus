#!/usr/bin/env python3
import argparse
import csv
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

TX_RE = re.compile(r'transcript_id "([^"]+)"')
GENE_RE = re.compile(r'gene_id "([^"]+)"')

FSM_ISM_CANON = {"full-splice_match", "incomplete-splice_match"}
FSM_ISM_ALIAS = {"fsm", "ism"}  # if someone exported short labels

def open_maybe_gzip(path: Path, mode: str = "rt"):
    p = str(path)
    if p.endswith(".gz"):
        return gzip.open(p, mode, encoding="utf-8", errors="replace", newline="")
    return open(p, mode, encoding="utf-8", errors="replace", newline="")

def norm_cat(x: str) -> str:
    return (x or "").strip()

def cat_is_fsm_ism(struct_cat: str) -> bool:
    c = norm_cat(struct_cat)
    if not c:
        return False
    c_low = c.lower()
    if c in FSM_ISM_CANON:
        return True
    if c_low in FSM_ISM_ALIAS:
        return True
    # tolerate minor variants
    c_low = c_low.replace("_", "-")
    if c_low in {"full-splice-match", "incomplete-splice-match"}:
        return True
    return False

def strip_version(tx: str) -> str:
    # ENST00000.2 -> ENST00000
    return re.sub(r"\.\d+$", "", tx)

def parse_manifest(manifest: Path) -> List[Tuple[str, str, Path]]:
    """
    Manifest TSV columns: Chemistry, Cell-line, GTF (path to classification.txt)
    """
    rows = []
    with open_maybe_gzip(manifest, "rt") as fh:
        # detect delimiter as tab (your file is TSV)
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit("ERROR: manifest has no header row.")

        # normalize column names
        cols = {c.strip(): c for c in reader.fieldnames}
        # required (your header exactly: Chemistry, Cell-line, GTF)
        req = ["Chemistry", "Cell-line", "GTF"]
        missing = [r for r in req if r not in cols]
        if missing:
            raise SystemExit(f"ERROR: manifest missing required columns: {', '.join(missing)}")

        for i, rec in enumerate(reader, start=2):
            chem = (rec.get(cols["Chemistry"]) or "").strip()
            cell = (rec.get(cols["Cell-line"]) or "").strip()
            p = (rec.get(cols["GTF"]) or "").strip()
            if not chem or not cell or not p:
                raise SystemExit(f"ERROR: manifest row {i} has empty Chemistry/Cell-line/GTF.")
            rows.append((chem, cell, Path(p)))
    return rows

def read_sqanti_fsm_ism_associated_transcripts(
    classification_path: Path,
    strip_tx_version: bool
) -> Set[str]:
    """
    Reads SQANTI3 classification.txt (tab-delimited).
    Keeps rows with structural_category in FSM/ISM,
    returns set of associated_transcript IDs.
    """
    keep: Set[str] = set()
    with open_maybe_gzip(classification_path, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit(f"ERROR: classification has no header: {classification_path}")
        # expected columns
        if "structural_category" not in reader.fieldnames or "associated_transcript" not in reader.fieldnames:
            raise SystemExit(
                f"ERROR: classification missing 'structural_category' or 'associated_transcript': {classification_path}"
            )

        for rec in reader:
            sc = rec.get("structural_category", "")
            if not cat_is_fsm_ism(sc):
                continue
            at = (rec.get("associated_transcript", "") or "").strip()
            if not at or at.upper() in {"NA", "NAN"} or at == ".":
                continue
            if strip_tx_version:
                at = strip_version(at)
            keep.add(at)
    return keep

def gtf_get_ids(attrs: str) -> Tuple[Optional[str], Optional[str]]:
    tx = None
    gene = None
    m = TX_RE.search(attrs)
    if m:
        tx = m.group(1)
    m = GENE_RE.search(attrs)
    if m:
        gene = m.group(1)
    return tx, gene

def build_cell_active_gtf(
    gtf_path: Path,
    out_path: Path,
    keep_tx: Set[str],
    strip_tx_version: bool,
    verbose: bool
) -> Dict[str, int]:
    """
    Streams GTF and writes subset:
      - any line with transcript_id in keep_tx
      - gene lines (no transcript_id) are kept if gene_id belongs to any kept transcript
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)

    written = 0
    kept_gene_ids: Set[str] = set()

    # First pass (light): collect gene_ids for kept transcripts
    if verbose:
        print(f"[INFO] First pass: collecting gene_ids from {gtf_path}", file=sys.stderr)

    with open_maybe_gzip(gtf_path, "rt") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            attrs = parts[8]
            tx, gene = gtf_get_ids(attrs)
            if tx is None or gene is None:
                continue
            tx_cmp = strip_version(tx) if strip_tx_version else tx
            if tx_cmp in keep_tx:
                kept_gene_ids.add(gene)

    # Second pass: write lines
    if verbose:
        print(f"[INFO] Second pass: writing subset to {out_path}", file=sys.stderr)

    kept_tx_seen: Set[str] = set()
    with open_maybe_gzip(gtf_path, "rt") as fh, open(out_path, "w", encoding="utf-8", newline="") as out:
        for line in fh:
            if not line:
                continue
            if line.startswith("#"):
                out.write(line)
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            attrs = parts[8]
            tx, gene = gtf_get_ids(attrs)

            keep_line = False
            if tx is not None:
                tx_cmp = strip_version(tx) if strip_tx_version else tx
                if tx_cmp in keep_tx:
                    keep_line = True
                    kept_tx_seen.add(tx_cmp)
            else:
                # gene line etc: keep if gene_id is in kept_gene_ids
                if gene is not None and gene in kept_gene_ids:
                    keep_line = True

            if keep_line:
                out.write(line if line.endswith("\n") else line + "\n")
                written += 1

    return {
        "lines_written": written,
        "kept_tx_count": len(keep_tx),
        "kept_tx_found_in_gtf": len(kept_tx_seen),
        "kept_gene_count": len(kept_gene_ids),
    }

def main():
    ap = argparse.ArgumentParser(
        description="Build chemistry×cell 'cell_active' reference GTFs from SQANTI3 classification.txt (FSM/ISM only) by subsetting a GENCODE v48 GTF."
    )
    ap.add_argument("-m", "--manifest", required=True, help="Manifest TSV with columns: Chemistry, Cell-line, GTF (path to SQANTI classification.txt).")
    ap.add_argument("-g", "--gencode-gtf", required=True, help="GENCODE v48 GTF (can be .gz).")
    ap.add_argument("-o", "--output", required=True, help="Output directory for chemistry×cell reference GTFs.")
    ap.add_argument("--strip-tx-version", action="store_true",
                    help="Strip transcript version suffix (e.g., ENST... .2) in BOTH SQANTI and GTF matching. Use only if your IDs don't match by version.")
    ap.add_argument("--keep-structure", action="store_true",
                    help="If set: write outputs under output/<Chemistry>/ ... Otherwise write all into output root.")
    ap.add_argument("--prefix", default="", help="Optional prefix for output filenames.")
    ap.add_argument("--suffix", default="cell_active.gtf", help="Output filename suffix (default: cell_active.gtf).")
    ap.add_argument("--strict", action="store_true",
                    help="Fail if a chemistry×cell group yields 0 kept transcripts (default: still writes an empty-ish GTF with headers).")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    manifest = Path(args.manifest).resolve()
    gtf = Path(args.gencode_gtf).resolve()
    out_root = Path(args.output).resolve()
    out_root.mkdir(parents=True, exist_ok=True)

    if not manifest.exists():
        raise SystemExit(f"ERROR: manifest not found: {manifest}")
    if not gtf.exists():
        raise SystemExit(f"ERROR: GTF not found: {gtf}")

    rows = parse_manifest(manifest)

    # check existence of ALL classification files first (strict)
    missing = [str(p) for _, _, p in rows if not Path(str(p).strip()).exists()]
    if missing:
        msg = "ERROR: Some classification files listed in manifest do not exist:\n" + "\n".join(missing[:50])
        if len(missing) > 50:
            msg += f"\n... plus {len(missing)-50} more"
        raise SystemExit(msg)

    # group rows by (chem, cell)
    groups: Dict[Tuple[str, str], List[Path]] = defaultdict(list)
    for chem, cell, p in rows:
        groups[(chem, cell)].append(Path(str(p).strip()))

    # process each group
    for (chem, cell), paths in sorted(groups.items()):
        keep_tx: Set[str] = set()
        for p in paths:
            keep_tx |= read_sqanti_fsm_ism_associated_transcripts(p, args.strip_tx_version)

        if args.verbose:
            print(f"[INFO] {chem}×{cell}: classification_files={len(paths)} FSM/ISM associated_transcripts={len(keep_tx)}", file=sys.stderr)

        if args.strict and len(keep_tx) == 0:
            raise SystemExit(f"ERROR: {chem}×{cell}: 0 FSM/ISM associated_transcripts found. Check SQANTI categories/columns.")

        out_dir = (out_root / chem) if args.keep_structure else out_root
        fname = f"{args.prefix}{chem}-{cell}-{args.suffix}"
        out_path = out_dir / fname

        stats = build_cell_active_gtf(
            gtf_path=gtf,
            out_path=out_path,
            keep_tx=keep_tx,
            strip_tx_version=args.strip_tx_version,
            verbose=args.verbose
        )

        if args.verbose:
            print(
                f"[OK] wrote {out_path}\n"
                f"     lines_written={stats['lines_written']} kept_tx={stats['kept_tx_count']} "
                f"kept_tx_found_in_gtf={stats['kept_tx_found_in_gtf']} kept_genes={stats['kept_gene_count']}",
                file=sys.stderr
            )

if __name__ == "__main__":
    main()