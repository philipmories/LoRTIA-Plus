#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compare_category_with_windows_clean.py

Match isoforms between two TSV files using TSS, TES, and intron boundary windows.

For each row in the base TSV:
- finds matching transcripts in the query TSV
- appends:
    - Matched_Transcript_IDs
    - Matched_Categories

Column names are case-insensitive and flexible.
"""

import argparse
import pandas as pd
from typing import List, Tuple, Dict


# ----------------------------
# Column handling
# ----------------------------

def _lower_cols(df: pd.DataFrame) -> Dict[str, str]:
    mapping = {}
    for col in df.columns:
        key = str(col).strip().lower()
        if key not in mapping:
            mapping[key] = col
    return mapping


def require_column(df: pd.DataFrame, candidates: List[str], label: str) -> str:
    mapping = _lower_cols(df)
    for cand in candidates:
        key = cand.lower().strip()
        if key in mapping:
            return mapping[key]

    raise ValueError(
        f"Missing required column ({label}). Allowed: {candidates}. Found: {list(df.columns)}"
    )


# ----------------------------
# Isoform matching logic
# ----------------------------

def parse_exon_comp(exon_comp: str) -> List[Tuple[int, int]]:
    exon_comp = str(exon_comp).strip()
    if exon_comp in {"", "nan", "None"}:
        return []

    exons = []
    for part in exon_comp.split(";"):
        part = part.strip()
        if not part:
            continue
        s, e = part.split("-")
        exons.append((int(s), int(e)))
    return exons


def get_tss(start: int, end: int, strand: str) -> int:
    return start if strand == "+" else end


def get_tes(start: int, end: int, strand: str) -> int:
    return end if strand == "+" else start


def get_internal_boundaries(exons: List[Tuple[int, int]]) -> List[int]:
    boundaries = []
    for i in range(len(exons) - 1):
        boundaries.append(exons[i][1])
        boundaries.append(exons[i + 1][0])
    return boundaries


def within(a: int, b: int, window: int) -> bool:
    return abs(a - b) <= window


def isoform_match(
    chrom1, start1, end1, strand1, excomp1,
    chrom2, start2, end2, strand2, excomp2,
    tss_w, tes_w, intron_w
) -> bool:

    if chrom1 != chrom2:
        return False
    if strand1 != strand2:
        return False

    # TSS / TES check
    if not within(get_tss(start1, end1, strand1), get_tss(start2, end2, strand2), tss_w):
        return False
    if not within(get_tes(start1, end1, strand1), get_tes(start2, end2, strand2), tes_w):
        return False

    # intron structure check
    ex1 = parse_exon_comp(excomp1)
    ex2 = parse_exon_comp(excomp2)

    b1 = get_internal_boundaries(ex1)
    b2 = get_internal_boundaries(ex2)

    if not b1 and not b2:
        return True
    if len(b1) != len(b2):
        return False

    for x, y in zip(b1, b2):
        if not within(x, y, intron_w):
            return False

    return True


# ----------------------------
# Main
# ----------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Compare isoforms between TSV files using TSS/TES/intron windows."
    )

    parser.add_argument("--base", required=True, help="Base TSV")
    parser.add_argument("--query", required=True, help="Query TSV")
    parser.add_argument("--out", required=True, help="Output TSV")

    parser.add_argument("--tss-window", type=int, default=10)
    parser.add_argument("--tes-window", type=int, default=10)
    parser.add_argument("--intron-window", type=int, default=5)

    parser.add_argument("--id-out-col", default="Matched_Transcript_IDs")
    parser.add_argument("--cat-out-col", default="Matched_Categories")

    args = parser.parse_args()

    base = pd.read_csv(args.base, sep="\t")
    query = pd.read_csv(args.query, sep="\t")

    # resolve columns
    b_chrom = require_column(base, ["chromosome", "chrom", "seqid"], "chromosome")
    b_start = require_column(base, ["start"], "start")
    b_end   = require_column(base, ["end", "stop"], "end")
    b_str   = require_column(base, ["strand"], "strand")
    b_exon  = require_column(base, ["exon_composition", "exon_comp"], "exon_composition")

    q_id    = require_column(query, ["transcript_id", "isoform", "id"], "Transcript_ID")
    q_chrom = require_column(query, ["chromosome", "chrom", "seqid"], "chromosome")
    q_start = require_column(query, ["start"], "start")
    q_end   = require_column(query, ["end", "stop"], "end")
    q_str   = require_column(query, ["strand"], "strand")
    q_exon  = require_column(query, ["exon_composition", "exon_comp"], "exon_composition")
    q_cat   = require_column(query, ["category"], "category")

    # convert coords
    for df, c in [(base, b_start), (base, b_end), (query, q_start), (query, q_end)]:
        df[c] = pd.to_numeric(df[c], errors="coerce").astype("Int64")

    base[args.id_out_col] = "0"
    base[args.cat_out_col] = "0"

    for idx, brow in base.iterrows():
        if pd.isna(brow[b_start]) or pd.isna(brow[b_end]):
            continue

        chrom1 = str(brow[b_chrom])
        strand1 = str(brow[b_str]).strip()
        start1 = int(brow[b_start])
        end1 = int(brow[b_end])
        excomp1 = str(brow[b_exon])

        matched_ids = []
        matched_cats = []

        for _, qrow in query.iterrows():
            if pd.isna(qrow[q_start]) or pd.isna(qrow[q_end]):
                continue

            if isoform_match(
                chrom1, start1, end1, strand1, excomp1,
                str(qrow[q_chrom]), int(qrow[q_start]), int(qrow[q_end]),
                str(qrow[q_str]).strip(), str(qrow[q_exon]),
                args.tss_window, args.tes_window, args.intron_window
            ):
                matched_ids.append(str(qrow[q_id]))
                matched_cats.append(str(qrow[q_cat]))

        if matched_ids:
            base.at[idx, args.id_out_col] = ",".join(matched_ids)
            base.at[idx, args.cat_out_col] = ",".join(matched_cats)

    base.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()