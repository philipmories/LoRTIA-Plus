#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
compare_isoforms_with_windows_clean.py

Compare two TSV files (base vs query) using TSS, TES and intron windows.

- Keeps all columns from the base file
- Adds numeric columns from the query file
- Values are summed for matching isoforms
"""

import argparse
import pandas as pd


# ----------------------------
# Helpers
# ----------------------------

def normalize_columns(df):
    """Convert column names to lowercase for internal consistency."""
    df.columns = [c.lower() for c in df.columns]
    return df


def parse_exon_comp(exon_comp):
    """Convert exon composition string to list of tuples."""
    exons = []
    for part in str(exon_comp).split(";"):
        if not part:
            continue
        s, e = part.split("-")
        exons.append((int(s), int(e)))
    return exons


def get_tss(start, end, strand):
    return start if strand == "+" else end


def get_tes(start, end, strand):
    return end if strand == "+" else start


def get_internal_boundaries(exons):
    """Return internal exon boundaries (splice sites)."""
    boundaries = []
    for i in range(len(exons) - 1):
        boundaries.append(exons[i][1])
        boundaries.append(exons[i + 1][0])
    return boundaries


def within(a, b, window):
    return abs(a - b) <= window


def isoform_match(row1, row2, tss_w, tes_w, intron_w):
    """Check if two isoforms match within given windows."""

    if row1.chromosome != row2.chromosome:
        return False

    if row1.strand != row2.strand:
        return False

    # TSS check
    if not within(
        get_tss(row1.start, row1.end, row1.strand),
        get_tss(row2.start, row2.end, row2.strand),
        tss_w,
    ):
        return False

    # TES check
    if not within(
        get_tes(row1.start, row1.end, row1.strand),
        get_tes(row2.start, row2.end, row2.strand),
        tes_w,
    ):
        return False

    # Intron structure check
    ex1 = parse_exon_comp(row1.exon_composition)
    ex2 = parse_exon_comp(row2.exon_composition)

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
        description="Compare isoforms between TSV files using TSS/TES/intron windows"
    )

    parser.add_argument("--base", required=True, help="Base TSV file")
    parser.add_argument("--query", required=True, help="Query TSV file")
    parser.add_argument("--out", required=True, help="Output TSV file")

    parser.add_argument("--tss-window", type=int, default=10)
    parser.add_argument("--tes-window", type=int, default=10)
    parser.add_argument("--intron-window", type=int, default=5)

    args = parser.parse_args()

    base = pd.read_csv(args.base, sep="\t")
    query = pd.read_csv(args.query, sep="\t")

    base = normalize_columns(base)
    query = normalize_columns(query)

    fixed_cols = {
        "isoform",
        "chromosome",
        "start",
        "end",
        "strand",
        "exon_composition",
    }

    query_value_cols = [c for c in query.columns if c not in fixed_cols]

    # initialize new columns
    for col in query_value_cols:
        base[col] = 0

    # matching loop
    for i, brow in base.iterrows():
        for _, qrow in query.iterrows():

            if isoform_match(
                brow,
                qrow,
                args.tss_window,
                args.tes_window,
                args.intron_window,
            ):
                for col in query_value_cols:
                    base.at[i, col] += qrow[col]

    base.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()