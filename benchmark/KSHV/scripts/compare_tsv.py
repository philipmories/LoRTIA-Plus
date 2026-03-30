#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compare TSV/GFF/GFF3 files by TSS or TES positions with strand-aware coordinate handling.

The script:
- accepts TSV or GFF/GFF3 as the first input
- accepts TSV or GFF/GFF3 as the second input
- derives TSS or TES positions from start/end and strand
- counts how many matching positions from the second file fall within a given window
- writes the result as a TSV with one additional count column
"""

import argparse
import os
import re
from bisect import bisect_left, bisect_right

import pandas as pd


def basename_no_ext(path: str) -> str:
    """
    Return the filename without common annotation/table extensions.
    """
    base = os.path.basename(path)
    for ext in [".gff3", ".gff", ".tsv", ".txt", ".csv"]:
        if base.lower().endswith(ext):
            return base[:-len(ext)]
    return os.path.splitext(base)[0]


def detect_format(path: str) -> str:
    """
    Detect whether the input should be treated as GFF/GFF3 or TSV-like text.
    """
    lower = path.lower()
    if lower.endswith(".gff") or lower.endswith(".gff3"):
        return "gff"
    return "tsv"


def parse_gff_attributes(attr: str) -> dict:
    """
    Parse GFF/GTF-style attribute strings into a dictionary.
    """
    result = {}
    if not isinstance(attr, str):
        return result

    parts = [part for part in attr.strip().strip(";").split(";") if part.strip()]
    for part in parts:
        part = part.strip()
        if "=" in part:
            key, value = part.split("=", 1)
            result[key.strip()] = value.strip()
        else:
            match = re.match(r'^\s*([A-Za-z0-9_]+)\s+"([^"]+)"\s*$', part)
            if match:
                result[match.group(1)] = match.group(2)

    return result


def tss_tes_from_start_end(start: int, end: int, strand: str):
    """
    Convert genomic start/end coordinates into strand-aware TSS/TES coordinates.
    """
    lo = min(int(start), int(end))
    hi = max(int(start), int(end))

    if strand == "+":
        return lo, hi
    if strand == "-":
        return hi, lo

    return lo, hi


def pick_site(start: int, end: int, strand: str, mode: str) -> int:
    """
    Return the TSS or TES coordinate for one record.
    """
    tss, tes = tss_tes_from_start_end(start, end, strand)
    return tss if mode.upper() == "TSS" else tes


def read_first_file(path: str, file_format: str) -> pd.DataFrame:
    """
    Read the first input file and preserve all its columns for output.
    """
    if file_format == "tsv":
        df = pd.read_csv(path, sep="\t", dtype=str)
        required = {"seqnames", "strand", "start", "end"}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(
                f"Missing required column(s) in the first TSV file: {sorted(missing)}"
            )
        return df

    df = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        header=None,
        dtype=str,
        names=["seqnames", "source", "type", "start", "end", "score", "strand", "phase", "attributes"],
    )

    ids = []
    for attr in df["attributes"].fillna(""):
        attr_dict = parse_gff_attributes(attr)
        ids.append(
            attr_dict.get("ID")
            or attr_dict.get("transcript_id")
            or attr_dict.get("Name")
            or attr_dict.get("gene_id")
            or ""
        )

    df.insert(0, "TR_ID", ids)

    if "All_reads" not in df.columns:
        df.insert(4, "All_reads", "")

    return df


def read_second_positions(path: str, file_format: str, mode: str):
    """
    Return a dictionary:
        (seqnames, strand) -> sorted list of TSS/TES positions
    """
    position_dict = {}

    if file_format == "tsv":
        df = pd.read_csv(path, sep="\t", dtype=str)
        required = {"seqnames", "strand", "start", "end"}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(
                f"Missing required column(s) in the second TSV file: {sorted(missing)}"
            )

        for _, row in df.iterrows():
            seqname = str(row["seqnames"])
            strand = str(row["strand"])

            if strand not in {"+", "-"}:
                continue

            try:
                start = int(row["start"])
                end = int(row["end"])
            except Exception:
                continue

            position = pick_site(start, end, strand, mode)
            position_dict.setdefault((seqname, strand), []).append(position)

    else:
        df = pd.read_csv(
            path,
            sep="\t",
            comment="#",
            header=None,
            dtype=str,
            names=["seqnames", "source", "type", "start", "end", "score", "strand", "phase", "attributes"],
        )

        for _, row in df.iterrows():
            seqname = str(row["seqnames"])
            strand = str(row["strand"])

            if strand not in {"+", "-"}:
                continue

            try:
                start = int(row["start"])
                end = int(row["end"])
            except Exception:
                continue

            position = pick_site(start, end, strand, mode)
            position_dict.setdefault((seqname, strand), []).append(position)

    for key in list(position_dict.keys()):
        position_dict[key].sort()

    return position_dict


def count_matches(sorted_positions, query_pos: int, window: int) -> int:
    """
    Count how many positions fall within [query_pos - window, query_pos + window].
    """
    if not sorted_positions:
        return 0

    left = query_pos - window
    right = query_pos + window
    i = bisect_left(sorted_positions, left)
    j = bisect_right(sorted_positions, right)
    return j - i


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Compare TSV/GFF files by TSS or TES positions using strand-aware start/end interpretation. "
            "The output column contains 0 if no match is found, otherwise the number of matches within the window."
        )
    )
    parser.add_argument("-i", "--input1", required=True, help="First input file: TSV or GFF/GFF3.")
    parser.add_argument("-g", "--input2", required=True, help="Second input file: GFF/GFF3 or TSV.")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file.")
    parser.add_argument("--mode", choices=["TSS", "TES"], required=True, help="Compare TSS or TES positions.")
    parser.add_argument(
        "--window",
        type=int,
        default=10,
        help="Window size in bp (±window). Default: 10.",
    )
    parser.add_argument(
        "--ignore_strand",
        action="store_true",
        help="Search for matches regardless of strand.",
    )

    args = parser.parse_args()

    format1 = detect_format(args.input1)
    format2 = detect_format(args.input2)

    df1 = read_first_file(args.input1, format1)

    if not {"seqnames", "strand", "start", "end"}.issubset(df1.columns):
        raise ValueError(
            "The first input must contain the following columns: seqnames, strand, start, end."
        )

    positions2 = read_second_positions(args.input2, format2, args.mode)
    new_column = f"{basename_no_ext(args.input2)}_{args.mode.upper()}"

    counts = []
    for _, row in df1.iterrows():
        seqname = str(row.get("seqnames", ""))
        strand = str(row.get("strand", ""))

        try:
            start = int(row.get("start"))
            end = int(row.get("end"))
        except Exception:
            counts.append(0)
            continue

        site = pick_site(start, end, strand, args.mode)

        if args.ignore_strand:
            match_count = (
                count_matches(positions2.get((seqname, "+"), []), site, args.window)
                + count_matches(positions2.get((seqname, "-"), []), site, args.window)
            )
        else:
            match_count = count_matches(positions2.get((seqname, strand), []), site, args.window)

        counts.append(match_count)

    df1[new_column] = counts
    df1.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
