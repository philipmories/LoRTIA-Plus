#!/usr/bin/env python3
"""Build an auditable six-annotator SQANTI3 basic-statistics manifest."""

import argparse
import os
import re
import sys

import pandas as pd


ANNOTATOR_PATTERNS = [
    (re.compile(r"stringtie2", re.I), "StringTie2"),
    (re.compile(r"(?:^|[/_.-])lortia(?:[/_.-]|$)", re.I), "LoRTIA-Plus"),
    (re.compile(r"(?:^|[/_.-])isoquant(?:[/_.-]|$)", re.I), "IsoQuant"),
    (re.compile(r"(?:^|[/_.-])flair(?:[/_.-]|$)", re.I), "FLAIR"),
    (re.compile(r"(?:^|[/_.-])nagata(?:[/_.-]|$)", re.I), "NAGATA"),
    (re.compile(r"(?:^|[/_.-])bambu(?:[/_.-]|$)", re.I), "bambu"),
]

EXPECTED_ANNOTATORS = [
    "LoRTIA-Plus", "bambu", "FLAIR", "IsoQuant", "NAGATA", "StringTie2"
]
EXPECTED_CHEMISTRIES = [
    "ONT-CapTrap", "ONT-cDNA", "ONT-dRNA", "PacBio", "PacBio-CapTrap"
]
EXPECTED_CELLS = ["H1", "H1-endo", "WTC11"]


def infer_annotator(path: str) -> str:
    for pattern, name in ANNOTATOR_PATTERNS:
        if pattern.search(path):
            return name
    raise ValueError(f"cannot infer annotator from path: {path}")


def normalize_cell(value: str) -> str:
    """Normalize known filename-derived cell labels without changing biology."""
    cell = value.strip()
    low = cell.lower()
    if low.startswith("h1-endo"):
        return "H1-endo"
    if low.startswith("wtc11"):
        return "WTC11"
    if low.startswith("h1"):
        return "H1"
    return cell


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Convert the transcript-recovery manifest to the format required "
            "by the SQANTI3 basic-statistics workflow."
        )
    )
    parser.add_argument("-i", "--input", required=True, help="Recovery manifest TSV.")
    parser.add_argument("-o", "--output", required=True, help="Output basic-stats manifest TSV.")
    parser.add_argument(
        "--no-check-files", action="store_true",
        help="Do not verify that every classification file exists."
    )
    parser.add_argument(
        "--allow-incomplete", action="store_true",
        help="Allow a design other than 6 annotators x 5 chemistries x 3 cell lines."
    )
    args = parser.parse_args()

    source = pd.read_csv(args.input, sep="\t", dtype=str)
    required = ["Chemistry", "Cell-line", "GTF"]
    missing = [column for column in required if column not in source.columns]
    if missing:
        sys.exit("[ERROR] Missing source column(s): " + ", ".join(missing))

    rows = []
    errors = []
    cell_normalizations = []
    for index, row in source.iterrows():
        path = str(row["GTF"]).strip()
        try:
            annotator = infer_annotator(path)
        except ValueError as exc:
            errors.append(f"row {index + 2}: {exc}")
            continue
        if not args.no_check_files and not os.path.isfile(path):
            errors.append(f"row {index + 2}: missing classification file: {path}")
        raw_cell = str(row["Cell-line"]).strip()
        cell = normalize_cell(raw_cell)
        if cell != raw_cell:
            cell_normalizations.append((index + 2, raw_cell, cell))
        rows.append({
            "Annotator": annotator,
            "Chemical": str(row["Chemistry"]).strip(),
            "Cell": cell,
            "Classification": path,
        })

    if errors:
        sys.exit("[ERROR]\n" + "\n".join(errors[:25]))

    output = pd.DataFrame(rows).drop_duplicates()
    duplicate_keys = output.duplicated(["Annotator", "Chemical", "Cell"], keep=False)
    if duplicate_keys.any():
        duplicated = output.loc[duplicate_keys, ["Annotator", "Chemical", "Cell"]]
        sys.exit("[ERROR] Duplicate design cells:\n" + duplicated.to_string(index=False))

    if not args.allow_incomplete:
        observed = set(
            map(tuple, output[["Annotator", "Chemical", "Cell"]].to_records(index=False))
        )
        expected = {
            (annotator, chemistry, cell)
            for annotator in EXPECTED_ANNOTATORS
            for chemistry in EXPECTED_CHEMISTRIES
            for cell in EXPECTED_CELLS
        }
        missing_design = sorted(expected - observed)
        extra_design = sorted(observed - expected)
        if missing_design or extra_design:
            details = []
            if missing_design:
                details.append("Missing design cells:\n" + "\n".join(map(str, missing_design[:25])))
            if extra_design:
                details.append("Unexpected design cells:\n" + "\n".join(map(str, extra_design[:25])))
            sys.exit("[ERROR] Incomplete or unexpected design.\n" + "\n".join(details))

    annotator_order = {name: i for i, name in enumerate(EXPECTED_ANNOTATORS)}
    chemistry_order = {name: i for i, name in enumerate(EXPECTED_CHEMISTRIES)}
    cell_order = {name: i for i, name in enumerate(EXPECTED_CELLS)}
    output["_a"] = output["Annotator"].map(annotator_order)
    output["_h"] = output["Chemical"].map(chemistry_order)
    output["_c"] = output["Cell"].map(cell_order)
    output = output.sort_values(["_h", "_c", "_a"]).drop(columns=["_a", "_h", "_c"])

    os.makedirs(os.path.dirname(os.path.abspath(args.output)) or ".", exist_ok=True)
    output.to_csv(args.output, sep="\t", index=False)

    counts = output.groupby("Annotator").size().to_dict()
    for row_number, before, after in cell_normalizations:
        print(f"[WARN] normalized Cell-line at row {row_number}: {before} -> {after}")
    print(f"[OK] wrote {args.output}")
    print(
        f"[OK] rows={len(output)} "
        f"groups={output[['Chemical', 'Cell']].drop_duplicates().shape[0]}"
    )
    for annotator in EXPECTED_ANNOTATORS:
        print(f"[OK] {annotator}: {counts.get(annotator, 0)}")


if __name__ == "__main__":
    main()
