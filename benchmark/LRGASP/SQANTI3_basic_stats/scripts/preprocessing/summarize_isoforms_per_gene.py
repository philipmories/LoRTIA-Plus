#!/usr/bin/env python3
"""
Compute gene-level isoform complexity statistics AND SQANTI3 structural category
statistics from SQANTI3 classification tables, based on a database file
(database_with_TES.txt) that lists all runs.

Input database format (tab-separated), columns (minimally):
    Annotator   Chemical   Cell   Classification   ...

For each row, the script:
  - reads the SQANTI3 classification file (Classification column),
  - extracts isoform, associated_gene, structural_category,
  - computes per Annotator × Chemical × Cell:

    1) Gene-level isoform complexity stats:
       - number of isoforms per gene
       - bins: 1, 2-3, 4-5, >=6 isoforms per gene
       - per-bin counts and fractions
       - mean / median isoforms per gene

    2) Structural category stats:
       - count of isoforms per SQANTI3 structural_category
       - fraction of isoforms per category

Outputs:
  --output           : TSV with isoform complexity stats (for stacked bar on isoforms/gene)
  --category-output  : TSV with structural_category stats (for SQANTI3 stacked bar)

Usage:
    python sqanti_stats.py \
        --database database_with_TES.txt \
        --output gene_isoform_complexity_stats.tsv \
        --category-output sqanti_category_stats.tsv
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute isoform complexity and SQANTI3 category statistics "
                    "from SQANTI3 classification tables."
    )
    parser.add_argument(
        "-d", "--database",
        required=True,
        help="Path to database_with_TES.txt (tab-separated)."
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Path to output TSV file for isoform complexity statistics."
    )
    parser.add_argument(
        "--category-output",
        required=True,
        help="Path to output TSV file for SQANTI3 structural category statistics."
    )
    parser.add_argument(
        "--skip-missing",
        action="store_true",
        help="Skip rows whose Classification file does not exist instead of failing."
    )
    return parser.parse_args()


def bin_isoform_count(n: int) -> str:
    """
    Bin isoform counts into ASCII-friendly categories:
        1, 2-3, 4-5, >=6
    """
    if n <= 1:
        return "1"
    elif 2 <= n <= 3:
        return "2-3"
    elif 4 <= n <= 5:
        return "4-5"
    else:
        return ">=6"


def main():
    args = parse_args()

    db_path = Path(args.database)
    if not db_path.is_file():
        sys.stderr.write(f"[ERROR] Database file not found: {db_path}\n")
        sys.exit(1)

    # Read database_with_TES.txt
    try:
        db = pd.read_csv(db_path, sep="\t", dtype=str)
    except Exception as e:
        sys.stderr.write(f"[ERROR] Failed to read database file {db_path}: {e}\n")
        sys.exit(1)

    required_cols = ["Annotator", "Chemical", "Cell", "Classification"]
    for col in required_cols:
        if col not in db.columns:
            sys.stderr.write(f"[ERROR] Required column '{col}' missing from database file.\n")
            sys.exit(1)

    all_isoforms = []

    for _, row in db.iterrows():
        annotator = row["Annotator"]
        chemical = row["Chemical"]
        cell = row["Cell"]
        class_path = Path(row["Classification"])

        if not class_path.is_file():
            msg = f"[WARNING] Classification file not found, skipping: {class_path}\n"
            if args.skip_missing:
                sys.stderr.write(msg)
                continue
            else:
                sys.stderr.write(msg.replace("[WARNING]", "[ERROR]"))
                sys.exit(1)

        sys.stderr.write(
            f"[INFO] Reading classification: {class_path} "
            f"({annotator}, {chemical}, {cell})\n"
        )

        try:
            cls = pd.read_csv(class_path, sep="\t", dtype=str)
        except Exception as e:
            sys.stderr.write(f"[ERROR] Failed to read classification file {class_path}: {e}\n")
            sys.exit(1)

        for col in ["isoform", "associated_gene", "structural_category"]:
            if col not in cls.columns:
                sys.stderr.write(
                    f"[ERROR] Column '{col}' not found in classification file {class_path}\n"
                )
                sys.exit(1)

        cls = cls[["isoform", "associated_gene", "structural_category"]].copy()

        # Drop isoforms without a valid gene ID
        cls = cls[
            cls["associated_gene"].notna()
            & (cls["associated_gene"] != ".")
            & (cls["associated_gene"] != "")
        ]

        if cls.empty:
            sys.stderr.write(
                f"[WARNING] No isoforms with valid associated_gene in {class_path}, skipping.\n"
            )
            continue

        cls["Annotator"] = annotator
        cls["Chemical"] = chemical
        cls["Cell"] = cell

        all_isoforms.append(cls)

    if not all_isoforms:
        sys.stderr.write("[ERROR] No isoform data collected from any classification file.\n")
        sys.exit(1)

    iso_df = pd.concat(all_isoforms, ignore_index=True)

    # ----------------------------------------------------------------------
    # 1) Gene-level isoform complexity stats (isoform_bin)
    # ----------------------------------------------------------------------
    gene_counts = (
        iso_df
        .groupby(["Annotator", "Chemical", "Cell", "associated_gene"], as_index=False)["isoform"]
        .nunique()
        .rename(columns={"isoform": "n_isoforms", "associated_gene": "gene_id"})
    )

    gene_counts["n_isoforms"] = gene_counts["n_isoforms"].astype(int)
    gene_counts["isoform_bin"] = gene_counts["n_isoforms"].apply(bin_isoform_count)
    gene_counts["isoform_bin"] = gene_counts["isoform_bin"].astype(str).str.strip()

    bin_order = ["1", "2-3", "4-5", ">=6"]

    stats_list = []
    group_cols = ["Annotator", "Chemical", "Cell"]

    for (annotator, chemical, cell), df_sub in gene_counts.groupby(group_cols):
        if df_sub.empty:
            continue

        n_genes_total = df_sub["gene_id"].nunique()
        mean_iso = df_sub["n_isoforms"].mean()
        median_iso = df_sub["n_isoforms"].median()

        bin_counts = (
            df_sub
            .groupby("isoform_bin")["gene_id"]
            .nunique()
            .reindex(bin_order, fill_value=0)
        )

        # ellenőrzés: a bin-ek összege egyezzen az összes génnel
        if bin_counts.sum() != n_genes_total:
            sys.stderr.write(
                f"[WARNING] Bin counts ({bin_counts.sum()}) != n_genes_total "
                f"({n_genes_total}) for {annotator}, {chemical}, {cell}\n"
            )

        for iso_bin, n_bin in bin_counts.items():
            frac_bin = n_bin / n_genes_total if n_genes_total > 0 else 0.0
            stats_list.append({
                "Annotator": annotator,
                "Chemical": chemical,
                "Cell": cell,
                "n_genes_total": n_genes_total,
                "isoform_bin": iso_bin,
                "n_genes_in_bin": int(n_bin),
                "frac_genes_in_bin": frac_bin,
                "mean_isoforms_per_gene": round(mean_iso, 3),
                "median_isoforms_per_gene": round(median_iso, 3),
            })

    iso_stats_df = pd.DataFrame(stats_list)

    out_path = Path(args.output)
    try:
        iso_stats_df.to_csv(out_path, sep="\t", index=False)
    except Exception as e:
        sys.stderr.write(f"[ERROR] Failed to write isoform complexity TSV {out_path}: {e}\n")
        sys.exit(1)

    sys.stderr.write(
        f"[INFO] Wrote isoform complexity statistics to {out_path}\n"
        f"[INFO] Rows: {len(iso_stats_df)}\n"
    )

    # ----------------------------------------------------------------------
    # 2) Structural category stats (SQANTI3 structural_category)
    # ----------------------------------------------------------------------
    cat_stats_list = []

    for (annotator, chemical, cell), df_sub in iso_df.groupby(group_cols):
        if df_sub.empty:
            continue

        n_iso_total = df_sub["isoform"].nunique()

        cat_counts = (
            df_sub
            .groupby("structural_category")["isoform"]
            .nunique()
            .sort_values(ascending=False)
        )

        for cat, n_cat in cat_counts.items():
            frac_cat = n_cat / n_iso_total if n_iso_total > 0 else 0.0
            cat_stats_list.append({
                "Annotator": annotator,
                "Chemical": chemical,
                "Cell": cell,
                "n_isoforms_total": n_iso_total,
                "structural_category": cat,
                "n_isoforms_in_cat": int(n_cat),
                "frac_isoforms_in_cat": frac_cat,
            })

    cat_stats_df = pd.DataFrame(cat_stats_list)

    cat_out_path = Path(args.category_output)
    try:
        cat_stats_df.to_csv(cat_out_path, sep="\t", index=False)
    except Exception as e:
        sys.stderr.write(f"[ERROR] Failed to write category stats TSV {cat_out_path}: {e}\n")
        sys.exit(1)

    sys.stderr.write(
        f"[INFO] Wrote SQANTI3 structural category statistics to {cat_out_path}\n"
        f"[INFO] Rows: {len(cat_stats_df)}\n"
    )


if __name__ == "__main__":
    main()
