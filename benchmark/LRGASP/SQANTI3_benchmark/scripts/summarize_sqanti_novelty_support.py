#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

NOVEL_CATEGORIES = {"novel_in_catalog", "novel_not_in_catalog"}
CANONICAL_SPLICE_SITES = {"GTAG", "GCAG", "ATAC"}


# =========================================================
# Helpers
# =========================================================

def infer_annotator_from_path(path: str) -> str:
    p = str(path)
    patterns = {
        "LoRTIA": r"(^|[\\/])LoRTIA([\\/_-]|$)",
        "BAMBU": r"(^|[\\/])BAMBU([\\/_-]|$)",
        "FLAIR": r"(^|[\\/])FLAIR([\\/_-]|$)",
        "IsoQuant": r"(^|[\\/])IsoQuant([\\/_-]|$)",
        "NAGATA": r"(^|[\\/])NAGATA([\\/_-]|$)",
    }
    for name, pat in patterns.items():
        if re.search(pat, p, flags=re.IGNORECASE):
            return name
    raise ValueError(f"Could not infer annotator from path: {path}")


def parse_bool_value(x):
    if pd.isna(x):
        return pd.NA
    v = str(x).strip().lower()
    if v in {"true", "t", "1", "yes", "y"}:
        return True
    if v in {"false", "f", "0", "no", "n", "na", "nan", "none", "", "????"}:
        return False
    return pd.NA


def parse_bool_series(s: pd.Series) -> pd.Series:
    return s.map(parse_bool_value).astype("boolean")


def normalize_structural_category(s: pd.Series) -> pd.Series:
    return (
        s.astype(str)
        .str.strip()
        .str.lower()
        .str.replace("-", "_", regex=False)
        .str.replace(" ", "_", regex=False)
    )


def normalize_all_canonical_state(s: pd.Series) -> pd.Series:
    v = s.astype(str).str.strip().str.lower()
    out = np.select(
        [
            v.eq("canonical"),
            v.eq("non_canonical"),
            v.eq("na") | v.eq("nan") | v.eq("") | v.eq("none"),
        ],
        [
            "canonical",
            "non_canonical",
            "NA",
        ],
        default="NA",
    )
    return pd.Series(out, index=s.index, dtype="object")


def classify_support_df(df: pd.DataFrame) -> pd.DataFrame:
    cond = [
        df["annotator_count"] >= 2,
        df["chemistry_count"] >= 2,
        df["cell_line_count"] >= 2,
    ]
    choices = [
        "cross_tool",
        "multi_chemistry_single_tool",
        "multi_cell_line_single_tool",
    ]
    df["support_class"] = np.select(cond, choices, default="single_condition")
    df["independent_support"] = df["support_class"] != "single_condition"
    return df


# =========================================================
# Loaders
# =========================================================

def load_classification_database_fast(manifest_path: str) -> pd.DataFrame:
    manifest = pd.read_csv(manifest_path, sep="\t").rename(
        columns={
            "Chemistry": "chemistry",
            "Cell-line": "cell_line",
            "GTF": "classification_path",
            "Annotator": "annotator",
        }
    )

    required = {"chemistry", "cell_line", "classification_path"}
    missing = required - set(manifest.columns)
    if missing:
        raise ValueError(f"Classification manifest missing columns: {sorted(missing)}")

    usecols = [
        "isoform",
        "structural_category",
        "associated_gene",
        "associated_transcript",
        "length",
        "exons",
        "all_canonical",
        "predicted_NMD",
        "bite",
    ]

    dfs = []
    for rec in manifest.to_dict("records"):
        path = rec["classification_path"]
        annotator = rec.get("annotator")
        if annotator is None or pd.isna(annotator):
            annotator = infer_annotator_from_path(path)

        df = pd.read_csv(path, sep="\t", usecols=lambda c: c in usecols)
        df["annotator"] = annotator
        df["chemistry"] = rec["chemistry"]
        df["cell_line"] = rec["cell_line"]
        dfs.append(df)

    out = pd.concat(dfs, ignore_index=True)

    out["structural_category_norm"] = normalize_structural_category(out["structural_category"])
    out["is_nic"] = out["structural_category_norm"].eq("novel_in_catalog")
    out["is_nnc"] = out["structural_category_norm"].eq("novel_not_in_catalog")
    out["is_novel_isoform"] = out["structural_category_norm"].isin(NOVEL_CATEGORIES)

    out["length_num"] = pd.to_numeric(out["length"], errors="coerce")
    out["exons_num"] = pd.to_numeric(out["exons"], errors="coerce")

    if "all_canonical" in out.columns:
        out["all_canonical_state"] = normalize_all_canonical_state(out["all_canonical"])
    else:
        out["all_canonical_state"] = "NA"

    if "predicted_NMD" in out.columns:
        out["predicted_NMD_bool"] = parse_bool_series(out["predicted_NMD"])
    else:
        out["predicted_NMD_bool"] = pd.Series(pd.NA, index=out.index, dtype="boolean")

    if "bite" in out.columns:
        out["bite_bool"] = parse_bool_series(out["bite"])
    else:
        out["bite_bool"] = pd.Series(pd.NA, index=out.index, dtype="boolean")

    for col in ["annotator", "chemistry", "cell_line", "structural_category_norm", "associated_gene"]:
        if col in out.columns:
            out[col] = out[col].astype("category")

    return out


def load_junction_database_fast(manifest_path: str) -> pd.DataFrame:
    manifest = pd.read_csv(manifest_path, sep="\t").rename(
        columns={
            "Annotator": "annotator",
            "Chemical": "chemistry",
            "Cell": "cell_line",
            "Junctions": "junctions_path",
        }
    )

    required = {"annotator", "chemistry", "cell_line", "junctions_path"}
    missing = required - set(manifest.columns)
    if missing:
        raise ValueError(f"Junction manifest missing columns: {sorted(missing)}")

    usecols = [
        "isoform",
        "chrom",
        "strand",
        "genomic_start_coord",
        "genomic_end_coord",
        "junction_category",
        "splice_site",
        "RTS_junction",
        "indel_near_junct",
    ]

    dfs = []
    for rec in manifest.to_dict("records"):
        df = pd.read_csv(rec["junctions_path"], sep="\t", usecols=lambda c: c in usecols)
        df["annotator"] = rec["annotator"]
        df["chemistry"] = rec["chemistry"]
        df["cell_line"] = rec["cell_line"]
        dfs.append(df)

    out = pd.concat(dfs, ignore_index=True)

    out["junction_category_norm"] = (
        out["junction_category"].astype(str).str.strip().str.lower().str.replace("-", "_", regex=False)
    )
    out["is_novel_junction"] = out["junction_category_norm"].ne("known")

    out["splice_site_upper"] = out["splice_site"].astype(str).str.upper().str.strip()
    out["canonical_from_splice_site"] = out["splice_site_upper"].isin(CANONICAL_SPLICE_SITES)

    if "RTS_junction" in out.columns:
        out["RTS_junction_bool"] = parse_bool_series(out["RTS_junction"])
    else:
        out["RTS_junction_bool"] = pd.Series(pd.NA, index=out.index, dtype="boolean")

    if "indel_near_junct" in out.columns:
        out["indel_near_junct_bool"] = parse_bool_series(out["indel_near_junct"])
    else:
        out["indel_near_junct_bool"] = pd.Series(pd.NA, index=out.index, dtype="boolean")

    out["junction_key"] = (
        out["chrom"].astype(str)
        + "|"
        + out["strand"].astype(str)
        + "|"
        + out["genomic_start_coord"].astype(str)
        + "|"
        + out["genomic_end_coord"].astype(str)
    )

    for col in ["annotator", "chemistry", "cell_line", "junction_category_norm", "splice_site_upper"]:
        out[col] = out[col].astype("category")

    return out


def merge_minimal(class_df: pd.DataFrame, junc_df: pd.DataFrame) -> pd.DataFrame:
    keep_cols = [
        "isoform",
        "annotator",
        "chemistry",
        "cell_line",
        "associated_gene",
        "structural_category_norm",
        "is_nic",
        "is_nnc",
        "is_novel_isoform",
    ]
    meta = class_df[keep_cols].drop_duplicates()

    merged = junc_df.merge(
        meta,
        on=["isoform", "annotator", "chemistry", "cell_line"],
        how="left",
    )
    return merged


# =========================================================
# Summaries
# =========================================================

def summarize_novelty_burden(class_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    by_condition = class_df.groupby(
        ["annotator", "chemistry", "cell_line"], observed=True, sort=False
    ).agg(
        n_total_isoforms=("isoform", "size"),
        n_nic=("is_nic", "sum"),
        n_nnc=("is_nnc", "sum"),
        n_novel_isoforms=("is_novel_isoform", "sum"),
    ).reset_index()

    by_condition["frac_nic"] = by_condition["n_nic"] / by_condition["n_total_isoforms"]
    by_condition["frac_nnc"] = by_condition["n_nnc"] / by_condition["n_total_isoforms"]
    by_condition["frac_novel_isoforms"] = by_condition["n_novel_isoforms"] / by_condition["n_total_isoforms"]
    by_condition["nnc_to_nic_ratio"] = by_condition["n_nnc"] / by_condition["n_nic"].replace(0, np.nan)

    by_annotator = class_df.groupby(
        ["annotator"], observed=True, sort=False
    ).agg(
        n_total_isoforms=("isoform", "size"),
        n_nic=("is_nic", "sum"),
        n_nnc=("is_nnc", "sum"),
        n_novel_isoforms=("is_novel_isoform", "sum"),
    ).reset_index()

    by_annotator["frac_nic"] = by_annotator["n_nic"] / by_annotator["n_total_isoforms"]
    by_annotator["frac_nnc"] = by_annotator["n_nnc"] / by_annotator["n_total_isoforms"]
    by_annotator["frac_novel_isoforms"] = by_annotator["n_novel_isoforms"] / by_annotator["n_total_isoforms"]
    by_annotator["nnc_to_nic_ratio"] = by_annotator["n_nnc"] / by_annotator["n_nic"].replace(0, np.nan)

    return by_condition, by_annotator


def summarize_junction_support(merged_junc: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    novel = merged_junc.loc[
        merged_junc["is_novel_junction"],
        ["junction_key", "annotator", "chemistry", "cell_line"]
    ].drop_duplicates()

    support = novel.groupby("junction_key", observed=True, sort=False).agg(
        annotator_count=("annotator", "nunique"),
        chemistry_count=("chemistry", "nunique"),
        cell_line_count=("cell_line", "nunique"),
    ).reset_index()

    support = classify_support_df(support)

    annotator_junc = novel[["annotator", "junction_key"]].drop_duplicates().merge(
        support[["junction_key", "support_class", "independent_support"]],
        on="junction_key",
        how="left",
    )

    by_annotator = annotator_junc.groupby(
        "annotator", observed=True, sort=False
    ).agg(
        n_unique_novel_junctions=("junction_key", "size"),
        frac_cross_tool=("support_class", lambda s: (s == "cross_tool").mean()),
        frac_multi_chemistry_single_tool=("support_class", lambda s: (s == "multi_chemistry_single_tool").mean()),
        frac_multi_cell_line_single_tool=("support_class", lambda s: (s == "multi_cell_line_single_tool").mean()),
        frac_single_condition=("support_class", lambda s: (s == "single_condition").mean()),
        frac_independent_support=("independent_support", "mean"),
    ).reset_index()

    annotator_chem_junc = novel[["annotator", "chemistry", "junction_key"]].drop_duplicates().merge(
        support[["junction_key", "support_class", "independent_support"]],
        on="junction_key",
        how="left",
    )

    by_annotator_chem = annotator_chem_junc.groupby(
        ["annotator", "chemistry"], observed=True, sort=False
    ).agg(
        n_unique_novel_junctions=("junction_key", "size"),
        frac_cross_tool=("support_class", lambda s: (s == "cross_tool").mean()),
        frac_multi_chemistry_single_tool=("support_class", lambda s: (s == "multi_chemistry_single_tool").mean()),
        frac_multi_cell_line_single_tool=("support_class", lambda s: (s == "multi_cell_line_single_tool").mean()),
        frac_single_condition=("support_class", lambda s: (s == "single_condition").mean()),
        frac_independent_support=("independent_support", "mean"),
    ).reset_index()

    return support, by_annotator, by_annotator_chem


def summarize_junction_plausibility(merged_junc: pd.DataFrame) -> pd.DataFrame:
    novel = merged_junc.loc[merged_junc["is_novel_junction"]].copy()
    novel["scope"] = "all_novel"

    nnc_only = novel.loc[novel["is_nnc"].fillna(False)].copy()
    nnc_only["scope"] = "nnc_only"

    both = pd.concat([novel, nnc_only], ignore_index=True)

    out = both.groupby(
        ["annotator", "chemistry", "scope"], observed=True, sort=False
    ).agg(
        n_junction_rows=("junction_key", "size"),
        canonical_fraction=("canonical_from_splice_site", "mean"),
        noncanonical_fraction=("canonical_from_splice_site", lambda s: (~s.fillna(False)).mean()),
        gtag_fraction=("splice_site_upper", lambda s: (s == "GTAG").mean()),
        gcag_fraction=("splice_site_upper", lambda s: (s == "GCAG").mean()),
        atac_fraction=("splice_site_upper", lambda s: (s == "ATAC").mean()),
        rts_fraction=("RTS_junction_bool", "mean"),
        indel_near_fraction=("indel_near_junct_bool", "mean"),
    ).reset_index()

    return out


def summarize_isoform_canonicality(class_df: pd.DataFrame) -> pd.DataFrame:
    novel = class_df.loc[class_df["is_novel_isoform"]].copy()
    novel["category"] = np.where(
        novel["is_nic"], "NIC",
        np.where(novel["is_nnc"], "NNC", "OTHER")
    )

    out = novel.groupby(
        ["annotator", "chemistry", "category"], observed=True, sort=False
    ).agg(
        n_isoforms=("isoform", "size"),
        canonical_isoform_fraction=("all_canonical_state", lambda s: (s == "canonical").mean()),
        noncanonical_isoform_fraction=("all_canonical_state", lambda s: (s == "non_canonical").mean()),
        na_all_canonical_fraction=("all_canonical_state", lambda s: (s == "NA").mean()),
        predicted_nmd_fraction=("predicted_NMD_bool", "mean"),
        bite_fraction=("bite_bool", "mean"),
        median_length=("length_num", "median"),
        median_exons=("exons_num", "median"),
    ).reset_index()

    return out


def summarize_locus_consistency(
    class_df: pd.DataFrame,
    merged_junc: pd.DataFrame,
    support_df: pd.DataFrame,
) -> pd.DataFrame:
    novel_iso = class_df.loc[
        class_df["is_novel_isoform"] & class_df["associated_gene"].notna(),
        ["isoform", "annotator", "chemistry", "cell_line", "associated_gene", "is_nic", "is_nnc"]
    ].copy()

    novel_iso["category"] = np.where(
        novel_iso["is_nic"], "NIC",
        np.where(novel_iso["is_nnc"], "NNC", "OTHER")
    )

    iso_gene_counts = novel_iso.groupby(
        ["annotator", "chemistry", "category", "associated_gene"],
        observed=True,
        sort=False
    ).size().reset_index(name="n_novel_isoforms")

    locus_base = iso_gene_counts.groupby(
        ["annotator", "chemistry", "category"],
        observed=True,
        sort=False
    ).agg(
        genes_with_novelty=("associated_gene", "nunique"),
        median_novel_isoforms_per_gene=("n_novel_isoforms", "median"),
        genes_with_ge_2_novel_isoforms=("n_novel_isoforms", lambda s: (s >= 2).sum()),
        genes_with_ge_3_novel_isoforms=("n_novel_isoforms", lambda s: (s >= 3).sum()),
        frac_genes_with_ge_2_novel_isoforms=("n_novel_isoforms", lambda s: (s >= 2).mean()),
    ).reset_index()

    novel_j = merged_junc.loc[
        merged_junc["is_novel_junction"] & merged_junc["associated_gene"].notna(),
        ["junction_key", "annotator", "chemistry", "associated_gene", "is_nic", "is_nnc"]
    ].drop_duplicates().merge(
        support_df[["junction_key", "support_class", "independent_support"]],
        on="junction_key",
        how="left",
    )

    novel_j["category"] = np.where(
        novel_j["is_nic"].fillna(False), "NIC",
        np.where(novel_j["is_nnc"].fillna(False), "NNC", "OTHER")
    )

    cross_tool_genes = novel_j.loc[novel_j["support_class"] == "cross_tool"].groupby(
        ["annotator", "chemistry", "category"],
        observed=True,
        sort=False
    )["associated_gene"].nunique().reset_index(name="genes_with_cross_tool_supported_novel_junction")

    ind_genes = novel_j.loc[novel_j["independent_support"].fillna(False)].groupby(
        ["annotator", "chemistry", "category"],
        observed=True,
        sort=False
    )["associated_gene"].nunique().reset_index(name="genes_with_independent_supported_novel_junction")

    out = locus_base.merge(
        cross_tool_genes,
        on=["annotator", "chemistry", "category"],
        how="left"
    ).merge(
        ind_genes,
        on=["annotator", "chemistry", "category"],
        how="left"
    )

    out["genes_with_cross_tool_supported_novel_junction"] = (
        out["genes_with_cross_tool_supported_novel_junction"].fillna(0).astype(int)
    )
    out["genes_with_independent_supported_novel_junction"] = (
        out["genes_with_independent_supported_novel_junction"].fillna(0).astype(int)
    )

    return out


# =========================================================
# Main
# =========================================================

def main():
    ap = argparse.ArgumentParser(
        description="Fast NIC/NNC support and plausibility summaries from SQANTI classification and junction databases."
    )
    ap.add_argument(
        "--classification-manifest",
        required=True,
        help="TSV with columns Chemistry, Cell-line, GTF[, Annotator]"
    )
    ap.add_argument(
        "--junction-manifest",
        required=True,
        help="TSV with columns Annotator, Chemical, Cell, Junctions"
    )
    ap.add_argument("--outdir", required=True, help="Output directory")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("Loading classification database...")
    class_df = load_classification_database_fast(args.classification_manifest)

    print("Loading junction database...")
    junc_df = load_junction_database_fast(args.junction_manifest)

    print("Merging minimal metadata...")
    merged_junc = merge_minimal(class_df, junc_df)

    print("Summarizing novelty burden...")
    burden_by_condition, burden_by_annotator = summarize_novelty_burden(class_df)
    burden_by_condition.to_csv(outdir / "nic_nnc_burden_by_condition.tsv", sep="\t", index=False)
    burden_by_annotator.to_csv(outdir / "nic_nnc_burden_by_annotator.tsv", sep="\t", index=False)

    print("Summarizing novel junction support...")
    support_df, support_by_annotator, support_by_annotator_chem = summarize_junction_support(merged_junc)
    support_by_annotator.to_csv(outdir / "novel_junction_support_by_annotator.tsv", sep="\t", index=False)
    support_by_annotator_chem.to_csv(outdir / "novel_junction_support_by_annotator_chemistry.tsv", sep="\t", index=False)

    print("Summarizing novel junction plausibility...")
    plausibility_df = summarize_junction_plausibility(merged_junc)
    plausibility_df.to_csv(outdir / "novel_junction_plausibility_by_annotator_chemistry.tsv", sep="\t", index=False)

    print("Summarizing novel isoform canonicality...")
    isoform_canonicality_df = summarize_isoform_canonicality(class_df)
    isoform_canonicality_df.to_csv(
        outdir / "novel_isoform_canonicality_by_annotator_chemistry_category.tsv",
        sep="\t",
        index=False,
    )

    print("Summarizing locus-level consistency...")
    locus_df = summarize_locus_consistency(class_df, merged_junc, support_df)
    locus_df.to_csv(
        outdir / "novel_locus_consistency_by_annotator_chemistry_category.tsv",
        sep="\t",
        index=False,
    )

    print(f"Done. Results written to: {outdir}")


if __name__ == "__main__":
    main()