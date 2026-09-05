#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
SQANTI3 classification.txt összefoglaló TSS/TES-végpont metrikákhoz (v2).

- Pontos SQANTI3 oszlopnevek:
    within_CAGE_peak
    within_polyA_site
    polyA_motif       (motif szekvencia, itt most NEM használjuk)
    polyA_dist
    polyA_motif_found (TRUE/FALSE)
"""

import argparse
import os
import sys
from typing import Dict, Any, List

import numpy as np
import pandas as pd

TRUE_SET = {"TRUE", "True", "true", "T", "1", 1, True}
FALSE_SET = {"FALSE", "False", "false", "F", "0", 0, False}

SQANTI_MAIN_CATS = {
    "full-splice_match": "FSM",
    "incomplete-splice_match": "ISM",
    "novel_in_catalog": "NIC",
    "novel_not_in_catalog": "NNC",
}

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Summarize SQANTI3 classification.txt TSS/TES metrics "
                    "by Annotator × Chemical × Cell × (FSM/ISM/NIC/NNC/ALL)."
    )
    p.add_argument("-d", "--database", required=True,
                   help="TSV like database_with_TES.txt (Annotator, Chemical, Cell, Classification).")
    p.add_argument("-o", "--output", required=True,
                   help="Output TSV.")
    p.add_argument("--sep", default="\t",
                   help="Separator for database file (default: TAB).")
    p.add_argument("--window_exact", type=int, default=5,
                   help="|diff_to_gene_TSS/TTS| <= window_exact counted as exact (default:5).")
    p.add_argument("--polyA_close_window", type=int, default=40,
                   help="polyA_dist in [-polyA_close_window,0] is 'close' (default:40).")
    return p.parse_args()


def as_bool_series(s: pd.Series) -> pd.Series:
    """TRUE/FALSE/NA -> pandas boolean + NA."""
    ss = s.astype("string")
    true_vals = {str(x) for x in TRUE_SET}
    false_vals = {str(x) for x in FALSE_SET}
    is_true = ss.isin(true_vals)
    is_false = ss.isin(false_vals)
    out = pd.Series(np.nan, index=ss.index, dtype="float64")
    out[is_true] = True
    out[is_false] = False
    return out.astype("boolean")


def safe_abs_median(series: pd.Series) -> float:
    if series is None:
        return np.nan
    s = pd.to_numeric(series, errors="coerce").dropna()
    if s.empty:
        return np.nan
    return float(np.median(np.abs(s.values)))


def exact_rate(series: pd.Series, window: int) -> float:
    if series is None:
        return np.nan
    s = pd.to_numeric(series, errors="coerce").dropna()
    if s.empty:
        return np.nan
    return float((np.abs(s.values) <= window).sum() / len(s))


def cage_metrics(df: pd.DataFrame) -> Dict[str, Any]:
    """TSS CAGE metrika: within_CAGE_peak."""
    out: Dict[str, Any] = {}
    n_iso = len(df)
    if n_iso == 0:
        out.update(
            within_cage_true=0,
            within_cage_false=0,
            within_cage_na=0,
            within_cage_rate=np.nan,
        )
        return out

    if "within_CAGE_peak" in df.columns:
        wp = as_bool_series(df["within_CAGE_peak"])
        n_true = int(wp.sum(skipna=True))
        n_false = int((~wp).sum(skipna=True))
        n_na = n_iso - n_true - n_false
        denom = n_true + n_false
        rate = n_true / denom if denom > 0 else np.nan
    else:
        n_true = n_false = n_na = 0
        rate = np.nan

    out["within_cage_true"] = n_true
    out["within_cage_false"] = n_false
    out["within_cage_na"] = n_na
    out["within_cage_rate"] = float(rate) if not np.isnan(rate) else np.nan
    return out


def polyA_metrics(df: pd.DataFrame, polyA_close_window: int) -> Dict[str, Any]:
    """
    TES polyA metrikák:

      within_polyA_site
      polyA_motif_found
      polyA_dist
    """
    out: Dict[str, Any] = {}
    n_iso = len(df)
    if n_iso == 0:
        out.update(
            polyA_site_true=0,
            polyA_site_false=0,
            polyA_site_na=0,
            polyA_site_rate=np.nan,
            polyA_motif_found_true=0,
            polyA_motif_found_false=0,
            polyA_motif_rate=np.nan,
            polyA_motif_close_rate=np.nan,
            polyA_support_true=0,
            polyA_support_rate=np.nan,
        )
        return out

    # -------- within_polyA_site --------
    if "within_polyA_site" in df.columns:
        wps = as_bool_series(df["within_polyA_site"])
        n_true = int(wps.sum(skipna=True))
        n_false = int((~wps).sum(skipna=True))
        n_na = n_iso - n_true - n_false
        denom = n_true + n_false
        rate = n_true / denom if denom > 0 else np.nan
    else:
        n_true = n_false = n_na = 0
        rate = np.nan
        wps = pd.Series([np.nan] * n_iso, index=df.index, dtype="float64")

    out["polyA_site_true"] = n_true
    out["polyA_site_false"] = n_false
    out["polyA_site_na"] = n_na
    out["polyA_site_rate"] = float(rate) if not np.isnan(rate) else np.nan

    # -------- polyA_motif_found --------
    if "polyA_motif_found" in df.columns:
        pmf = as_bool_series(df["polyA_motif_found"])
        n_motif_true = int(pmf.sum(skipna=True))
        n_motif_false = int((~pmf).sum(skipna=True))
        motif_rate = n_motif_true / (n_motif_true + n_motif_false) if (n_motif_true + n_motif_false) > 0 else np.nan
    else:
        pmf = pd.Series([np.nan] * n_iso, index=df.index, dtype="float64")
        n_motif_true = 0
        n_motif_false = 0
        motif_rate = np.nan

    out["polyA_motif_found_true"] = n_motif_true
    out["polyA_motif_found_false"] = n_motif_false
    out["polyA_motif_rate"] = float(motif_rate) if not np.isnan(motif_rate) else np.nan

    # -------- motif_close: polyA_dist in [-W, 0] & motif_found --------
    if "polyA_dist" in df.columns and "polyA_motif_found" in df.columns:
        dist = pd.to_numeric(df["polyA_dist"], errors="coerce")
        pmf_bool = as_bool_series(df["polyA_motif_found"]).fillna(False)
        motif_close_mask = pmf_bool & dist.ge(-polyA_close_window) & dist.le(0)
        n_close = int(motif_close_mask.sum())
        close_rate = n_close / n_motif_true if n_motif_true > 0 else np.nan
    else:
        motif_close_mask = pd.Series([False] * n_iso, index=df.index)
        n_close = 0
        close_rate = np.nan

    out["polyA_motif_close_rate"] = float(close_rate) if not np.isnan(close_rate) else np.nan

    # -------- polyA_support: site OR motif_close --------
    if "within_polyA_site" in df.columns:
        wps_bool = as_bool_series(df["within_polyA_site"]).fillna(False)
    else:
        wps_bool = pd.Series([False] * n_iso, index=df.index)

    support_mask = wps_bool | motif_close_mask
    n_support = int(support_mask.sum())
    support_rate = n_support / n_iso if n_iso > 0 else np.nan

    out["polyA_support_true"] = n_support
    out["polyA_support_rate"] = float(support_rate) if not np.isnan(support_rate) else np.nan

    return out


def compute_metrics_for_subset(df: pd.DataFrame,
                               category_label: str,
                               window_exact: int,
                               polyA_close_window: int) -> Dict[str, Any]:
    n_iso = len(df)
    base: Dict[str, Any] = {"category": category_label, "n_isoforms": int(n_iso)}

    if n_iso == 0:
        base.update(
            within_cage_true=0,
            within_cage_false=0,
            within_cage_na=0,
            within_cage_rate=np.nan,
            polyA_site_true=0,
            polyA_site_false=0,
            polyA_site_na=0,
            polyA_site_rate=np.nan,
            polyA_motif_found_true=0,
            polyA_motif_found_false=0,
            polyA_motif_rate=np.nan,
            polyA_motif_close_rate=np.nan,
            polyA_support_true=0,
            polyA_support_rate=np.nan,
            abs_dtss_median=np.nan,
            tss_exact_rate=np.nan,
            abs_dtts_median=np.nan,
            tes_exact_rate=np.nan,
        )
        return base

    base.update(cage_metrics(df))
    base.update(polyA_metrics(df, polyA_close_window=polyA_close_window))

    if "diff_to_gene_TSS" in df.columns:
        base["abs_dtss_median"] = safe_abs_median(df["diff_to_gene_TSS"])
        base["tss_exact_rate"] = exact_rate(df["diff_to_gene_TSS"], window_exact)
    else:
        base["abs_dtss_median"] = np.nan
        base["tss_exact_rate"] = np.nan

    if "diff_to_gene_TTS" in df.columns:
        base["abs_dtts_median"] = safe_abs_median(df["diff_to_gene_TTS"])
        base["tes_exact_rate"] = exact_rate(df["diff_to_gene_TTS"], window_exact)
    else:
        base["abs_dtts_median"] = np.nan
        base["tes_exact_rate"] = np.nan

    return base


def summarize_one_classification(classif_path: str,
                                 annotator: str,
                                 chemical: str,
                                 cell: str,
                                 window_exact: int,
                                 polyA_close_window: int) -> List[Dict[str, Any]]:
    if not os.path.exists(classif_path):
        sys.stderr.write(f"[WARN] Missing classification: {classif_path}\n")
        return []

    try:
        df = pd.read_csv(classif_path, sep="\t", low_memory=False)
    except Exception as e:
        sys.stderr.write(f"[ERROR] Failed to read {classif_path}: {e}\n")
        return []

    if "structural_category" not in df.columns:
        sys.stderr.write(f"[ERROR] structural_category missing in {classif_path}\n")
        return []

    df_main = df[df["structural_category"].isin(SQANTI_MAIN_CATS.keys())].copy()
    results: List[Dict[str, Any]] = []

    subsets = [("ALL", df_main)]
    for sq_cat, short in SQANTI_MAIN_CATS.items():
        sub = df_main[df_main["structural_category"] == sq_cat].copy()
        subsets.append((short, sub))

    for cat_label, subdf in subsets:
        metrics = compute_metrics_for_subset(
            subdf, category_label=cat_label,
            window_exact=window_exact,
            polyA_close_window=polyA_close_window,
        )
        metrics["Annotator"] = annotator
        metrics["Chemical"] = chemical
        metrics["Cell"] = cell
        metrics["source_classif"] = classif_path
        results.append(metrics)

    return results


def main():
    args = parse_args()

    try:
        db = pd.read_csv(args.database, sep=args.sep)
    except Exception as e:
        sys.stderr.write(f"[ERROR] Failed to read database {args.database}: {e}\n")
        sys.exit(1)

    required = {"Annotator", "Chemical", "Cell", "Classification"}
    missing = required - set(db.columns)
    if missing:
        sys.stderr.write(f"[ERROR] Missing columns in database: {', '.join(sorted(missing))}\n")
        sys.exit(1)

    all_rows: List[Dict[str, Any]] = []
    for _, row in db.iterrows():
        annot = str(row["Annotator"])
        chem = str(row["Chemical"])
        cell = str(row["Cell"])
        classif = str(row["Classification"])

        sys.stderr.write(f"[INFO] {annot} / {chem} / {cell}\n")

        res = summarize_one_classification(
            classif_path=classif,
            annotator=annot,
            chemical=chem,
            cell=cell,
            window_exact=args.window_exact,
            polyA_close_window=args.polyA_close_window,
        )
        all_rows.extend(res)

    if not all_rows:
        sys.stderr.write("[WARN] No metrics collected.\n")
        sys.exit(0)

    out_df = pd.DataFrame(all_rows)

    col_order = [
        "Annotator", "Chemical", "Cell", "category",
        "n_isoforms",
        "within_cage_true", "within_cage_false", "within_cage_na", "within_cage_rate",
        "polyA_site_true", "polyA_site_false", "polyA_site_na", "polyA_site_rate",
        "polyA_motif_found_true", "polyA_motif_found_false", "polyA_motif_rate",
        "polyA_motif_close_rate",
        "polyA_support_true", "polyA_support_rate",
        "abs_dtss_median", "tss_exact_rate",
        "abs_dtts_median", "tes_exact_rate",
        "source_classif",
    ]
    col_order = [c for c in col_order if c in out_df.columns]
    out_df = out_df[col_order]

    out_df.to_csv(args.output, sep="\t", index=False)
    sys.stderr.write(f"[OK] Written: {args.output}\n")


if __name__ == "__main__":
    main()
