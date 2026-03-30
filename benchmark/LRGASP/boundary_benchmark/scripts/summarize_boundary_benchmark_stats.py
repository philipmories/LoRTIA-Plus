#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path
from typing import List, Dict, Tuple
import pandas as pd
import numpy as np
from scipy import stats
import statsmodels.api as sm
import statsmodels.formula.api as smf

DEFAULT_LORTIA = "LoRTIA"

def bh_fdr(pvals: List[float]) -> List[float]:
    """Benjamini–Hochberg FDR adjusted p-values."""
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    if n == 0:
        return []
    order = np.argsort(pvals)
    ranks = np.arange(1, n + 1)
    adj = np.empty(n, dtype=float)
    adj[order] = pvals[order] * n / ranks
    # enforce monotonicity
    adj_sorted = adj[order]
    adj_sorted = np.minimum.accumulate(adj_sorted[::-1])[::-1]
    adj[order] = np.clip(adj_sorted, 0.0, 1.0)
    return adj.tolist()

def holm(pvals: List[float]) -> List[float]:
    """Holm step-down adjusted p-values."""
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    if n == 0:
        return []
    order = np.argsort(pvals)
    adj = np.empty(n, dtype=float)
    for i, idx in enumerate(order):
        adj[idx] = (n - i) * pvals[idx]
    # step-down monotonicity
    ordered_adj = adj[order]
    ordered_adj = np.maximum.accumulate(ordered_adj)
    adj[order] = np.clip(ordered_adj, 0.0, 1.0)
    return adj.tolist()

def adjust_pvals(pvals: List[float], method: str) -> List[float]:
    method = method.lower()
    if method in ("bh", "fdr", "bh-fdr", "benjamini-hochberg"):
        return bh_fdr(pvals)
    if method in ("holm", "holm-bonferroni"):
        return holm(pvals)
    if method in ("none", "raw"):
        return pvals
    raise ValueError(f"Unknown p-adjust method: {method}")

def cohen_dz_paired(diffs: np.ndarray) -> float:
    """Cohen's dz for paired designs: mean(diff)/sd(diff)."""
    diffs = np.asarray(diffs, dtype=float)
    sd = diffs.std(ddof=1) if len(diffs) > 1 else np.nan
    if sd == 0 or np.isnan(sd):
        return np.nan
    return diffs.mean() / sd

def read_inputs(paths: List[Path]) -> pd.DataFrame:
    dfs = []
    for p in paths:
        df = pd.read_csv(p, sep="\t")
        df["__source_file__"] = str(p)
        dfs.append(df)
    out = pd.concat(dfs, ignore_index=True)
    # normalize types
    if "window" in out.columns:
        out["window"] = pd.to_numeric(out["window"], errors="coerce").astype("Int64")
    return out

def sanity_check(df: pd.DataFrame):
    needed = {"mode", "window", "chem", "cell", "annotator", "f1"}
    missing = needed - set(df.columns)
    if missing:
        raise SystemExit(f"Missing columns in input TSV(s): {', '.join(sorted(missing))}")

def run_anova_blocked(df_sub: pd.DataFrame) -> Dict[str, float]:
    """
    Two-way additive ANOVA: F1 ~ C(annotator) + C(cell)
    Returns p-values and effect sizes for annotator and cell.
    """
    # Drop missing F1
    d = df_sub.dropna(subset=["f1"]).copy()
    # OLS with categorical factors
    model = smf.ols("f1 ~ C(annotator) + C(cell)", data=d).fit()
    aov = sm.stats.anova_lm(model, typ=2)  # Type II SS, good for additive blocked design
    # aov index names: C(annotator), C(cell), Residual
    res = {}
    for term in ["C(annotator)", "C(cell)"]:
        if term in aov.index:
            res[f"{term}_df"] = float(aov.loc[term, "df"])
            res[f"{term}_F"] = float(aov.loc[term, "F"])
            res[f"{term}_p"] = float(aov.loc[term, "PR(>F)"])
            # partial eta^2
            ss_term = float(aov.loc[term, "sum_sq"])
            ss_res = float(aov.loc["Residual", "sum_sq"]) if "Residual" in aov.index else np.nan
            res[f"{term}_partial_eta2"] = ss_term / (ss_term + ss_res) if (ss_term + ss_res) > 0 else np.nan
        else:
            res[f"{term}_df"] = np.nan
            res[f"{term}_F"] = np.nan
            res[f"{term}_p"] = np.nan
            res[f"{term}_partial_eta2"] = np.nan
    # residual df
    if "Residual" in aov.index:
        res["Residual_df"] = float(aov.loc["Residual", "df"])
        res["Residual_ss"] = float(aov.loc["Residual", "sum_sq"])
    else:
        res["Residual_df"] = np.nan
        res["Residual_ss"] = np.nan
    return res

def planned_lortia_tests(df_sub: pd.DataFrame, lortia_name: str, padj: str) -> pd.DataFrame:
    """
    For a given (mode, window, chem):
    paired t-tests across cell lines: LoRTIA vs each other annotator.
    """
    d = df_sub.copy()
    # pivot to cell x annotator
    piv = d.pivot_table(index="cell", columns="annotator", values="f1", aggfunc="mean")
    if lortia_name not in piv.columns:
        return pd.DataFrame()

    results = []
    lvals = piv[lortia_name].dropna()
    for other in piv.columns:
        if other == lortia_name:
            continue
        # matched cells present in both
        joined = pd.concat([piv[lortia_name], piv[other]], axis=1, join="inner").dropna()
        if joined.shape[0] < 2:
            continue
        diffs = (joined[lortia_name] - joined[other]).values
        t, p = stats.ttest_rel(joined[lortia_name].values, joined[other].values, nan_policy="omit")
        results.append({
            "other_annotator": other,
            "n_cells": int(joined.shape[0]),
            "mean_delta_f1": float(np.mean(diffs)),
            "sd_delta_f1": float(np.std(diffs, ddof=1)) if joined.shape[0] > 1 else np.nan,
            "cohen_dz": float(cohen_dz_paired(diffs)),
            "t": float(t),
            "p_raw": float(p),
        })

    if not results:
        return pd.DataFrame()

    res_df = pd.DataFrame(results)
    res_df["p_adj"] = adjust_pvals(res_df["p_raw"].tolist(), padj)
    res_df = res_df.sort_values(["p_adj", "p_raw", "other_annotator"]).reset_index(drop=True)
    return res_df

def cellline_pair_tests(df_sub: pd.DataFrame, padj: str) -> pd.DataFrame:
    """
    chemistry-n belül: cell line pairwise paired t-test annotátorokon párosítva (n≈5).
    Ez kiegészítő jellegű, de a régi pipeline-odhoz illeszkedik.
    """
    d = df_sub.copy()
    piv = d.pivot_table(index="annotator", columns="cell", values="f1", aggfunc="mean")

    cells = [c for c in piv.columns if pd.notna(c)]
    cells = sorted(cells)
    if len(cells) < 2:
        return pd.DataFrame()

    pairs = []
    for i in range(len(cells)):
        for j in range(i + 1, len(cells)):
            c1, c2 = cells[i], cells[j]
            joined = pd.concat([piv[c1], piv[c2]], axis=1, join="inner").dropna()
            if joined.shape[0] < 2:
                continue
            diffs = (joined[c1] - joined[c2]).values
            t, p = stats.ttest_rel(joined[c1].values, joined[c2].values, nan_policy="omit")
            pairs.append({
                "cell_a": c1,
                "cell_b": c2,
                "n_annotators": int(joined.shape[0]),
                "mean_delta_f1": float(np.mean(diffs)),
                "sd_delta_f1": float(np.std(diffs, ddof=1)) if joined.shape[0] > 1 else np.nan,
                "cohen_dz": float(cohen_dz_paired(diffs)),
                "t": float(t),
                "p_raw": float(p),
            })

    if not pairs:
        return pd.DataFrame()

    out = pd.DataFrame(pairs)
    out["p_adj"] = adjust_pvals(out["p_raw"].tolist(), padj)
    out = out.sort_values(["p_adj", "p_raw", "cell_a", "cell_b"]).reset_index(drop=True)
    return out

def means_table(df: pd.DataFrame) -> pd.DataFrame:
    return (df.groupby(["mode", "window", "chem", "cell", "annotator"], as_index=False)
              .agg(mean_f1=("f1", "mean"),
                   mean_precision=("precision", "mean"),
                   mean_recall=("recall", "mean"),
                   tp=("tp", "sum"), fp=("fp", "sum"), fn=("fn", "sum"),
                   n_pred=("n_pred", "mean"), n_ref=("n_ref", "mean")))

def main():
    ap = argparse.ArgumentParser(
        description="ANOVA + paired t-tests from benchmark TSV outputs (TSS/TES)."
    )
    ap.add_argument("-i", "--inputs", nargs="+", required=True,
                    help="Input benchmark TSV(s), e.g. TSS_03.05.tsv TES_03.05.tsv")
    ap.add_argument("-o", "--outprefix", required=True,
                    help="Output prefix (folder/prefix). Example: stats/out -> out_anova.tsv, out_lortia_vs_others.tsv, ...")
    ap.add_argument("--lortia-name", default=DEFAULT_LORTIA,
                    help='Annotator name used for planned comparisons (default: "LoRTIA")')
    ap.add_argument("--padj", choices=["bh", "holm", "none"], default="bh",
                    help="Multiple-testing correction within each chemistry for planned tests (default: bh)")
    ap.add_argument("--include-cellline-pairs", action="store_true",
                    help="Also do cell-line pairwise paired t-tests within chemistry (paired over annotators).")
    args = ap.parse_args()

    in_paths = [Path(p).resolve() for p in args.inputs]
    for p in in_paths:
        if not p.exists():
            raise SystemExit(f"Missing input: {p}")

    df = read_inputs(in_paths)
    sanity_check(df)

    # Aggregate means just in case there are duplicates
    df_use = (df.groupby(["mode","window","chem","cell","annotator"], as_index=False)
                .agg(f1=("f1","mean"),
                     precision=("precision","mean"),
                     recall=("recall","mean"),
                     tp=("tp","sum"), fp=("fp","sum"), fn=("fn","sum"),
                     n_pred=("n_pred","mean"), n_ref=("n_ref","mean")))

    # Output: means table
    means_out = Path(f"{args.outprefix}_means.tsv")
    means_table(df).to_csv(means_out, sep="\t", index=False)

    # ANOVA per (mode, window, chem)
    anova_rows = []
    lortia_rows = []
    cellpair_rows = []

    for (mode, window, chem), sub in df_use.groupby(["mode","window","chem"], sort=True):
        # ANOVA blocked by cell line
        an = run_anova_blocked(sub)
        anova_rows.append({
            "mode": mode,
            "window": int(window) if pd.notna(window) else window,
            "chem": chem,
            "n_rows": int(sub.shape[0]),
            "n_cells": int(sub["cell"].nunique()),
            "n_annotators": int(sub["annotator"].nunique()),
            "p_annotator": an["C(annotator)_p"],
            "F_annotator": an["C(annotator)_F"],
            "df_annotator": an["C(annotator)_df"],
            "partial_eta2_annotator": an["C(annotator)_partial_eta2"],
            "p_cell": an["C(cell)_p"],
            "F_cell": an["C(cell)_F"],
            "df_cell": an["C(cell)_df"],
            "partial_eta2_cell": an["C(cell)_partial_eta2"],
            "df_residual": an["Residual_df"],
        })

        # Planned: LoRTIA vs others
        post = planned_lortia_tests(sub, args.lortia_name, args.padj)
        if not post.empty:
            post.insert(0, "chem", chem)
            post.insert(0, "window", int(window) if pd.notna(window) else window)
            post.insert(0, "mode", mode)
            post.insert(3, "lortia", args.lortia_name)
            post["padj_method"] = args.padj
            lortia_rows.append(post)

        # Optional: cell-line pairwise comparisons (paired across annotators)
        if args.include_cellline_pairs:
            cp = cellline_pair_tests(sub, args.padj)
            if not cp.empty:
                cp.insert(0, "chem", chem)
                cp.insert(0, "window", int(window) if pd.notna(window) else window)
                cp.insert(0, "mode", mode)
                cp["padj_method"] = args.padj
                cellpair_rows.append(cp)

    # Write outputs
    anova_out = Path(f"{args.outprefix}_anova.tsv")
    pd.DataFrame(anova_rows).sort_values(["mode","window","chem"]).to_csv(anova_out, sep="\t", index=False)

    lortia_out = Path(f"{args.outprefix}_lortia_vs_others.tsv")
    if lortia_rows:
        pd.concat(lortia_rows, ignore_index=True).to_csv(lortia_out, sep="\t", index=False)
    else:
        pd.DataFrame().to_csv(lortia_out, sep="\t", index=False)

    if args.include_cellline_pairs:
        cell_out = Path(f"{args.outprefix}_cellline_pairs.tsv")
        if cellpair_rows:
            pd.concat(cellpair_rows, ignore_index=True).to_csv(cell_out, sep="\t", index=False)
        else:
            pd.DataFrame().to_csv(cell_out, sep="\t", index=False)

    print(f"[OK] wrote:\n  {means_out}\n  {anova_out}\n  {lortia_out}", file=sys.stderr)
    if args.include_cellline_pairs:
        print(f"  {Path(f'{args.outprefix}_cellline_pairs.tsv')}", file=sys.stderr)

if __name__ == "__main__":
    main()