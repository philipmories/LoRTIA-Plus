#!/usr/bin/env python3
import argparse
import csv
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
import statsmodels.formula.api as smf

# ---------------------------
# Helpers
# ---------------------------

TX_RE = re.compile(r'transcript_id "([^"]+)"')
GENE_RE = re.compile(r'gene_id "([^"]+)"')

FSM_CANON = "full-splice_match"
ISM_CANON = "incomplete-splice_match"

DEFAULT_ANNOTATORS = ["LoRTIA", "bambu", "FLAIR", "IsoQuant", "NAGATA", "BAMBU", "Flair", "Bambu"]

def open_maybe_gzip(path: Path, mode: str = "rt"):
    p = str(path)
    if p.endswith(".gz"):
        return gzip.open(p, mode, encoding="utf-8", errors="replace", newline="")
    return open(p, mode, encoding="utf-8", errors="replace", newline="")

def strip_version(tx: str) -> str:
    return re.sub(r"\.\d+$", "", tx)

def norm_annot(name: str) -> str:
    if name is None:
        return "."
    if name.upper() == "BAMBU":
        return "bambu"
    if name.lower() == "flair":
        return "FLAIR"
    if name.lower() == "lortia":
        return "LoRTIA"
    return name

def infer_annotator_from_path(p: Path) -> str:
    s = str(p)
    # try folder tokens first
    for tok in ["LoRTIA", "IsoQuant", "FLAIR", "NAGATA", "bambu", "BAMBU"]:
        if re.search(rf"(^|[/_\\-]){re.escape(tok)}($|[/_\\-])", s, flags=re.IGNORECASE):
            return norm_annot(tok)
    # fallback: filename contains
    base = p.name
    for tok in ["LoRTIA", "IsoQuant", "FLAIR", "NAGATA", "bambu", "BAMBU"]:
        if re.search(tok, base, flags=re.IGNORECASE):
            return norm_annot(tok)
    return "."

def cat_is_fsm(sc: str) -> bool:
    if sc is None:
        return False
    c = sc.strip()
    if not c:
        return False
    c_low = c.lower().replace("_", "-")
    return c == FSM_CANON or c_low == "full-splice-match"

def cat_is_ism(sc: str) -> bool:
    if sc is None:
        return False
    c = sc.strip()
    if not c:
        return False
    c_low = c.lower().replace("_", "-")
    return c == ISM_CANON or c_low == "incomplete-splice-match"

def bh_fdr(pvals: List[float]) -> List[float]:
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    if n == 0:
        return []
    order = np.argsort(p)
    ranks = np.arange(1, n + 1)
    adj = np.empty(n, dtype=float)
    adj[order] = p[order] * n / ranks
    # monotone
    adj_sorted = adj[order]
    adj_sorted = np.minimum.accumulate(adj_sorted[::-1])[::-1]
    adj[order] = np.clip(adj_sorted, 0.0, 1.0)
    return adj.tolist()

def cohen_dz_paired(diffs: np.ndarray) -> float:
    diffs = np.asarray(diffs, dtype=float)
    if diffs.size < 2:
        return np.nan
    sd = diffs.std(ddof=1)
    if sd == 0:
        return np.nan
    return float(diffs.mean() / sd)

# ---------------------------
# I/O parsing
# ---------------------------

def parse_manifest(manifest: Path) -> pd.DataFrame:
    df = pd.read_csv(manifest, sep="\t")
    required = {"Chemistry", "Cell-line", "GTF", "Reference"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"ERROR: manifest missing columns: {', '.join(sorted(missing))}")
    # Normalize whitespace / CRLF
    for c in ["Chemistry", "Cell-line", "GTF", "Reference"]:
        df[c] = df[c].astype(str).str.strip()
    return df

def read_reference_transcripts(ref_gtf: Path, strip_tx_version: bool) -> Set[str]:
    txs: Set[str] = set()
    with open_maybe_gzip(ref_gtf, "rt") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            attrs = parts[8]
            m = TX_RE.search(attrs)
            if not m:
                continue
            tx = m.group(1)
            if strip_tx_version:
                tx = strip_version(tx)
            txs.add(tx)
    return txs

def read_sqanti_sets_and_counts(classif: Path, strip_tx_version: bool):
    """
    Returns:
      fsm_set, ism_set, any_set (associated_transcript IDs),
      category_counts (structural_category isoform counts)
    """
    fsm: Set[str] = set()
    ism: Set[str] = set()
    cat_counts = defaultdict(int)

    with open_maybe_gzip(classif, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit(f"ERROR: SQANTI classification has no header: {classif}")
        if "structural_category" not in reader.fieldnames or "associated_transcript" not in reader.fieldnames:
            raise SystemExit(f"ERROR: missing required SQANTI columns in {classif} (need structural_category, associated_transcript)")

        for rec in reader:
            sc = (rec.get("structural_category") or "").strip()
            cat_counts[sc] += 1

            at = (rec.get("associated_transcript") or "").strip()
            if at.upper() in {"", "NA", "NAN"} or at == ".":
                continue
            if strip_tx_version:
                at = strip_version(at)

            if cat_is_fsm(sc):
                fsm.add(at)
            elif cat_is_ism(sc):
                ism.add(at)

    any_set = fsm | ism
    return fsm, ism, any_set, dict(cat_counts)

# ---------------------------
# Stats
# ---------------------------

def anova_blocked(df: pd.DataFrame, value_col: str) -> Dict[str, float]:
    """
    ANOVA per chemistry, blocking by cell line:
      value ~ C(annotator) + C(cell)
    """
    d = df.dropna(subset=[value_col]).copy()
    model = smf.ols(f"{value_col} ~ C(annotator) + C(cell)", data=d).fit()
    aov = sm.stats.anova_lm(model, typ=2)
    out = {}
    if "C(annotator)" in aov.index:
        out["F_annotator"] = float(aov.loc["C(annotator)", "F"])
        out["p_annotator"] = float(aov.loc["C(annotator)", "PR(>F)"])
        out["df_annotator"] = float(aov.loc["C(annotator)", "df"])
    else:
        out["F_annotator"] = np.nan
        out["p_annotator"] = np.nan
        out["df_annotator"] = np.nan

    if "C(cell)" in aov.index:
        out["F_cell"] = float(aov.loc["C(cell)", "F"])
        out["p_cell"] = float(aov.loc["C(cell)", "PR(>F)"])
        out["df_cell"] = float(aov.loc["C(cell)", "df"])
    else:
        out["F_cell"] = np.nan
        out["p_cell"] = np.nan
        out["df_cell"] = np.nan

    if "Residual" in aov.index:
        out["df_residual"] = float(aov.loc["Residual", "df"])
    else:
        out["df_residual"] = np.nan
    return out

def planned_lortia_paired(df: pd.DataFrame, value_col: str, lortia_name: str = "LoRTIA") -> pd.DataFrame:
    """
    For a chemistry: paired t-tests across cell lines comparing LoRTIA vs others.
    """
    piv = df.pivot_table(index="cell", columns="annotator", values=value_col, aggfunc="mean")
    if lortia_name not in piv.columns:
        return pd.DataFrame()

    rows = []
    for other in piv.columns:
        if other == lortia_name:
            continue
        joined = pd.concat([piv[lortia_name], piv[other]], axis=1, join="inner").dropna()
        if joined.shape[0] < 2:
            continue
        diffs = (joined[lortia_name] - joined[other]).values
        t, p = stats.ttest_rel(joined[lortia_name].values, joined[other].values, nan_policy="omit")
        rows.append({
            "other_annotator": other,
            "n_cells": int(joined.shape[0]),
            "mean_delta": float(np.mean(diffs)),
            "sd_delta": float(np.std(diffs, ddof=1)) if joined.shape[0] > 1 else np.nan,
            "cohen_dz": cohen_dz_paired(diffs),
            "t": float(t),
            "p_raw": float(p),
        })
    if not rows:
        return pd.DataFrame()

    out = pd.DataFrame(rows)
    out["p_adj_BH"] = bh_fdr(out["p_raw"].tolist())
    out = out.sort_values(["p_adj_BH", "p_raw", "other_annotator"]).reset_index(drop=True)
    return out

# ---------------------------
# Main
# ---------------------------

def main():
    ap = argparse.ArgumentParser(
        description="SQANTI3 FSM/ISM transcript-level benchmark against chemistry×cell 'cell_active' GTF references, with ANOVA + paired tests and plot-ready TSVs."
    )
    ap.add_argument("-m", "--manifest", required=True, help="TSV with columns: Chemistry, Cell-line, GTF (classification.txt), Reference (cell_active.gtf).")
    ap.add_argument("-o", "--outdir", required=True, help="Output directory.")
    ap.add_argument("--strip-tx-version", action="store_true",
                    help="Strip transcript version suffix (ENST... .2) in BOTH SQANTI and GTF matching.")
    ap.add_argument("--lortia-name", default="LoRTIA", help='Annotator label for planned comparisons (default: "LoRTIA").')
    ap.add_argument("--write-composition", action="store_true",
                    help="Also write structural_category composition counts (incl. NIC/NNC etc.) for supplementary plots.")
    ap.add_argument("--strict", action="store_true",
                    help="If any file missing -> error (default: always error for missing). If a group has empty reference -> error.")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    manifest = Path(args.manifest).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    if not manifest.exists():
        raise SystemExit(f"ERROR: manifest not found: {manifest}")

    mf = parse_manifest(manifest)

    # infer annotator from classification path
    mf["annotator"] = mf["GTF"].apply(lambda x: norm_annot(infer_annotator_from_path(Path(x))))

    # strict existence checks
    missing = []
    for _, r in mf.iterrows():
        if not Path(r["GTF"]).exists():
            missing.append(("classification", r["GTF"]))
        if not Path(r["Reference"]).exists():
            missing.append(("reference", r["Reference"]))
    if missing:
        msg = "ERROR: Missing files listed in manifest:\n" + "\n".join([f"[{k}] {p}" for k, p in missing[:50]])
        if len(missing) > 50:
            msg += f"\n... plus {len(missing)-50} more"
        raise SystemExit(msg)

    # cache reference transcript sets per (chem, cell, refpath)
    ref_cache: Dict[Tuple[str, str, str], Set[str]] = {}

    metrics_rows = []
    comp_rows = []

    for i, r in mf.iterrows():
        chem = r["Chemistry"]
        cell = r["Cell-line"]
        classif = Path(r["GTF"])
        ref = Path(r["Reference"])
        annot = r["annotator"]

        ref_key = (chem, cell, str(ref))
        if ref_key not in ref_cache:
            U = read_reference_transcripts(ref, args.strip_tx_version)
            if args.strict and len(U) == 0:
                raise SystemExit(f"ERROR: Empty reference transcript set for {chem}×{cell}: {ref}")
            ref_cache[ref_key] = U
        else:
            U = ref_cache[ref_key]

        fsm, ism, any_set, cat_counts = read_sqanti_sets_and_counts(classif, args.strip_tx_version)

        # intersect with reference universe
        fsm_in = fsm & U
        ism_in = ism & U
        any_in = any_set & U

        n_ref = len(U)
        n_fsm_in = len(fsm_in)
        n_ism_in = len(ism_in)
        n_any_in = len(any_in)

        recall_fsm = (n_fsm_in / n_ref) if n_ref else np.nan
        recall_ism = (n_ism_in / n_ref) if n_ref else np.nan
        recall_any = (n_any_in / n_ref) if n_ref else np.nan

        denom = (n_fsm_in + n_ism_in)
        fl_ratio = (n_fsm_in / denom) if denom else np.nan
        ism_share = (n_ism_in / denom) if denom else np.nan

        metrics_rows.append({
            "Chemistry": chem,
            "Cell-line": cell,
            "Annotator": annot,
            "classification_path": str(classif),
            "reference_path": str(ref),
            "n_ref_tx": n_ref,
            "n_fsm_tx": len(fsm),
            "n_ism_tx": len(ism),
            "n_any_tx": len(any_set),
            "n_fsm_in_ref": n_fsm_in,
            "n_ism_in_ref": n_ism_in,
            "n_any_in_ref": n_any_in,
            "Recall_FSM": recall_fsm,
            "Recall_ISM": recall_ism,
            "Recall_FSM_ISM": recall_any,
            "FL_ratio_FSM_over_FSM_ISM": fl_ratio,
            "ISM_share": ism_share,
        })

        if args.write_composition:
            for cat, cnt in cat_counts.items():
                comp_rows.append({
                    "Chemistry": chem,
                    "Cell-line": cell,
                    "Annotator": annot,
                    "structural_category": cat if cat else "(empty)",
                    "n_isoforms": int(cnt),
                })

        if args.verbose and (i % 10 == 0):
            print(f"[INFO] processed {i+1}/{len(mf)}: {chem} {cell} {annot}", file=sys.stderr)

    metrics = pd.DataFrame(metrics_rows)

    # ---------------------------
    # Stats (ANOVA + post hoc), per metric
    # ---------------------------
    stats_anova = []
    stats_posthoc = []

    # tidy df for stats
    d = metrics.rename(columns={"Cell-line": "cell", "Annotator": "annotator"}).copy()

    metric_cols = ["Recall_FSM_ISM", "FL_ratio_FSM_over_FSM_ISM", "Recall_FSM", "Recall_ISM", "ISM_share"]

    for metric in metric_cols:
        for chem, sub in d.groupby("Chemistry", sort=True):
            sub2 = sub[["annotator", "cell", metric]].dropna()
            if sub2.empty:
                continue

            a = anova_blocked(sub2.assign(chem=chem), metric)
            stats_anova.append({
                "Metric": metric,
                "Chemistry": chem,
                "n_rows": int(sub2.shape[0]),
                "n_cells": int(sub2["cell"].nunique()),
                "n_annotators": int(sub2["annotator"].nunique()),
                **a
            })

            post = planned_lortia_paired(sub2, metric, lortia_name=args.lortia_name)
            if not post.empty:
                post.insert(0, "Chemistry", chem)
                post.insert(0, "Metric", metric)
                post.insert(2, "Annotator_1", args.lortia_name)
                post.rename(columns={"other_annotator": "Annotator_2"}, inplace=True)
                stats_posthoc.append(post)

    anova_df = pd.DataFrame(stats_anova).sort_values(["Metric", "Chemistry"])
    posthoc_df = pd.concat(stats_posthoc, ignore_index=True) if stats_posthoc else pd.DataFrame()

    # ---------------------------
    # Plot-ready TSVs (compatible-ish with your old R)
    # ---------------------------
    # Boxplot inputs: name the value column as F1 for reuse (even though it's not F1)
    recall_plot = metrics[["Chemistry", "Cell-line", "Annotator", "Recall_FSM_ISM"]].rename(
        columns={"Cell-line": "Cell line", "Recall_FSM_ISM": "F1"}
    )
    fl_plot = metrics[["Chemistry", "Cell-line", "Annotator", "FL_ratio_FSM_over_FSM_ISM"]].rename(
        columns={"Cell-line": "Cell line", "FL_ratio_FSM_over_FSM_ISM": "F1"}
    )

    # Heatmap inputs: LoRTIA vs others from posthoc table (mean_delta + BH)
    # Use the same column names your heatmap script expects
    pairwise_recall = pd.DataFrame()
    pairwise_fl = pd.DataFrame()
    if not posthoc_df.empty:
        def make_pairwise(metric_name: str) -> pd.DataFrame:
            x = posthoc_df[posthoc_df["Metric"] == metric_name].copy()
            if x.empty:
                return x
            return pd.DataFrame({
                "Chemistry": x["Chemistry"],
                "Annotator_1": x["Annotator_1"],
                "Annotator_2": x["Annotator_2"],
                "mean_diff_F1_(a1_minus_a2)": x["mean_delta"],
                "p_adj_BH": x["p_adj_BH"],
            })

        pairwise_recall = make_pairwise("Recall_FSM_ISM")
        pairwise_fl = make_pairwise("FL_ratio_FSM_over_FSM_ISM")

    # ---------------------------
    # Write outputs
    # ---------------------------
    metrics_path = outdir / "sqanti_fsmism_metrics.tsv"
    metrics.to_csv(metrics_path, sep="\t", index=False)

    anova_path = outdir / "sqanti_fsmism_anova.tsv"
    anova_df.to_csv(anova_path, sep="\t", index=False)

    posthoc_path = outdir / "sqanti_fsmism_lortia_vs_others.tsv"
    posthoc_df.to_csv(posthoc_path, sep="\t", index=False)

    recall_plot_path = outdir / "SQANTI_Recall_FSM_ISM.tsv"
    fl_plot_path = outdir / "SQANTI_FL_ratio.tsv"
    recall_plot.to_csv(recall_plot_path, sep="\t", index=False)
    fl_plot.to_csv(fl_plot_path, sep="\t", index=False)

    if not pairwise_recall.empty:
        pairwise_recall.to_csv(outdir / "SQANTI_Recall_FSM_ISM__pairwise.tsv", sep="\t", index=False)
    if not pairwise_fl.empty:
        pairwise_fl.to_csv(outdir / "SQANTI_FL_ratio__pairwise.tsv", sep="\t", index=False)

    if args.write_composition:
        comp_df = pd.DataFrame(comp_rows)
        comp_df.to_csv(outdir / "sqanti_structural_category_counts.tsv", sep="\t", index=False)

    print("[OK] wrote:", file=sys.stderr)
    print(f"  {metrics_path}", file=sys.stderr)
    print(f"  {anova_path}", file=sys.stderr)
    print(f"  {posthoc_path}", file=sys.stderr)
    print(f"  {recall_plot_path}", file=sys.stderr)
    print(f"  {fl_plot_path}", file=sys.stderr)
    if (outdir / "SQANTI_Recall_FSM_ISM__pairwise.tsv").exists():
        print(f"  {outdir / 'SQANTI_Recall_FSM_ISM__pairwise.tsv'}", file=sys.stderr)
    if (outdir / "SQANTI_FL_ratio__pairwise.tsv").exists():
        print(f"  {outdir / 'SQANTI_FL_ratio__pairwise.tsv'}", file=sys.stderr)
    if args.write_composition:
        print(f"  {outdir / 'sqanti_structural_category_counts.tsv'}", file=sys.stderr)

if __name__ == "__main__":
    main()