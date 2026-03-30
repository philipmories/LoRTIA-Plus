#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import numpy as np


REQUIRED_COLS = {"structural_category", "length", "coding", "polyA_motif_found"}

# opcionális rövidítések (ha a structural_category-t röviden is szeretnéd látni)
CAT_SHORT = {
    "full-splice_match": "FSM",
    "incomplete-splice_match": "ISM",
    "novel_in_catalog": "NIC",
    "novel_not_in_catalog": "NNC",
    "fusion": "FUS",
    "antisense": "AS",
    "genic": "GENIC",
    "genic_intron": "GENIC_INTRON",
}


def detect_delimiter(header_line: str) -> str:
    # SQANTI classification jellemzően TSV
    if "\t" in header_line:
        return "\t"
    if "," in header_line:
        return ","
    return "\t"


def try_read_header(path: Path) -> Tuple[Optional[List[str]], Optional[str]]:
    try:
        with path.open("r", encoding="utf-8", errors="replace") as f:
            line = f.readline().strip("\n\r")
        if not line:
            return None, None
        delim = detect_delimiter(line)
        cols = [c.strip() for c in line.split(delim)]
        return cols, delim
    except Exception:
        return None, None


def find_sqanti_classification_like_files(root: Path) -> List[Tuple[Path, str]]:
    hits: List[Tuple[Path, str]] = []
    for p in root.rglob("*"):
        if not p.is_file():
            continue
        # kicsit szűkítünk: tipikus kiterjesztések
        if p.suffix.lower() not in {".txt", ".tsv", ".csv"}:
            continue

        cols, delim = try_read_header(p)
        if cols is None or delim is None:
            continue

        colset = set(cols)
        if REQUIRED_COLS.issubset(colset):
            hits.append((p, delim))

    return sorted(hits, key=lambda x: str(x[0]))


def parse_bool(x) -> Optional[bool]:
    if x is None:
        return None
    s = str(x).strip().lower()
    if s in {"", "nan", "none"}:
        return None
    if s in {"true", "t", "1", "yes", "y"}:
        return True
    if s in {"false", "f", "0", "no", "n"}:
        return False
    return None


def classify_coding(x) -> Optional[str]:
    """
    SQANTI 'coding' oszlop gyakran: 'coding' vagy 'non_coding'.
    Itt robusztusan próbáljuk besorolni.
    """
    if x is None:
        return None
    s = str(x).strip().lower()
    if s in {"", "nan", "none"}:
        return None

    # noncoding előbb, hogy a 'non_coding' ne essen bele a 'coding' feltételbe
    if s.startswith("non") or "non_coding" in s or "noncoding" in s or "non-coding" in s:
        return "noncoding"
    if "coding" in s or s in {"c"}:
        return "coding"
    return None


def summarize_one_df(df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    # csak a releváns oszlopokkal dolgozunk (gyorsabb és nem húzzuk be az ORF_seq-et stb.)
    df = df.copy()

    # hossz
    df["length_num"] = pd.to_numeric(df["length"], errors="coerce")

    # kategória
    df["cat_full"] = df["structural_category"].astype(str).str.strip()
    df["cat_norm"] = df["cat_full"].str.lower()

    # rövidítés (ha nincs mapping, marad az eredeti)
    df["cat_short"] = df["cat_norm"].map(CAT_SHORT).fillna(df["cat_full"])

    # coding
    df["coding_class"] = df["coding"].apply(classify_coding)

    # polyA
    df["polya_bool"] = df["polyA_motif_found"].apply(parse_bool)

    # group summary
    def group_stats(g: pd.DataFrame) -> Dict[str, object]:
        n = int(len(g))
        mean_len = float(np.nanmean(g["length_num"].values)) if n > 0 else np.nan

        coding_n = int((g["coding_class"] == "coding").sum())
        noncoding_n = int((g["coding_class"] == "noncoding").sum())
        coding_denom = coding_n + noncoding_n
        coding_ratio = (coding_n / coding_denom) if coding_denom > 0 else np.nan

        polya_true = int((g["polya_bool"] == True).sum())
        polya_false = int((g["polya_bool"] == False).sum())
        polya_denom = polya_true + polya_false
        polya_ratio = (polya_true / polya_denom) if polya_denom > 0 else np.nan

        return {
            "n": n,
            "mean_length": mean_len,
            "coding_n": coding_n,
            "noncoding_n": noncoding_n,
            "coding_ratio": coding_ratio,
            "polya_true_n": polya_true,
            "polya_false_n": polya_false,
            "polya_ratio": polya_ratio,
        }

    rows = []
    for cat, g in df.groupby("cat_short", dropna=False):
        st = group_stats(g)
        rows.append({
            "structural_category": cat,
            "structural_category_long": g["cat_full"].iloc[0],
            **st
        })

    # ALL sor
    all_st = group_stats(df)
    rows.append({
        "structural_category": "ALL",
        "structural_category_long": "ALL",
        **all_st
    })

    summary_df = pd.DataFrame(rows)

    # boxplot statok (fájl+ kategória): min,q1,median,q3,max,mean,n
    box_rows = []
    for cat, g in df.groupby("cat_short", dropna=False):
        x = g["length_num"].dropna().values
        if len(x) == 0:
            box_rows.append({
                "structural_category": cat,
                "n": int(len(g)),
                "min": np.nan,
                "q1": np.nan,
                "median": np.nan,
                "q3": np.nan,
                "max": np.nan,
                "mean": np.nan,
            })
            continue

        q1, med, q3 = np.quantile(x, [0.25, 0.5, 0.75])
        box_rows.append({
            "structural_category": cat,
            "n": int(len(x)),
            "min": float(np.min(x)),
            "q1": float(q1),
            "median": float(med),
            "q3": float(q3),
            "max": float(np.max(x)),
            "mean": float(np.mean(x)),
        })

    # ALL box stat is
    x_all = df["length_num"].dropna().values
    if len(x_all) > 0:
        q1, med, q3 = np.quantile(x_all, [0.25, 0.5, 0.75])
        box_rows.append({
            "structural_category": "ALL",
            "n": int(len(x_all)),
            "min": float(np.min(x_all)),
            "q1": float(q1),
            "median": float(med),
            "q3": float(q3),
            "max": float(np.max(x_all)),
            "mean": float(np.mean(x_all)),
        })
    else:
        box_rows.append({
            "structural_category": "ALL",
            "n": 0,
            "min": np.nan,
            "q1": np.nan,
            "median": np.nan,
            "q3": np.nan,
            "max": np.nan,
            "mean": np.nan,
        })

    box_df = pd.DataFrame(box_rows)
    return summary_df, box_df


def main():
    ap = argparse.ArgumentParser(
        description="Recursively summarize SQANTI classification-like TSVs into category stats + boxplot stats."
    )
    ap.add_argument("-i", "--input", required=True, help="Root folder to search recursively.")
    ap.add_argument("-o", "--output", required=True, help="Output folder for TSV results.")
    args = ap.parse_args()

    in_root = Path(args.input).resolve()
    out_dir = Path(args.output).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    files = find_sqanti_classification_like_files(in_root)
    if not files:
        raise SystemExit(
            f"No SQANTI classification-like files found under: {in_root}\n"
            f"Needed header columns: {sorted(REQUIRED_COLS)}"
        )

    all_summary = []
    all_box = []

    for path, delim in files:
        rel = path.relative_to(in_root).as_posix()

        # csak a szükséges oszlopokat olvassuk (különösen fontos az ORF_seq miatt!)
        df = pd.read_csv(
            path,
            sep=delim,
            usecols=lambda c: c in REQUIRED_COLS,
            dtype=str,
            low_memory=False,
        )

        summary_df, box_df = summarize_one_df(df)

        summary_df.insert(0, "file", rel)
        box_df.insert(0, "file", rel)

        all_summary.append(summary_df)
        all_box.append(box_df)

    summary_out = pd.concat(all_summary, ignore_index=True)
    box_out = pd.concat(all_box, ignore_index=True)

    # Rendezés: file, majd ismert kategóriák preferált sorrendben, aztán többi
    preferred = ["FSM", "ISM", "NIC", "NNC", "FUS", "AS", "GENIC", "GENIC_INTRON", "ALL"]
    cat_order = {c: i for i, c in enumerate(preferred)}
    summary_out["_cat_ord"] = summary_out["structural_category"].map(lambda x: cat_order.get(str(x), 999))
    box_out["_cat_ord"] = box_out["structural_category"].map(lambda x: cat_order.get(str(x), 999))

    summary_out = summary_out.sort_values(["file", "_cat_ord", "structural_category"]).drop(columns=["_cat_ord"])
    box_out = box_out.sort_values(["file", "_cat_ord", "structural_category"]).drop(columns=["_cat_ord"])

    summary_path = out_dir / "sqanti_category_summary.tsv"
    box_path = out_dir / "sqanti_length_boxplot_stats.tsv"

    summary_out.to_csv(summary_path, sep="\t", index=False)
    box_out.to_csv(box_path, sep="\t", index=False)

    print(f"Wrote: {summary_path}")
    print(f"Wrote: {box_path}")
    print(f"Processed files: {len(files)}")


if __name__ == "__main__":
    main()
