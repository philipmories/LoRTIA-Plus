#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Sum_gffs_fast.py (FULL: intron + tes + tss)

- list_file: each line is a prefix like:
    .../Guinea_Pig_4hpi_1/Guinea_Pig_4hpi_1
  and points to files like:
    <prefix>_intron.gff3
    <prefix>_tss.gff3
    <prefix>_tes.gff3 / <prefix>_not_ts_tes.gff3 / <prefix>_ts_tes.gff3

- infer_method(): derives method/sample name from prefix path (your Guinea_Pig_* stays intact)

- intron: aggregates per method across prefixes (e.g. chromosomes) into ONE column/method.

- tes/tss: wobble clustering, representative site chosen per cluster, then per-method score sum per cluster.
  IMPORTANT FIX: we compute every cluster column for ALL methods, so columns never “disappear”.

- missing values in tes/tss matrix can be written as 0 or NA (choose via --missing).
"""

import os
import re
import sys
import pandas as pd
from argparse import ArgumentParser
from typing import List, Dict, Optional

COLS = ["contig", "source", "feature", "start", "end", "score", "strand", "frame", "info"]


# ---------------- CLI ----------------

def parsing():
    p = ArgumentParser(description="Merge LoRTIA feature GFF3s from multiple prefixes.")
    p.add_argument("list_file", help="Prefixes list file (one prefix per line).")
    p.add_argument("outpath", help="Output folder.")
    p.add_argument("feature", help="intron/tes/tss.")
    p.add_argument("-b", "--wobble", type=int, default=10, help="Wobble for merging TSS/TES (default +/-10).")

    p.add_argument(
        "--intron_agg", choices=["sum", "max", "first"], default="sum",
        help="Aggregation for identical intron keys within one method (default: sum)."
    )

    p.add_argument(
        "--tes_input", choices=["auto", "not_ts", "tes", "ts"], default="auto",
        help="TES input: auto=prefer *_not_ts_tes.gff3 then fallback to *_tes.gff3 (default: auto)."
    )

    p.add_argument(
        "--missing", choices=["0", "na"], default="0",
        help="How to encode missing method support in a TES/TSS cluster: 0 or NA (default: 0)."
    )

    p.add_argument(
        "--show-methods", action="store_true",
        help="Print inferred method for each prefix and exit."
    )
    return p.parse_args()


# ---------------- helpers ----------------

def read_prefixes(list_file: str) -> List[str]:
    out = []
    with open(list_file, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip()
            if s:
                out.append(s.rstrip("/"))
    return out


def _is_contig_token(tok: str) -> bool:
    return bool(re.fullmatch(r"chr([0-9]{1,3}|[XYM]|MT|Un.*|UN.*|GL.*|KI.*|JH.*)", tok))


def infer_method(prefix: str) -> str:
    """
    Walk up from end and skip contig-ish parts and contigs folder.
    In your case (.../Guinea_Pig_4hpi_1/Guinea_Pig_4hpi_1) -> returns Guinea_Pig_4hpi_1
    """
    parts = prefix.rstrip("/").split("/")
    if not parts:
        return "NA"

    i = len(parts) - 1
    while i >= 0 and (_is_contig_token(parts[i]) or parts[i] in {"contig", "contigs"}):
        i -= 1
    if i < 0:
        return "NA"

    tok = parts[i]
    # replicate folder pattern (digits)
    if tok.isdigit() and i - 1 >= 0:
        return f"{parts[i-1]}_rep{tok}"
    return tok


def group_prefixes_by_method(prefixlist: List[str]) -> Dict[str, List[str]]:
    d: Dict[str, List[str]] = {}
    for pref in prefixlist:
        m = infer_method(pref)
        d.setdefault(m, []).append(pref)
    return d


def has_data_lines(path: str) -> bool:
    """True if file has at least one non-empty, non-comment line."""
    try:
        with open(path, "rt", encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.strip() and not line.startswith("#"):
                    return True
    except FileNotFoundError:
        return False
    return False


def pick_gff_path(pref: str, feature: str, tes_input: str = "auto") -> str:
    if feature != "tes":
        candidates = [f"{pref}_{feature}.gff3"]
    else:
        if tes_input == "auto":
            candidates = [f"{pref}_not_ts_tes.gff3", f"{pref}_tes.gff3"]
        elif tes_input == "not_ts":
            candidates = [f"{pref}_not_ts_tes.gff3"]
        elif tes_input == "ts":
            candidates = [f"{pref}_ts_tes.gff3"]
        else:
            candidates = [f"{pref}_tes.gff3"]

    for p in candidates:
        if os.path.exists(p):
            return p
    raise FileNotFoundError(f"Missing all candidates for {pref} ({feature}, tes_input={tes_input}): {candidates}")


def safe_read_gff_any(path: str) -> pd.DataFrame:
    """Read full 9-column GFF3 (skip comments). Returns empty df if comment-only/empty."""
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    if not has_data_lines(path):
        return pd.DataFrame(columns=COLS)

    try:
        df = pd.read_csv(path, sep="\t", names=COLS, comment="#", header=None, engine="c")
    except (pd.errors.EmptyDataError, StopIteration) as e:
        return pd.DataFrame(columns=COLS)

    if df is None or df.empty:
        return pd.DataFrame(columns=COLS)
    return df


def safe_read_intron_cols(path: str) -> pd.DataFrame:
    """
    Read intron minimal columns: contig, start, end, score, strand.
    Returns empty df if comment-only/empty.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    if not has_data_lines(path):
        return pd.DataFrame(columns=["contig", "start", "end", "score", "strand"])

    try:
        df = pd.read_csv(
            path,
            sep="\t",
            names=COLS,
            comment="#",
            header=None,
            usecols=[0, 3, 4, 5, 6],  # contig, start, end, score, strand
            engine="c",
        )
    except (pd.errors.EmptyDataError, StopIteration, IndexError) as e:
        return pd.DataFrame(columns=["contig", "start", "end", "score", "strand"])

    if df is None or df.empty:
        return pd.DataFrame(columns=["contig", "start", "end", "score", "strand"])

    df.columns = ["contig", "start", "end", "score", "strand"]
    return df


# ---------------- intron merge ----------------

def fast_merge_introns(prefixlist: List[str], outpath: str, intron_agg: str = "sum"):
    keys = ["contig", "start", "end", "strand"]
    method_map = group_prefixes_by_method(prefixlist)
    summary = None

    # üres (de nevezett) MultiIndex, hogy a join sose omoljon el
    empty_idx = pd.MultiIndex.from_tuples([], names=keys)

    for method, prefs in sorted(method_map.items()):
        series_total = None

        for pref in prefs:
            path = f"{pref}_intron.gff3"
            if not os.path.exists(path):
                raise FileNotFoundError(f"Missing: {path}")

            raw = safe_read_intron_cols(path)
            if raw.empty:
                continue

            raw["start"] = pd.to_numeric(raw["start"], errors="coerce").fillna(0).astype("int32")
            raw["end"] = pd.to_numeric(raw["end"], errors="coerce").fillna(0).astype("int32")
            raw["score"] = pd.to_numeric(raw["score"], errors="coerce").fillna(0.0).astype("float32")

            gb = raw.groupby(keys, sort=False, observed=True)["score"]
            if intron_agg == "max":
                s = gb.max()
            elif intron_agg == "first":
                s = gb.first()
            else:
                s = gb.sum()

            # biztosítsuk az index nevet (pandas néha eldobhatja bizonyos műveleteknél)
            s.index = s.index.set_names(keys)

            if series_total is None:
                series_total = s
            else:
                if intron_agg == "max":
                    series_total = series_total.combine(s, func=max, fill_value=0.0)
                elif intron_agg == "first":
                    series_total = series_total.combine_first(s)
                else:
                    series_total = series_total.add(s, fill_value=0.0)

        # 🔧 FIX: ha nincs intron ebben a mintában, akkor is legyen oszlop — nevezett MultiIndex-szel
        if series_total is None:
            series_total = pd.Series([], index=empty_idx, dtype="float32", name=method)
        else:
            series_total.index = series_total.index.set_names(keys)
            series_total.name = method

        summary = series_total.to_frame() if summary is None else summary.join(series_total, how="outer")

    if summary is None:
        summary = pd.DataFrame(index=empty_idx)

    summary = summary.fillna(0.0).reset_index()
    os.makedirs(outpath, exist_ok=True)
    summary.to_csv(os.path.join(outpath, "intron.txt"), sep="\t", index=False)

# ---------------- TES/TSS merge ----------------

def merge_tss_tes(prefixlist: List[str], outpath: str, feature: str, wobble: int, tes_input: str, missing_mode: str):
    method_map = group_prefixes_by_method(prefixlist)
    all_methods = sorted(method_map.keys())

    frames = []
    for method, prefs in sorted(method_map.items()):
        for pref in prefs:
            path = pick_gff_path(pref, feature, tes_input=tes_input)
            df = safe_read_gff_any(path)
            if df.empty:
                continue
            df["method"] = method
            frames.append(df)

    gffs = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=COLS + ["method"])
    if gffs.empty:
        os.makedirs(outpath, exist_ok=True)
        gffs.to_csv(os.path.join(outpath, f"{feature}.txt"), sep="\t", index=False)
        return

    gffs["feature"] = gffs["feature"].astype(str).str.lower()
    gffs["start"] = pd.to_numeric(gffs["start"], errors="coerce").fillna(0).astype(int)
    gffs["end"] = pd.to_numeric(gffs["end"], errors="coerce").fillna(0).astype(int)
    gffs["score"] = pd.to_numeric(gffs["score"], errors="coerce").fillna(0.0).astype(float)

    # keep only requested feature (tes/tss)
    gffs = gffs.loc[gffs["feature"] == feature].copy()
    if gffs.empty:
        os.makedirs(outpath, exist_ok=True)
        gffs.to_csv(os.path.join(outpath, f"{feature}.txt"), sep="\t", index=False)
        return

    summary_chunks = []
    fill_value = 0.0 if missing_mode == "0" else pd.NA

    for contig in gffs["contig"].dropna().unique():
        d = gffs.loc[gffs["contig"] == contig]
        for strand in d["strand"].dropna().unique():
            stranddf = d.loc[d["strand"] == strand].sort_values(by="start").copy()
            if stranddf.empty:
                continue

            # cluster by walking sorted starts (chain within wobble)
            ids = []
            counter = 0
            prev = -10**12
            for st in stranddf["start"].tolist():
                if st not in range(prev - wobble, prev + wobble + 1):
                    counter += 1
                ids.append(counter)
                prev = st
            stranddf["ID"] = ids

            # max score per cluster
            stranddf["is_greatest"] = stranddf["score"].eq(
                stranddf.groupby("ID", sort=False, observed=True)["score"].transform("max")
            )

            # tie-break: leftmost for (TES on '-') or (TSS on '+'), else rightmost
            leftmost = (strand == "-" and feature == "tes") or (strand == "+" and feature == "tss")
            sub = stranddf.loc[stranddf["is_greatest"]].copy()
            if sub.empty:
                continue

            chosen = sub.groupby("ID", sort=False, observed=True)["start"].transform("min" if leftmost else "max")

            stranddf["is_picked"] = False
            stranddf.loc[sub.index, "is_picked"] = sub["start"].eq(chosen)

            chart = stranddf.loc[stranddf["is_picked"]].copy()
            chart.drop_duplicates(subset=["contig", "start", "strand"], inplace=True)

            # FIX: compute columns for ALL methods (not only those present in this chunk)
            for m in all_methods:
                sums = stranddf.loc[stranddf["method"] == m].groupby("ID", sort=False, observed=True)["score"].sum()
                mapped = chart["ID"].map(sums)
                if missing_mode == "0":
                    chart[m] = mapped.fillna(0.0)
                else:
                    chart[m] = mapped.astype("Float64")  # allow NA

            summary_chunks.append(chart)

    out = pd.concat(summary_chunks, ignore_index=True) if summary_chunks else pd.DataFrame()

    # Column order: meta first, then all sample columns
    meta_cols = [c for c in (COLS + ["method", "ID", "is_greatest", "is_picked"]) if c in out.columns]
    sample_cols = [c for c in all_methods if c in out.columns]
    rest = [c for c in out.columns if c not in set(meta_cols + sample_cols)]
    out = out[meta_cols + rest + sample_cols]

    os.makedirs(outpath, exist_ok=True)
    out.to_csv(os.path.join(outpath, f"{feature}.txt"), sep="\t", index=False)


# ---------------- main ----------------

def main():
    args = parsing()
    prefixlist = read_prefixes(args.list_file)
    if not prefixlist:
        print("[ERROR] list_file is empty.", file=sys.stderr)
        raise SystemExit(2)

    if args.show_methods:
        for p in prefixlist[:2000]:
            print(f"{p}\t=>\t{infer_method(p)}")
        return

    feature = args.feature.lower().strip()
    if feature not in {"intron", "tes", "tss"}:
        raise ValueError("feature must be one of: intron / tes / tss")

    os.makedirs(args.outpath, exist_ok=True)

    if feature == "intron":
        fast_merge_introns(prefixlist, args.outpath, intron_agg=args.intron_agg)
    else:
        merge_tss_tes(prefixlist, args.outpath, feature=feature, wobble=args.wobble,
                      tes_input=args.tes_input, missing_mode=args.missing)


if __name__ == "__main__":
    main()
