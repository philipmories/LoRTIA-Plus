#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
SQANTI classification összefoglaló egy database/manifest TSV alapján.

Database TSV (pl. database_with_TES.txt) soronként:
    Annotator   Chemical    Cell    GTF   Classification   TSS   TES

Kimenet: chemistry × cell × structural_category × annotator összesítés
- coding_rate                      = coding == "coding"
- within_cage_rate                 = within_CAGE_peak == TRUE
- polyA_site_rate                  = within_polyA_site == TRUE
- polyA_motif_rate                 = polyA_motif != NA
- polyA_site_and_motif_rate        = polyA_motif != NA & within_polyA_site == TRUE
(+ n_isoforms, weight_sum, források)
"""

import argparse, os, sys
import pandas as pd
import numpy as np

CAT_MAP = {
    "full-splice_match"      : "FSM",
    "incomplete-splice_match": "ISM",
    "novel_in_catalog"       : "NIC",
    "novel_not_in_catalog"   : "NNC",
}

TRUE_SET = {"TRUE", "True", "true", True, 1, "1"}


def to_bool(x):
    if pd.isna(x):
        return False
    return str(x) in TRUE_SET


def parse_args():
    ap = argparse.ArgumentParser(
        description="Egyszerű SQANTI summarizer (CAGE/polyA/coding ráták) "
                    "database/manifest TSV alapján."
    )
    ap.add_argument("-i", "--index", required=True,
                    help="Database/manifest TSV (pl. database_with_TES.txt).")
    ap.add_argument("-o", "--out", required=True,
                    help="Output TSV path.")

    # oszlopnevek a database-ben
    ap.add_argument("--annot-col", default="Annotator",
                    help="Annotátor oszlop neve a database-ben (alap: Annotator).")
    ap.add_argument("--chem-col", default="Chemical",
                    help="Kémia oszlop neve (alap: Chemical).")
    ap.add_argument("--cell-col", default="Cell",
                    help="Sejtvonal oszlop neve (alap: Cell).")
    ap.add_argument("--classif-col", default="Classification",
                    help="SQANTI classification oszlop neve (alap: Classification).")
    ap.add_argument("--gtf-col", default="GTF",
                    help="(opcionális) GTF path oszlop neve (alap: GTF).")
    ap.add_argument("--tss-col", default="TSS",
                    help="(opcionális) TSS GTF path oszlop neve (alap: TSS).")
    ap.add_argument("--tes-col", default="TES",
                    help="(opcionális) TES GTF path oszlop neve (alap: TES).")

    ap.add_argument("--weight", choices=["none", "iso_exp", "gene_exp"],
                    default="none",
                    help="Súlyozás iso_exp/gene_exp alapján (alap: none=1.0 minden sor).")

    return ap.parse_args()


def main():
    args = parse_args()

    # database_with_TES.txt beolvasása
    try:
        db = pd.read_csv(args.index, sep="\t", na_values=["NA"], low_memory=False)
    except Exception as e:
        sys.exit(f"[ERROR] Nem tudtam beolvasni az index/database fájlt ({args.index}): {e}")

    for col in [args.annot_col, args.chem_col, args.cell_col, args.classif_col]:
        if col not in db.columns:
            sys.exit(f"[ERROR] Hiányzó oszlop a database-ben: {col}")

    rows = []
    scanned = parsed = 0

    for _, row in db.iterrows():
        scanned += 1
        annot = str(row[args.annot_col])
        chem  = str(row[args.chem_col])
        cell  = str(row[args.cell_col])
        classif_path = str(row[args.classif_col])

        gtf_path = row[args.gtf_col] if args.gtf_col in db.columns else np.nan
        tss_path = row[args.tss_col] if args.tss_col in db.columns else np.nan
        tes_path = row[args.tes_col] if args.tes_col in db.columns else np.nan

        if not classif_path or classif_path == "nan":
            continue
        if not os.path.exists(classif_path):
            print(f"[WARN] Classification fájl nem létezik: {classif_path}", file=sys.stderr)
            continue

        # SQANTI classification beolvasása
        try:
            df = pd.read_csv(
                classif_path,
                sep="\t",
                na_values=["NA"],
                low_memory=False,
                dtype={
                    "within_CAGE_peak": "string",
                    "within_polyA_site": "string",
                    "polyA_motif_found": "string",
                    "coding": "string"
                }
            )
        except Exception as e:
            print(f"[WARN] Nem tudtam beolvasni a classification-t ({classif_path}): {e}", file=sys.stderr)
            continue

        if "structural_category" not in df.columns:
            print(f"[WARN] Missing 'structural_category' in {classif_path}; skipping.", file=sys.stderr)
            continue

        parsed += 1

        # biztosítsuk, hogy a szükséges oszlopok meglegyenek
        for col in [
            "within_CAGE_peak", "within_polyA_site",
            "polyA_motif", "polyA_motif_found", "coding"
        ]:
            if col not in df.columns:
                df[col] = np.nan

        # Súlyok
        if args.weight == "iso_exp" and "iso_exp" in df.columns:
            w = pd.to_numeric(df["iso_exp"], errors="coerce").fillna(1.0)
        elif args.weight == "gene_exp" and "gene_exp" in df.columns:
            w = pd.to_numeric(df["gene_exp"], errors="coerce").fillna(1.0)
        else:
            w = pd.Series(1.0, index=df.index)

        # Bool flag-ek
        coding_flag = (df["coding"].astype(str) == "coding")

        within_cage  = df["within_CAGE_peak"].apply(to_bool)
        within_polyA = df["within_polyA_site"].apply(to_bool)

        # polyA_motif != NA / nem üres
        motif_non_na = (~df["polyA_motif"].isna()) & (df["polyA_motif"].astype(str) != "NA") \
                       & (df["polyA_motif"].astype(str) != "")

        motif_and_site = motif_non_na & within_polyA

        # aggregálás structural_category szerint
        tmp = pd.DataFrame({
            "structural_category": df["structural_category"].astype(str),
            "coding_flag": coding_flag.astype(int),
            "within_cage": within_cage.astype(int),
            "within_polyA": within_polyA.astype(int),
            "motif_non_na": motif_non_na.astype(int),
            "motif_and_site": motif_and_site.astype(int),
            "w": w
        })

        for cat, g in tmp.groupby("structural_category"):
            cat3 = CAT_MAP.get(cat, cat)
            total_w = g["w"].sum()
            if total_w <= 0:
                continue

            def wrate(col):
                return float((g[col].astype(float) * g["w"]).sum() / total_w)

            rows.append({
                "chemistry": chem,
                "cell": cell,
                "category": cat3,
                "annotator": annot,

                # kért ráták
                "coding_rate": wrate("coding_flag"),          # coding / összes
                "within_cage_rate": wrate("within_cage"),     # TRUE / összes
                "polyA_site_rate": wrate("within_polyA"),     # TRUE / összes
                "polyA_motif_rate": wrate("motif_non_na"),    # polyA_motif != NA / összes
                "polyA_site_and_motif_rate": wrate("motif_and_site"),

                "n_isoforms": int(len(g)),
                "weight_sum": float(total_w),

                "source_classif": classif_path,
                "source_gtf": gtf_path,
                "source_tss": tss_path,
                "source_tes": tes_path,
            })

    if not rows:
        sys.exit("[ERROR] Nem jött létre egyetlen sor sem. "
                 "Ellenőrizd a database-t és a classification fájlokat!")

    out = pd.DataFrame(rows)

    cols = [
        "chemistry", "cell", "category", "annotator",
        "coding_rate",
        "within_cage_rate",
        "polyA_site_rate",
        "polyA_motif_rate",
        "polyA_site_and_motif_rate",
        "n_isoforms", "weight_sum",
        "source_classif", "source_gtf", "source_tss", "source_tes"
    ]
    for c in cols:
        if c not in out.columns:
            out[c] = np.nan
    out = out.reindex(columns=cols)

    out.sort_values(["cell", "chemistry", "category", "annotator"], inplace=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)

    print(f"[INFO] Database sorok: {len(db)}, feldolgozott classification-ok: {parsed}")
    print(f"[INFO] Wrote: {args.out}")


if __name__ == "__main__":
    main()
