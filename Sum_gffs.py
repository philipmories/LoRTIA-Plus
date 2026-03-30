#!/usr/bin/env python3

import pandas as pd
from argparse import ArgumentParser

def ends(gffs, summary, args, prefixlist):
    for contig in set(gffs["contig"]):
        d = gffs.loc[gffs.contig == contig]
        for strand in set(d["strand"]):
            stranddf = d.loc[d.strand == strand]
            stranddf = stranddf.sort_values(by="start")
            counter = 0
            ID_list = []
            prev = -999
            for index, row in stranddf.iterrows():
                if prev not in range(row["start"] - args.wobble, row["start"]
                                 + args.wobble + 1):
                    counter += 1
                ID_list.append(strand+str(counter))
                prev = row["start"]
            stranddf["ID"] = ID_list

            is_greatest_list = []
            for index, row in stranddf.iterrows():
                maximum = stranddf.loc[stranddf.ID == row["ID"]]["score"].max()
                is_greatest = row["score"] == maximum
                is_greatest_list.append(is_greatest)
            stranddf["is_greatest"] = is_greatest_list

            if (strand == "-" and args.feature == ("tes")) or (strand == "+"
               and args.feature == ("tss")):
                leftmost = True
            else:
                leftmost = False

            is_picked_list = []
            if leftmost:
                for index, row in stranddf.iterrows():
                    if row["is_greatest"]:
                        subdf = stranddf.loc[(stranddf.ID == row["ID"]) & (
                                stranddf.is_greatest == True)]
                        is_picked = row["start"] == subdf["start"].min()
                    else:
                        is_picked = False
                    is_picked_list.append(is_picked)
            else:
                for index, row in stranddf.iterrows():
                    if row["is_greatest"]:
                        subdf = stranddf.loc[(stranddf.ID == row["ID"]) & (
                                stranddf.is_greatest == True)]
                        is_picked = row["start"] == subdf["start"].max()
                    else:
                        is_picked = False
                    is_picked_list.append(is_picked)
            stranddf["is_picked"] = is_picked_list

            chart = stranddf.loc[stranddf.is_picked == True].copy()
            chart.drop_duplicates(subset=["start", "strand"], inplace=True)

            for prefix in prefixlist:
                method = prefix.split("/")[-1]
                df = stranddf.loc[stranddf.method == method]
                got_list = []
                for index, row in chart.iterrows():
                    got_list.append(df.loc[df.ID == row["ID"]]["score"].sum())
                chart[method] = got_list

            summary = pd.concat([summary, chart], ignore_index=True)
    return summary

def main():
    args = parsing()
    prefixlist = []
    with open(args.list_file) as file:
        for line in file:
            prefixlist.append(line.strip())

    col = ["contig","source","feature","start","end","score","strand","frame","info"]

    gffs = pd.DataFrame()
    for prefix in prefixlist:
        d = pd.read_csv(
            prefix + "_" + args.feature + ".gff3",
            sep="\t",
            names=col,
            comment="#",           # ✅ gff3 header/comment sorok kihagyása
            dtype=str
        )
        d["method"] = prefix.split("/")[-1]
        gffs = pd.concat([gffs, d], ignore_index=True)

    # ✅ típusok rendberakása (score lehet '.' is)
    gffs["start"] = pd.to_numeric(gffs["start"], errors="coerce").fillna(0).astype(int)
    gffs["end"]   = pd.to_numeric(gffs["end"],   errors="coerce").fillna(0).astype(int)
    gffs["score"] = pd.to_numeric(gffs["score"], errors="coerce").fillna(0.0)

    summary = pd.DataFrame()

    if args.feature != "intron":
        summary = ends(gffs, summary, args, prefixlist)
    else:
        summary = gffs.copy()
        summary.drop_duplicates(subset=["contig", "start", "end", "strand"], inplace=True)

        for prefix in prefixlist:
            method = prefix.split("/")[-1]
            df = gffs.loc[gffs.method == method].copy()

            got_list = []
            for index, row in summary.iterrows():
                scores = df.loc[
                    (df.start == row["start"]) &
                    (df.end == row["end"]) &
                    (df.contig == row["contig"]) &
                    (df.strand == row["strand"]),
                    "score"
                ].to_numpy()

                # ✅ FIX: numpy array-ra nem lehet 'if scores:'-t írni
                got = float(scores.sum()) if scores.size > 0 else 0.0
                got_list.append(got)

            summary[method] = got_list

    summary.to_csv("{}/{}.txt".format(args.outpath, args.feature), sep="\t", index=False)

def parsing():
    parser = ArgumentParser(description="This is a module that merges gffs.")
    parser.add_argument("list_file", help="Folders list file", metavar="list_file")
    parser.add_argument("outpath", help="Output folder", metavar="outpath")
    parser.add_argument("feature", help="intron/tes/tss", metavar="feature")
    parser.add_argument("-b", "--wobble", dest="wobble", type=int, default=10,
                        help="Wobble for merging TSS/TES (default +/-10).", metavar="[integer]")
    return parser.parse_args()

if __name__== "__main__":
    main()
