#!/usr/bin/env python3

from argparse import ArgumentParser, ArgumentTypeError
import pandas as pd

def str2bool(value):
    """Parse an explicit command-line boolean value."""
    if isinstance(value, bool):
        return value
    normalized = str(value).strip().lower()
    if normalized in {"true", "t", "1", "yes", "y"}:
        return True
    if normalized in {"false", "f", "0", "no", "n"}:
        return False
    raise ArgumentTypeError("expected a boolean value: True or False")



def line_end(df, new_df, feature, sign):
    """
    Prepares gff line for transcript ends (TSS and TES).
    """
    for index, row in df.iterrows():
        i = len(new_df)
        if row["is_qualified"]:    
            new_df.loc[i] = [row["contig"],
                             "LoRTIA", 
                             feature,
                             row["pos"],
                             row["pos"],
                             row["count"],
                             sign,
                             ".",
                             row["poisp"]]
            i += 1

def line_intron(df, new_df, feature):
    """
    Prepares gff line for introns.
    """
    for index, row in df.iterrows():
        i = len(new_df)
        if row["is_qualified"]:    
            new_df.loc[i] = [row["contig"],
                             "LoRTIA", 
                             feature,
                             row["left"] + 1,
                             row["right"] - 1,
                             row["count"],
                             row["strand"],
                             ".",
                             row["consensus"]]
            i += 1

def Gff_creator(args):
    """
    Creates gffs from stats files.
    """
    print("Creating {} gff file...".format(args.feature))
    if not args.output_gff:
        args.output_gff = "{}_{}.gff3".format(args.prefix, args.feature)
    cols = ["contig",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "info"]
    new_df = pd.DataFrame(columns=cols)
    if args.feature == "intron":
        file = "{}_{}.tsv".format(args.prefix, args.feature)
        df = pd.read_csv(file, sep = "\t")
        line_intron(df, new_df, args.feature)
        if args.force_consensus:
            # Keep only the canonical splice-site classes documented by LoRTIA.
            # Do not rely on the string "None": pandas may parse it as NaN,
            # which made this filter environment/version dependent.
            allowed_consensus = {"GT/AG", "GC/AG", "AT/AC"}
            new_df = new_df.loc[new_df["info"].isin(allowed_consensus)]
    else:
        if args.feature == "tss":
            filepos = "{}_l5_{}.tsv".format(args.prefix, args.feature)
            fileneg = "{}_r5_{}.tsv".format(args.prefix, args.feature)
        else:
            filepos = "{}_r3_{}.tsv".format(args.prefix, args.feature)
            fileneg = "{}_l3_{}.tsv".format(args.prefix, args.feature)
        dfpos = pd.read_csv(filepos, sep = "\t")
        dfneg = pd.read_csv(fileneg, sep = "\t")
    
        line_end(dfpos, new_df, args.feature, "+")
        line_end(dfneg, new_df, args.feature, "-")
        if args.significance == "poisson":
            try:
                alpha = 0.05 / len(new_df)
                new_df = new_df.loc[new_df["info"] < alpha]
            except:
                print("WARNING! There were no significant {} features in this\
                      sample.".format(args.feature))
    new_df = new_df.sort_values(by=['end'])
    new_df = new_df.sort_values(by=['start'])
    new_df = new_df.sort_values(by=['contig'])
    new_df.to_csv(args.output_gff, sep="\t", index = False, header = False)


###############################################################################
###                             Main function                               ###
###############################################################################

def main():
    args = parsing()
    Gff_creator(args)

def parsing():
    parser = ArgumentParser(description="This is the third module of \
                            LoRTIA, a Long-read RNA-Seq Transcript Isofom \
                            Annotator. This module creates gff files from\
                            the feature statistics.")
    parser.add_argument("prefix",
                        help="The path and the prefix of the statistics, \
                        which are to be used for the gff.",
                        metavar="prefix")
    parser.add_argument("feature",
                        help="The type of feature for which the gff is \
                        generated. Options include 'tss' for transcriptional \
                        start sites, 'tes' for transcriptional end sites and \
                        'intron' for introns.",
                        metavar="feature")
    parser.add_argument("-o", "--output_gff",
                        help="The output file that is to be generated.",
                        default=False,
                        metavar="[file]")
    parser.add_argument("-s", "--significance",
                        dest="significance",
                        help="Optional Poisson significance filter [poisson]. "
                             "The top-level LoRTIA Plus workflow applies it to "
                             "TSS generation only; the standalone module applies "
                             "it to the requested endpoint feature. Default: "
                             "disabled.",
                        default=False,
                        metavar="[string]")
    parser.add_argument("-f", "--force_consensus",
                        dest="force_consensus",
                        help="Require GT-AG, GC-AG or AT-AC consensus splice-site "
                             "pairs for intron output. Pass True or False. "
                             "Default: False.",
                        type=str2bool,
                        default=False,
                        metavar="[bool]")
    return parser.parse_args()


if __name__== "__main__":
    main()
