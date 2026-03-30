#!/usr/bin/env python3

import pysam
import polars as pl
from argparse import ArgumentParser
from ast import literal_eval
import pandas as pd
import numpy as np


def build_feature_index(gff_df):
    """
    Build an index for TSS/TES features:
    (contig, strand) -> sorted numpy array of start positions.
    """
    index = {}
    if gff_df is None:
        return index
    grouped = gff_df.group_by(["contig", "strand"]).agg(
        pl.col("start").sort().alias("starts")
    )
    for row in grouped.iter_rows(named=True):
        key = (row["contig"], row["strand"])
        index[key] = np.array(row["starts"], dtype=int)
    return index


def build_intron_index(intron_gff):
    """
    Build a fast lookup set for introns:
    (contig, strand, intron_start, intron_end)
    where intron_start = gff_start - 1 and intron_end = gff_end + 1
    to match the coordinates stored in the SAM tags.
    """
    intron_index = set()
    if intron_gff is None:
        return intron_index
    for row in intron_gff.iter_rows(named=True):
        contig = row["contig"]
        strand = row["strand"]
        intron_start = row["start"] - 1
        intron_end = row["end"] + 1
        intron_index.add((contig, strand, intron_start, intron_end))
    return intron_index


def read_ends_indexed(feature_index, read, contig, strand, end, args, feature):
    """
    Checks whether the read has accepted TES or TSS features using the
    pre-built feature_index.
    feature: "tss" or "tes" -> we set tag "ts"/"te".
    """
    key = (contig, strand)
    starts = feature_index.get(key)
    tag_name = feature[:-1]  # "tss" -> "ts", "tes" -> "te"
    if starts is None or len(starts) == 0:
        read.set_tag(tag_name, -1, "i")
        return read, False

    left = end - args.wobble
    right = end + args.wobble

    i = int(np.searchsorted(starts, left, side="left"))
    while i < len(starts) and starts[i] <= right:
        position = int(starts[i])
        read.set_tag(tag_name, position, "i")
        return read, True

    read.set_tag(tag_name, -1, "i")
    return read, False


def read_introns_fast(read, introns, args, contig, strand, intron_index):
    """
    STRICT intron handling:
      - if ANY intron from the read is not found in the accepted intron GFF,
        tag it as 'ga' and DO NOT let the read contribute to transcript calling.
      - also allow matching against '.' strand in the intron GFF (if present).
    """
    if not introns:
        read.set_tag("in", None)
        return read, False

    if isinstance(introns[0], int):
        if len(introns) == 2:
            pairs = [(introns[0], introns[1])]
        else:
            pairs = list(zip(introns[0::2], introns[1::2]))
    else:
        pairs = list(introns)

    real_introns = ()
    missing_any = False

    for start_i, end_i in pairs:
        key = (contig, strand, start_i, end_i)
        key_any = (contig, ".", start_i, end_i)

        if (key in intron_index) or (key_any in intron_index):
            real_introns += (start_i, end_i)
        else:
            missing_any = True
            read.set_tag("ga", str((start_i, end_i)), "Z")

    if missing_any:
        read.set_tag("in", None)
        return read, False

    read.set_tag("in", str(real_introns), "Z")
    return read, True


def shifted_tss_match(tss_index, read, contig, strand, tss_pos, args):
    """
    dRNA TSS SNAP:
    1) try shifted position (strand-aware) within wobble
    2) if not found -> fallback to raw position
    """
    shift = int(args.drna_tss_shift)
    if shift <= 0:
        return read_ends_indexed(tss_index, read, contig, strand, tss_pos, args, "tss")

    if strand == "+":
        shifted = tss_pos - shift
    else:
        shifted = tss_pos + shift

    if shifted < 1:
        shifted = 1

    read, found = read_ends_indexed(tss_index, read, contig, strand, shifted, args, "tss")
    if found:
        return read, True

    # fallback
    return read_ends_indexed(tss_index, read, contig, strand, tss_pos, args, "tss")


def bam_iterator(bam, tr_dict, tss_index, tes_index, intron_index, outbam, args):
    for read in bam:
        tss_found = tes_found = intron_found = gap = False
        contig = read.reference_name

        if read.is_reverse:
            strand = "-"
            tss_tag = "r5"
            tss_pos = read.reference_end
            tes_tag = "l3"
            tes_pos = read.reference_start + 1
        else:
            strand = "+"
            tss_tag = "l5"
            tss_pos = read.reference_start + 1
            tes_tag = "r3"
            tes_pos = read.reference_end

        # --- adapter/QC tags (may be missing) ---
        try:
            tss_info = read.get_tag(tss_tag).split(",")[3]
        except KeyError:
            tss_info = "missing"
        except Exception:
            tss_info = "missing"

        try:
            tes_info = read.get_tag(tes_tag).split(",")[3]
        except KeyError:
            tes_info = "missing"
        except Exception:
            tes_info = "missing"

        # --- TSS ---
        # non-dRNA: keep original gating
        # dRNA (--drna): always attempt TSS matching with shift+fallback
        if args.drna:
            read, tss_found = shifted_tss_match(tss_index, read, contig, strand, tss_pos, args)
        else:
            if (
                tss_info == "correct"
                or tss_info == "potential template switching"
                or tes_info == "out of place"
            ):
                read, tss_found = read_ends_indexed(
                    tss_index, read, contig, strand, tss_pos, args, "tss"
                )

        # --- TES (unchanged) ---
        if (
            tes_info == "correct"
            or tes_info == "potential template switching"
            or tes_info == "out of place"
        ):
            read, tes_found = read_ends_indexed(
                tes_index, read, contig, strand, tes_pos, args, "tes"
            )

        # --- INTRON ---
        try:
            intron_tag = read.get_tag("in")
        except KeyError:
            intron_tag = None

        intron_found = False
        if intron_tag:
            introns = literal_eval(intron_tag)
            read, intron_found = read_introns_fast(
                read, introns, args, contig, strand, intron_index
            )

        # GAP-like artefacts tagged as "ga"
        try:
            if read.get_tag("ga"):
                gap = True
        except KeyError:
            gap = False

        if tss_found and tes_found and not gap:
            leftend = read.get_tag("ts")
            rightend = read.get_tag("te")
            if intron_found:
                tr = (contig, strand, leftend, read.get_tag("in"), rightend)
            else:
                tr = (contig, strand, leftend, rightend)

            if tr in tr_dict:
                tr_dict[tr] += 1
            else:
                tr_dict[tr] = 1
            read.set_tag("tr", str(tr), "Z")

        outbam.write(read)

    bam.close()
    outbam.close()


def create_gff(tr_dict, tr_gff, args):
    """
    Pandas output to GFF3 + TSV.
    """
    c = 1
    PAR = ";Parent="
    ID = "ID="
    TRID = ";transcript_id="
    x = 1
    exon = "exon" + str(x)
    for key, value in tr_dict.items():
        if value >= args.mintr_count:
            tr = "tr" + str(c)
            if len(key) == 4:
                contig, strand, start, end = key
                l = len(tr_gff)
                if strand == "+":
                    tr_gff.loc[l] = [contig, "LoRTIA", "mRNA", start, end, value, strand, ".", ID + tr + TRID + tr]
                    tr_gff.loc[l+1] = [contig, "LoRTIA", "exon", start, end, value, strand, ".", ID + exon + PAR + tr + TRID + tr]
                else:
                    tr_gff.loc[l] = [contig, "LoRTIA", "mRNA", end, start, value, strand, ".", ID + tr + TRID + tr]
                    tr_gff.loc[l+1] = [contig, "LoRTIA", "exon", end, start, value, strand, ".", ID + exon + PAR + tr + TRID + tr]
            elif len(key) == 5:
                contig, strand, start, intron, end = key
                l = len(tr_gff)
                intron = literal_eval(intron)
                if strand == "+":
                    tr_gff.loc[l] = [contig, "LoRTIA", "mRNA", start, end, value, strand, ".", ID + tr + TRID + tr]
                    tr_gff.loc[l+1] = [contig, "LoRTIA", "exon", start, intron[0], value, strand, ".", ID + exon + PAR + tr + TRID + tr]
                    x += 1
                    exon = "id" + str(x)
                    if len(intron) > 2:
                        for i in range(1, len(intron)-1)[::2]:
                            l = len(tr_gff)
                            tr_gff.loc[l] = [contig, "LoRTIA", "exon", intron[i], intron[i+1], value, strand, ".", ID + exon + PAR + tr + TRID + tr]
                            x += 1
                            exon = "id" + str(x)
                    l = len(tr_gff)
                    tr_gff.loc[l] = [contig, "LoRTIA", "exon", intron[-1], end, value, strand, ".", ID + exon + PAR + tr + TRID + tr]
                    x += 1
                    exon = "id" + str(x)
                else:
                    tr_gff.loc[l] = [contig, "LoRTIA", "mRNA", end, start, value, strand, ".", ID + tr + TRID + tr]
                    tr_gff.loc[l+1] = [contig, "LoRTIA", "exon", end, intron[0], value, strand, ".", ID + exon + PAR + tr + TRID + tr]
                    x += 1
                    exon = "id" + str(x)
                    if len(intron) > 2:
                        for i in range(1, len(intron)-1)[::2]:
                            l = len(tr_gff)
                            tr_gff.loc[l] = [contig, "LoRTIA", "exon", intron[i], intron[i+1], value, strand, ".", ID + exon + PAR + tr + TRID + tr]
                            x += 1
                            exon = "id" + str(x)
                    l = len(tr_gff)
                    tr_gff.loc[l] = [contig, "LoRTIA", "exon", intron[-1], start, value, strand, ".", ID + exon + PAR + tr + TRID + tr]
                    x += 1
                    exon = "id" + str(x)

            tr_dict[key] = (tr, value)
            c += 1
        else:
            tr_dict[key] = ("not_qualified", value)

    tr_gff = tr_gff.sort_values(by=["contig", "start", "end"])
    tr_gff.to_csv(args.output_prefix + ".gff3", sep="\t", index=False, header=False)
    tr_tsv = pd.DataFrame.from_dict(tr_dict, orient="index")
    tr_tsv.to_csv(args.output_prefix + ".tsv", sep="\t", header=False)


def annotate_tr(args):
    print("Annotating transcripts...")
    cols = ["contig", "source", "feature", "start", "end", "score", "strand", "frame", "info"]

    # Load GFFs with Polars
    tss_gff = pl.read_csv(
        args.gff_prefix + "_tss.gff3",
        separator="\t",
        has_header=False,
        new_columns=cols,
    )
    tes_gff = pl.read_csv(
        args.gff_prefix + "_tes.gff3",
        separator="\t",
        has_header=False,
        new_columns=cols,
    )
    intron_gff = pl.read_csv(
        args.gff_prefix + "_intron.gff3",
        separator="\t",
        has_header=False,
        new_columns=cols,
    )

    # Build fast indexes
    tss_index = build_feature_index(tss_gff)
    tes_index = build_feature_index(tes_gff)
    intron_index = build_intron_index(intron_gff)

    bam = pysam.AlignmentFile(args.in_bam, "rb")
    tr_dict = {}
    outbam = pysam.AlignmentFile(args.output_prefix + ".bam", "wb", template=bam)

    bam_iterator(bam, tr_dict, tss_index, tes_index, intron_index, outbam, args)
    pysam.index(args.output_prefix + ".bam")

    tr_gff = pd.DataFrame(columns=cols)
    print("Generating gff...")
    create_gff(tr_dict, tr_gff, args)


def main():
    args = parsing()
    annotate_tr(args)


def parsing():
    parser = ArgumentParser(description="Annotate transcripts using accepted TSS/TES/intron features and BAM tags.")
    parser.add_argument("in_bam", help="Input .bam with tags (ts/te/in etc).", metavar="input_file")
    parser.add_argument("gff_prefix", help="Prefix for *_tss.gff3 *_tes.gff3 *_intron.gff3", metavar="gff_prefix")
    parser.add_argument("output_prefix", help="Output prefix for .bam and .gff3", metavar="output_prefix")

    parser.add_argument("--wobble", dest="wobble", type=int, default=10,
                        help="Window for feature matching (+/- wobble). Default 10.")

    parser.add_argument("--gap", dest="gap", type=int, default=1,
                        help="Largest allowed gap (unused in strict intron mode). Default 1.")

    parser.add_argument("--mintr_count", dest="mintr_count", type=int, default=1,
                        help="Minimum reads supporting a transcript to report it. Default 1.")

    # ✅ dRNA mode + shift
    parser.add_argument("--drna", action="store_true",
                        help="Enable dRNA-specific TSS handling: always attempt TSS, "
                             "try shifted match first then fallback to raw.")
    parser.add_argument("--drna-tss-shift", dest="drna_tss_shift", type=int, default=30,
                        help="dRNA TSS shift in nt (strand-aware) used in shifted match. Default 30. Set 0 to disable.")

    return parser.parse_args()


if __name__ == "__main__":
    main()