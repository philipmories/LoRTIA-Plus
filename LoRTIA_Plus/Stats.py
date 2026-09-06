#!/usr/bin/env python3

import pandas as pd
import os
import sys
from scipy.stats import poisson
from scipy.spatial import cKDTree
from argparse import ArgumentParser
from ast import literal_eval
import pysam
from bisect import bisect_left, bisect_right
from collections import deque



###############################################################################
# Coverage cache (fast): build per-contig dicts by streaming samtools depth once
###############################################################################

def _build_cov_by_contig(coverage_file, contigs_needed):
    """
    Stream samtools depth TSV (contig\\tpos\\tdepth) ONCE and keep only contigs in contigs_needed.
    Returns: {contig: {pos:int -> depth:int}}
    """
    need = set(contigs_needed)
    cov_by = {c: {} for c in need}
    with open(coverage_file, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            contig = parts[0]
            if contig not in need:
                continue
            try:
                pos = int(parts[1])
                depth = int(parts[2])
            except ValueError:
                continue
            cov_by[contig][pos] = depth
    return cov_by

###############################################################################
###                           Common functions                              ###
###############################################################################

def get_cov(position, cov_dict):
    """
    Calculates coverage based on a cov_dict
    """
    if position in cov_dict:
        return cov_dict.get(position)
    else:
        return 0

def coverage(pos_list, args, contig):
    """
    Calculates average coverages from a given distance for a list of positions.
    Uses args._cov_by_contig[contig] built once per Stats() call.
    """
    cov_dict = {}
    if hasattr(args, "_cov_by_contig"):
        cov_dict = args._cov_by_contig.get(contig, {})

    coverages = []
    for pos in pos_list:
        to_avg = []
        # original window logic kept, but cov_dict lookup is O(1)
        if args.distance > args.cov_sample:
            for position in range(pos + args.distance - args.cov_sample,
                                  pos + args.distance):
                to_avg.append(get_cov(position, cov_dict))
        elif args.distance < args.cov_sample:
            for position in range(pos + args.distance,
                                  pos + args.distance - args.cov_sample):
                to_avg.append(get_cov(position, cov_dict))
        else:
            to_avg.append(get_cov(pos + args.distance, cov_dict))
        coverages.append(sum(to_avg) / len(to_avg) if to_avg else 0)
    return coverages

def check_if_qualified(df, minimum, ratio):
    """
    Checks whether the feature position satisfies the minimum count and minimum
    ratio of coverage requirements.
    """
    qual_list = []
    for index, row in df.iterrows():
        is_qual = (row["count"] >= minimum
                   and row["ratio"] >= ratio
                   and row["is_picked"])
        qual_list.append(is_qual)
    return qual_list

def find_features(args):
    """
    Reads in a dataframe from csv, chops it and processes it on contigs.
    """
    df = pd.read_csv(args.feature_file, sep = "\t", names = ["pos", "count"])
    df["pos"] = df["pos"].apply(literal_eval)
    df[["contig", "pos"]] = df["pos"].apply(pd.Series)
    contig_set = list(set(df.contig))
    # Build per-contig coverage dicts once (FAST)
    args._cov_by_contig = _build_cov_by_contig(args.coverage_file, contig_set)

    # Introns also need the reference sequence. Open the indexed FASTA once
    # for the complete Stats() call instead of rescanning the whole FASTA
    # separately for every contig.
    fasta_handle = None
    if args.feature == "in":
        try:
            fasta_handle = pysam.FastaFile(args.reference)
        except (OSError, ValueError):
            # Create the standard .fai index if it is missing, then retry.
            try:
                pysam.faidx(args.reference)
                fasta_handle = pysam.FastaFile(args.reference)
            except Exception as exc:
                raise RuntimeError(
                    "Could not open or index the reference FASTA: "
                    f"{args.reference}"
                ) from exc
        args._fasta_handle = fasta_handle

    new_df = ()
    missing_reference_contigs = []
    try:
        for contig in contig_set:
            current_df = df.loc[df.contig == contig].copy()
            if args.feature == "in":
                # Human BAM headers frequently contain ALT/patch/decoy contigs
                # that are absent from a primary-assembly FASTA. Intron
                # sequence validation cannot be performed for such contigs,
                # so skip them instead of aborting the complete run.
                if contig not in fasta_handle.references:
                    missing_reference_contigs.append(
                        (str(contig), len(current_df))
                    )
                    continue
                current_df = contig_introns(current_df, args, contig)
            else:
                current_df = contig_ends(current_df, args, contig)
            if len(new_df) != 0:
                new_df = pd.concat([new_df, current_df])
            else:
                new_df = current_df
    finally:
        if fasta_handle is not None:
            fasta_handle.close()
            del args._fasta_handle

    if missing_reference_contigs:
        if args.feature_file.endswith("_in.tsv"):
            log_file = args.feature_file[:-7] + "_missing_reference_contigs.log"
        else:
            log_file = args.feature_file + ".missing_reference_contigs.log"

        with open(log_file, "w") as log:
            log.write("contig\tintron_candidates_skipped\n")
            for missing_contig, skipped_count in sorted(
                missing_reference_contigs
            ):
                log.write(f"{missing_contig}\t{skipped_count}\n")

        skipped_total = sum(
            skipped_count
            for _, skipped_count in missing_reference_contigs
        )
        print(
            "WARNING: Skipped "
            f"{skipped_total} intron candidates on "
            f"{len(missing_reference_contigs)} contig(s) absent from the "
            "reference FASTA. "
            f"Details: {log_file}"
        )
    if args.feature[1] == "3":
        feat = "_tes"
    elif args.feature[1] == "5":
        feat = "_tss"
    else:
        feat = "tron"
    new_df.to_csv(args.feature_file.replace(".tsv", "{}.tsv".format(feat)),
                  index=False,
                  sep="\t")

def Stats(args):
    """
    Sets argument types and runs stat functions for features.
    """
    if not args.feature:
        args.feature = args.feature_file[-6:-4]
    if args.feature == "r5" or args.feature == "r3":
        multiplier = -1
    else:
        multiplier = 1
    print("Calculating {} feature statistics...".format(args.feature))
    args.distance = abs(args.distance) * multiplier
    args.cov_sample = abs(args.cov_sample) * multiplier
    if os.stat(args.feature_file).st_size == 0:
        print("Feature file {} is empty. There is nothing to do here.".format(
              args.feature_file))
    else:
        find_features(args)


###############################################################################
###                             End functions                               ###
###############################################################################

def pick_from_greatests(dictionary, wobble):
    """
    Picks the left- or rightmost positions of the greatests list in a window
    determined by the wobble size. Whether the left or the rightmost positions
    are desired can be set by the user, and the list is ordered accordingly.
    """
    previous = -100
    is_picked_list = []
    for pos, is_greatest in dictionary.items():
        is_picked = False
        if is_greatest:
            if previous not in range(pos - wobble, pos + wobble + 1):
                is_picked = True
            previous = pos
        is_picked_list.append(is_picked)
    return is_picked_list

def check_if_greatest(tuples, wobble):
    """
    Finds the feature position with the highest read support in a window
    determined by the wobble size.

    This preserves the original LoRTIA rule exactly: a feature is marked
    greatest if no feature within [pos - wobble, pos + wobble] has a strictly
    larger count. Equal-count ties are therefore all retained as greatests.

    The old implementation scanned every tuple for every feature (O(N^2)).
    Here, a position-sorted copy and a monotonic deque maintain the maximum
    count in the moving genomic window in O(N) after sorting. Results are
    returned in the original tuple order.
    """
    if not tuples:
        return []

    # Preserve the old behaviour for an invalid negative wobble: Python's
    # range(pos - wobble, pos + wobble + 1) is empty, so every feature was
    # effectively considered a greatest.
    if wobble < 0:
        return [True] * len(tuples)

    indexed = sorted(
        enumerate(tuples),
        key=lambda item: item[1][0]
    )
    positions = [int(item[1][0]) for item in indexed]
    counts = [item[1][1] for item in indexed]

    result = [True] * len(tuples)
    max_deque = deque()
    right = 0
    n = len(indexed)

    for original_index, (pos, count) in indexed:
        pos = int(pos)
        upper = pos + wobble

        # Add every point that has entered the right edge of the current
        # symmetric window. Keep counts in decreasing order in the deque.
        while right < n and positions[right] <= upper:
            while max_deque and counts[max_deque[-1]] <= counts[right]:
                max_deque.pop()
            max_deque.append(right)
            right += 1

        lower = pos - wobble

        # Remove points that have left the left edge of the window.
        while max_deque and positions[max_deque[0]] < lower:
            max_deque.popleft()

        max_count = counts[max_deque[0]] if max_deque else count
        result[original_index] = count >= max_count

    return result

def count_average(tuples, window):
    """
    Counts the average values needed for the Poisson and the Pólya-Aeppli
    distributions.

    This preserves the original LoRTIA definition exactly: for each feature,
    sum the counts of every tuple whose position lies within
    [pos - window, pos + window], then divide by the full window width
    (2 * window + 1), including genomic positions with zero observed counts.

    The old implementation scanned all tuples for every feature (O(N^2)).
    Here, an indexed sorted copy plus prefix sums allows the same window sums
    to be obtained in O(N log N), while returning values in the original
    tuple order (important because r3/r5 are sorted descending upstream).
    """
    if not tuples:
        return []

    indexed = sorted(
        enumerate(tuples),
        key=lambda item: item[1][0]
    )
    positions = [int(item[1][0]) for item in indexed]

    prefix = [0]
    for _, (_, count) in indexed:
        prefix.append(prefix[-1] + count)

    denominator = 2 * window + 1
    result = [0.0] * len(tuples)

    for sorted_index, (original_index, (pos, count)) in enumerate(indexed):
        pos = int(pos)
        left = bisect_left(positions, pos - window)
        right = bisect_right(positions, pos + window)
        window_sum = prefix[right] - prefix[left]
        result[original_index] = window_sum / denominator

    return result

def _choose_jitter_peak(position, peak_indices, peak_positions, peak_raw_counts,
                        jitter, prefer_left=True):
    """
    Assign one genomic position to at most one already-selected legacy peak.

    Candidate peaks must lie within +/- ``jitter`` nt. Assignment priority is:
      1. smallest genomic distance to the peak;
      2. larger RAW peak support;
      3. legacy directional preference (leftmost for l5/l3, rightmost for r5/r3).

    The function never creates a new peak; it only chooses among peaks that
    were already retained by the unmodified legacy wobble logic.
    """
    if not peak_indices or jitter <= 0:
        return None

    left = bisect_left(peak_positions, int(position) - jitter)
    right = bisect_right(peak_positions, int(position) + jitter)
    if left >= right:
        return None

    candidates = peak_indices[left:right]
    return min(
        candidates,
        key=lambda idx: (
            abs(int(position) - int(peak_raw_counts[idx][0])),
            -float(peak_raw_counts[idx][1]),
            int(peak_raw_counts[idx][0]) if prefer_left
            else -int(peak_raw_counts[idx][0]),
        ),
    )


def _aggregate_support_onto_legacy_peaks(df, jitter, prefer_left=True):
    """
    Add mapper-jittered raw endpoint support to FIXED legacy wobble peaks.

    ``df`` must already contain the legacy ``is_greatest`` and ``is_picked``
    columns computed from the exact-coordinate raw endpoint counts. Therefore
    this function cannot create, remove, or move a peak. It changes only the
    support attached to rows with ``is_picked=True``.

    Raw endpoint coordinates that fall within the jitter radius of more than
    one legacy peak are assigned once using nearest-peak / raw-support /
    directional tie-breaking. Genomic positions in overlapping jitter windows
    are partitioned with the same rule so that the Poisson expectation can be
    scaled to the effective number of nucleotide positions represented by each
    aggregated support value without double-counting overlapping windows.
    """
    out = df.copy()
    out["raw_count"] = pd.to_numeric(out["count"], errors="coerce").fillna(0)
    out["jitter_count"] = out["raw_count"]
    out["jitter_width"] = 1
    out["jitter_members"] = 1

    if out.empty or jitter <= 0:
        return out

    peak_indices = list(out.index[out["is_picked"].astype(bool)])
    if not peak_indices:
        return out

    # Assignment lookup must be in ascending genomic order irrespective of
    # l5/l3 versus r5/r3 display order.
    peak_indices = sorted(peak_indices, key=lambda idx: int(out.at[idx, "pos"]))
    peak_positions = [int(out.at[idx, "pos"]) for idx in peak_indices]
    peak_raw_counts = {
        idx: (int(out.at[idx, "pos"]), float(out.at[idx, "raw_count"]))
        for idx in peak_indices
    }

    support = {idx: 0.0 for idx in peak_indices}
    member_count = {idx: 0 for idx in peak_indices}

    # Each observed exact endpoint coordinate contributes to at most one peak.
    for idx, row in out.iterrows():
        chosen = _choose_jitter_peak(
            int(row["pos"]),
            peak_indices,
            peak_positions,
            peak_raw_counts,
            jitter,
            prefer_left=prefer_left,
        )
        if chosen is None:
            continue
        support[chosen] += float(row["raw_count"])
        member_count[chosen] += 1

    # Partition the full nucleotide-level jitter windows by the same rule.
    # This gives the effective number of positions over which support was
    # pooled and therefore the width used to scale the raw per-base Poisson
    # background. Overlapping peak windows are never double-counted.
    width = {idx: 0 for idx in peak_indices}
    # Iterate only the UNION of local jitter windows. Never traverse the full
    # genomic span between distant peaks (which could be hundreds of Mb on a
    # human chromosome). Endpoint coordinates are 1-based, so clip the lower
    # edge at 1.
    jitter_positions = set()
    for peak_pos in peak_positions:
        jitter_positions.update(
            range(max(1, peak_pos - jitter), peak_pos + jitter + 1)
        )
    for pos in jitter_positions:
        chosen = _choose_jitter_peak(
            pos,
            peak_indices,
            peak_positions,
            peak_raw_counts,
            jitter,
            prefer_left=prefer_left,
        )
        if chosen is not None:
            width[chosen] += 1

    for idx in peak_indices:
        # A picked peak always assigns its own exact coordinate to itself, so
        # support should be >0. Keep a defensive fallback to raw_count.
        agg = support[idx]
        if agg <= 0:
            agg = float(out.at[idx, "raw_count"])
        out.at[idx, "count"] = agg
        out.at[idx, "jitter_count"] = agg
        out.at[idx, "jitter_width"] = max(1, int(width[idx]))
        out.at[idx, "jitter_members"] = max(1, int(member_count[idx]))

    return out


def _contig_ends_with_jitter(df, args, contig, jitter):
    """
    Endpoint-jitter mode with LEGACY-WOBBLE-FIRST peak selection.

    Processing order:
      raw exact-coordinate endpoints
        -> unchanged legacy local-maximum / wobble peak selection
        -> support aggregation onto those fixed peaks only
        -> minimum-count / coverage-ratio qualification
        -> Poisson using the unchanged legacy local background.

    Thus jitter can rescue an already-existing legacy peak by restoring support
    scattered over nearby mapper coordinates, but it cannot create a new peak,
    delete a legacy peak, or shift a peak coordinate.
    """
    df = df.copy()
    df["pos"] = df["pos"].astype(int)
    prefer_left = args.feature in ("l5", "l3")
    if prefer_left:
        df = df.sort_values(by="pos")
    else:
        df = df.sort_values(by="pos", ascending=False)

    # ------------------------------------------------------------------
    # 1) EXACT legacy peak selection on raw counts. Do not alter this.
    # ------------------------------------------------------------------
    raw_pos_count = list(zip(df["pos"], df["count"]))
    raw_average_values = count_average(raw_pos_count, args.window)
    raw_average = dict(zip(df.index, raw_average_values))

    df["is_greatest"] = check_if_greatest(raw_pos_count, args.wobble)
    greatests_dict = dict(zip(df["pos"], df["is_greatest"]))
    df["is_picked"] = pick_from_greatests(greatests_dict, args.wobble)

    # ------------------------------------------------------------------
    # 2) Jitter changes SUPPORT ONLY on those fixed legacy peaks.
    # ------------------------------------------------------------------
    df = _aggregate_support_onto_legacy_peaks(
        df,
        jitter,
        prefer_left=prefer_left,
    )

    # Coverage is still evaluated at the fixed legacy peak coordinate.
    df["coverage"] = coverage(df["pos"], args, contig)

    # Keep the legacy LoRTIA local Poisson background unchanged. Endpoint
    # jitter is interpreted as coordinate uncertainty of support for the SAME
    # already-selected biological end feature, not as 2*jitter+1 independent
    # genomic tests. Therefore only the observed support is collapsed onto the
    # fixed legacy peak; the local expected value remains exactly the value
    # calculated from the raw endpoint landscape. This also guarantees that
    # enabling jitter cannot make the per-peak Poisson p-value worse solely
    # because the jitter radius is larger.
    df["raw_average"] = [raw_average[i] for i in df.index]
    df["average"] = df["raw_average"]
    picked_mask = df["is_picked"].astype(bool)

    df["ratio"] = df["count"] / df["coverage"]
    df["is_qualified"] = check_if_qualified(df, args.minimum, args.ratio)
    # Non-peaks remain visible in the audit TSV but can never become final
    # features merely because jitter is enabled.
    df.loc[~picked_mask, "is_qualified"] = False

    # Preserve LoRTIA's historical tail definition and local expectation
    # exactly; only the observed count can increase through jitter collapsing.
    df["poisp"] = 1 - poisson.cdf(df["count"], df["average"])
    return df

def contig_ends(df, args, contig):
    """
    Processes the dataframe for one contig looking for TSSs or TESs.

    By default (``--end_jitter 0``) this function executes the validated
    legacy LoRTIA logic below without endpoint-support aggregation.
    """
    end_jitter = int(getattr(args, "end_jitter", 0) or 0)
    if end_jitter > 0:
        return _contig_ends_with_jitter(df, args, contig, end_jitter)

    # ---- validated legacy path: do not alter ----
    df["pos"] = df["pos"].astype(int)
    if args.feature == "l5" or args.feature == "l3":
    # This makes sure that the leftmost position is taken for each left feature
        df = df.sort_values(by = "pos")
    else:
        df = df.sort_values(by = "pos", ascending = False)
    df["coverage"] = coverage(df["pos"], args, contig)
    pos_count = list(zip(df["pos"], df["count"]))
    df["average"] = count_average(pos_count, args.window)
    df["is_greatest"] = check_if_greatest(pos_count, args.wobble)
    greatests_dict = dict(zip(df["pos"], df["is_greatest"]))
    df["is_picked"] = pick_from_greatests(greatests_dict, args.wobble)
    df["ratio"] = df["count"] / df["coverage"]
    df["is_qualified"] = check_if_qualified(df, args.minimum, args.ratio)
    # Preserve the original LoRTIA Poisson significance definition exactly.
    # Legacy LoRTIA uses P(X > observed_count) = 1 - CDF(observed_count, mu).
    # Do not change this to sf(count - 1): that would alter historical TSS/TES
    # calls, especially for low-count end features.
    df["poisp"] = 1 - poisson.cdf(df["count"], df["average"])
    return df


###############################################################################
###                          Intronic functions                             ###
###############################################################################

def intron_picker(tuples, args):
    """
    Checks whether there is a more frequent intron near the intron and if there
    is, it only picks the intron if it is more abundant than the set 'rare
    intron' limit.

    The original implementation compared every intron with every other intron
    (O(N^2)). Here the (left, right) splice-junction coordinates are indexed
    with a cKDTree and queried with Chebyshev distance (p=inf), which is exactly
    the old rectangular rule:

        abs(left2 - left) <= intron_wobble
        AND
        abs(right2 - right) <= intron_wobble

    The downstream abundance rule is intentionally unchanged.
    """
    if not tuples:
        return []

    wobble = args.intron_wobble

    # Preserve the old behaviour for a negative wobble: both Python ranges
    # are empty, so no neighbours are collected and every intron is retained.
    if wobble < 0:
        return [True] * len(tuples)

    points = [(float(left), float(right)) for left, right, count in tuples]
    counts = [count for left, right, count in tuples]

    tree = cKDTree(points)
    neighbourhoods = tree.query_ball_point(
        points,
        r=wobble,
        p=float("inf"),
    )

    is_picked_list = []
    for count, neighbours in zip(counts, neighbourhoods):
        if len(neighbours) > 1:
            max_count = max(counts[index] for index in neighbours)
            is_picked = count / args.rare_intron >= max_count
        else:
            is_picked = True
        is_picked_list.append(is_picked)

    return is_picked_list
        
def check_consensus(left2, right2):
    """
    Finds consensus splice junctions and sets the strand of the intron
    accordingly.
    """
    consensus_list = []
    strand_list = []
    for i in range(len(left2)):
        if left2[i].lower() == "gt" and right2[i].lower() == "ag":
            consensus_list.append("GT/AG")
            strand_list.append("+")
        elif left2[i].lower() == "gc" and right2[i].lower() == "ag":
            consensus_list.append("GC/AG")
            strand_list.append("+")
        elif left2[i].lower() == "at" and right2[i].lower() == "ac":
            consensus_list.append("AT/AC")
            strand_list.append("+")
        elif left2[i].lower() == "ct" and right2[i].lower() == "ac":
            consensus_list.append("GT/AG")
            strand_list.append("-")
        elif left2[i].lower() == "ct" and right2[i].lower() == "gc":
            consensus_list.append("GC/AG")        
            strand_list.append("-")
        elif left2[i].lower() == "gt" and right2[i].lower() == "at":
            consensus_list.append("AT/AC")        
            strand_list.append("-")
        else:
            consensus_list.append("None")        
            strand_list.append(".")
    return consensus_list, strand_list


def get_score(scores, args):
    """
    Gets the alignment scores from the aligner and calculates the maximum
    alignment score between the two exon ends that could have lead to template
    switching.
    """
    pos = args.shs_for_ts + 1
    rscores = []
    score = 0
    while pos <= 2 * args.shs_for_ts:
        score += scores[pos - 1]
        rscores.append(score)
        pos += 1
    rmax = max(rscores)
    if rmax < 0:
        rmax = 0
    pos = args.shs_for_ts
    lscores = []
    score = 0
    while pos >= 1:
        score += scores[pos - 1]
        lscores.append(score)
        pos += -1
    lmax = max(lscores)
    if lmax < 0:
        lmax = 0
    return rmax + lmax

def align(lseq, rseq, args):
    """
    This is a special local aligner that does not allow gaps nor shifts. It
    reports the highest alignment score that borders or involves the centre.
    """
    scores = []
    for i in range(len(lseq)):
        if lseq[i] == rseq[i]:
            scores.append(args.match_score)
        else:
            scores.append(args.mismatch_score)
    return get_score(scores, args)

def _fetch_like_python_slice(fasta, contig, start, end):
    """
    Fetches a FASTA interval with the same boundary behaviour as the old
    SeqIO/Bio.Seq expression: sequence[start:end].

    pysam.fetch() uses 0-based, end-exclusive coordinates, but unlike normal
    Python slicing it does not accept all out-of-range/negative coordinates.
    slice.indices() normalizes them first, preserving the old behaviour even
    for features very close to contig boundaries.
    """
    contig_length = fasta.get_reference_length(contig)
    norm_start, norm_end, _ = slice(int(start), int(end)).indices(contig_length)
    if norm_start >= norm_end:
        return ""
    return fasta.fetch(contig, norm_start, norm_end)


def intron_seq(df, args, contig):
    """
    Gets exon-intron border sequences and sends them to the aligner.

    The reference FASTA is indexed and opened once per intron Stats() run
    (in find_features), so this function performs direct random access instead
    of rescanning the complete FASTA for every contig.
    """
    ts = args.shs_for_ts
    fasta = getattr(args, "_fasta_handle", None)
    if fasta is None:
        raise RuntimeError("Reference FASTA handle was not initialized.")

    if contig not in fasta.references:
        raise RuntimeError(
            "Internal reference-contig routing error: "
            f"'{contig}' was passed to intron_seq() although it is absent "
            "from the supplied reference FASTA."
        )

    is_ts_list = []
    lseq_list = []
    rseq_list = []
    l2_list = []
    r2_list = []
    ts_list = []

    for index, row in df.iterrows():
        left = int(row["left"])
        right = int(row["right"]) - 1

        leftseq = _fetch_like_python_slice(
            fasta, contig, left - ts, left + ts
        )
        rightseq = _fetch_like_python_slice(
            fasta, contig, right - ts, right + ts
        )

        ts_score = align(leftseq, rightseq, args)
        is_ts = ts_score >= args.match_score * ts
        is_ts_list.append(is_ts)
        ts_list.append(ts_score)
        lseq_list.append(leftseq)
        rseq_list.append(rightseq)
        l2_list.append(leftseq[ts:ts + 2])
        r2_list.append(rightseq[ts - 2:ts])

    df["is_potential_ts"] = is_ts_list
    df["leftseq"] = lseq_list
    df["rightseq"] = rseq_list
    df["left2"] = l2_list
    df["right2"] = r2_list
    df["ts_score"] = ts_list
    return df

def contig_introns(df, args, contig):
    """
    Processes the dataframe for one contig looking for introns.
    """
    df[["left", "right"]] = df["pos"].apply(pd.Series)
    df["left"] = df["left"].astype(int)
    df["right"] = df["right"].astype(int)
    df["rcov"] = coverage(df["right"], args, contig)
    args.distance = abs(args.distance) * -1
    args.cov_sample = abs(args.cov_sample) * -1
    df["lcov"] = coverage(df["left"], args, contig)
    df["coverage"] = (df["lcov"] + df["rcov"]) / 2
    df["ratio"] = df["count"] / df["coverage"]
    df = intron_seq(df, args, contig)
    left_right_count = list(zip(df["left"], df["right"], df["count"]))
    df["is_picked"] = intron_picker(left_right_count, args)
    df["is_qualified"] = check_if_qualified(df, args.minimum, args.ratio)
    df["consensus"], df["strand"] = check_consensus(list(df["left2"]),
                                                    list(df["right2"]))
    args.distance = args.distance * -1
    args.cov_sample = args.cov_sample * -1
    return df


###############################################################################
###                             Main function                               ###
###############################################################################

def main():
    args = parsing()
    Stats(args)


def parsing():
    parser = ArgumentParser(description="This is the second module of \
                            LoRTIA, a Long-read RNA-Seq Transcript Isofom \
                            Annotator. This module calculates the statistics \
                            of transcript features.")
    parser.add_argument("coverage_file",
                        help="The tsv file which contains the coverages.\
                        The tsv file should contain 3 columns: contig, \
                        position and coverage.",
                        metavar="coverage_file")
    parser.add_argument("feature_file",
                        help="A tab-separated values file containing feature\
                        statistics produced by the Samprocessor.",
                        metavar="feature_file")
    parser.add_argument("-w", "--window",
                        dest="window",
                        help="The window that is examined when calculating \
                        the Poisson distribution. Setting low values finds \
                        false positives in a noisy data, while setting high \
                        values leads to false negatives due to the different \
                        transcriptional activity of different genomic regions.\
                        The default value is 50, which translates to a 101 nt\
                        bin (examined nucleotide +/- 50 nucleotides).",
                        type=int,
                        default=50,
                        metavar="[integer]")
    parser.add_argument("-r", "--reference",
                        dest="reference",
                        help="Reference FASTA used for intron splice-site and "
                             "junction-boundary sequence evaluation. Required "
                             "for standalone intron analysis.",
                        default=None,
                        metavar="[reference_fasta]")
    parser.add_argument("-m", "--minimum",
                        dest="minimum", 
                        help="The minimal number of reads for the feature to\
                        be accepted.",
                        type=int,
                        default=2,
                        metavar="[integer]")
    parser.add_argument("-f", "--feature",
                        dest="feature",
                        help="The feature that is examined. Options are \
                        'r5' for reverse strand 5' ends, 'l3' for \
                        reverse strand 3' ends, 'l5' for forward strand 5'\
                        ends, 'r3' for forward strand 3' ends and 'in' for \
                        introns. By default the tsv file's last two characters\
                        before the .tsv extension are considered.",
                        default=False,
                        metavar="[string]")
    parser.add_argument("-b", "--wobble",
                        dest="wobble",
                        help="The window, in which only one of each feature \
                        is expected, and locations with lesser support are \
                        considered to be derivatives of the major. The default\
                        value is 10, which means that only one feature of a \
                        kind can be described in a 21 nt bin (location +/-10 \
                        nt). This only applies to TSSs and TESs.",
                        type=int,
                        default=10,
                        metavar="[integer]")
    parser.add_argument("--end_jitter", "--endpoint_jitter",
                        dest="end_jitter",
                        help="Optional TSS/TES endpoint-support aggregation window. "
                        "Exact end coordinates within +/- N nt of a selected "
                        "representative are summed before feature qualification. "
                        "This is independent of --wobble and is disabled by "
                        "default (0) to preserve legacy LoRTIA output.",
                        type=int,
                        default=0,
                        metavar="[integer]")
    parser.add_argument("-i", "--intron_wobble",
                        dest="intron_wobble",
                        help="This option is only important for error-prone \
                        reads. Sequencing errors can disrupt the mapping of \
                        introns. Rare splice juntions can be detected in the \
                        close vicinity of more frequently utilized splice\
                        junctions. The rare splice junctions are likely to \
                        be results of sequencing errors of the more frequent \
                        version. This option regulates the window in which a \
                        rare intron will be considered to have stemmed from a \
                        sequencing error. The default value is 15 nt. That \
                        means that the rare introns which are no further than \
                        15 nt away from more frequent introns, will be \
                        considered to be sequencing errors.",
                        type=int,
                        default=15,
                        metavar="[integer]")
    parser.add_argument("--rare_intron",
                        dest="rare_intron",
                        help="This option is only important for error-prone \
                        reads. Sequencing errors can disrupt the mapping of \
                        introns. Rare splice juntions can be detected in the \
                        close vicinity of more frequently utilized splice\
                        junctions. The rare splice junctions are likely to \
                        be results of sequencing errors of the more frequent \
                        version. This option determines how much rarer an \
                        should be than the most frequent intron in its +/- \
                        'intron_wobble' vicinity, in orderd to be discarded as\
                        a sequencing error. The default value is 0.05.",
                        type=float,
                        default=0.05,
                        metavar="[float]")
    parser.add_argument("--match_score", 
                        dest="match_score",
                        help="The alignment scores for each match when \
                        searching for adapters. Penalty scores should be \
                        supplied as negative vaules. The default is: 2",
                        type=float,
                        metavar="[float]", 
                        default=2.0,
                        required=False)
    parser.add_argument("--mismatch_score", 
                        dest="mismatch_score",
                        help="The alignment scores for each mismatch when \
                        searching for adapters. Penalty scores should be \
                        supplied as negative vaules. The default is: -3",
                        type=float,
                        metavar="[float]", 
                        default=-3.0,
                        required=False)
    parser.add_argument("--shs_for_ts",
                        dest="shs_for_ts",
                        help="Half-window size used for intron-boundary short-"
                             "homology scoring between donor- and acceptor-side "
                             "reference sequence. The resulting score is retained "
                             "as junction-level template-switching QC information. "
                             "Default: 3 nt.",
                        type=int,
                        default=3,
                        metavar="[integer]")
    parser.add_argument("-t", "--ratio",
                        dest = "ratio",
                        help = "The minimal ratio of the coverage that a \
                        feature has to reach to be accepted. The default value\
                        is 0.001.",
                        type=float,
                        default=0.001,
                        metavar="[float]")
    parser.add_argument("-d", "--distance",
                        dest="distance",
                        help="The distance from the feature position where \
                        coverage should be calculated. The default value is \
                        15. A positive value should be given, the inward \
                        direction is calculated by the program automatically.",
                        type=int,
                        default=15,
                        metavar="[integer]")
    parser.add_argument("-s", "--cov_sample",
                        dest="cov_sample",
                        help="The number of nucleotides where the coverage \
                        should be averaged. This many consecutive nucleotides\
                        will be considered from the 'distance' towards the \
                        feature. Its absolute value has to be smaller than \
                        or equal to the value of 'distance'. The default value\
                        is 5.",
                        type=int,
                        default=5,
                        metavar="[integer]")
    args = parser.parse_args()
    if not args.feature:
        args.feature = args.feature_file[-6:-4]
    if args.feature == "in" and not args.reference:
        parser.error("--reference is required for standalone intron analysis.")
    return args
    

if __name__== "__main__":
    main()
