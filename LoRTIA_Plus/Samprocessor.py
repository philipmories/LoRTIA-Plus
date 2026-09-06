#!/usr/bin/env python3

###############################################################################
###                            SAMprocessor                                 ###
###                 This script handles the .sam files.                     ###
###                And oh boy, the way it handles them...                   ###
###############################################################################

import pysam
import operator
import os
import subprocess
import pandas as pd
from argparse import ArgumentParser
from Bio.Seq import Seq
from Bio import Align

###############################################################################
###                               ADAPTERS                                  ###
###############################################################################
#TeloPrime 5' TGGATTGATATGTAATACGACTCACTATAG
#PacBio 5' AGAGTACATGGG
#MinION 5' TTTCTGTTGGTGCTGATATTGCTGCCATTAGGCCGGG
#MinION 3' ACTTGCCTGTCGCTCTATCTTC


###############################################################################
###                 Bedtools and Samtools-like functions                    ###
###############################################################################

def input_sorter(in_file, prefix):
    """
    Calls the pysam sort function on the input.
    """
    print("Sorting {}...".format(in_file))
    pysam.sort("-n", "-O" "SAM", "-o", "{}_temp.sam".format(prefix), in_file)


# The code below was meant to run bedtools and get output line-by-line and
# only write out non0 coverage in order to buy space, but it runs incredibly
# slowly.
#
#def runner(bam, tsv, strand):
#    cmd = ["genomeCoverageBed", "-ibam {}".format(bam), "-d", "-split {}".format(strand)]
#    cmd = ["genomeCoverageBed -ibam {} -d -split {}".format(bam, strand)]
#    btools = Popen(cmd, stdout=PIPE, shell = True)
#    while True:
#        line = btools.stdout.readline().decode('utf-8')
#        if not line:
#            break
#        if line[-3:] != "\t0\n":
#            out_appender(line, tsv)
#    btools.stdout.close()

def _samtools_supports_threads():
    """
    Best-effort detection whether `samtools depth` supports -@ (threads).
    If detection fails, we return False to avoid breaking older samtools.
    """
    try:
        import subprocess
        p = subprocess.run(["samtools", "depth", "--help"], capture_output=True, text=True)
        return "-@" in (p.stdout + p.stderr)
    except Exception:
        return False


def coverage_counter(outbam, out_stranded_bam, threads=1):
    """
    Create coverage TSVs using **samtools depth** (fast).

    Outputs (contig\tpos\tdepth):
      - *allcov.tsv    from outbam
      - *pluscov.tsv   from out_stranded_bam, forward reads only (FLAG !16)
      - *minuscov.tsv  from out_stranded_bam, reverse reads only (FLAG 16)

    NOTE:
      - We intentionally do NOT output zero-depth positions to keep files small.
        Missing positions are treated as depth=0 by Stats.py.
      - Uses -d 0 (no max-depth cap).
    """
    import subprocess

    print("Counting coverage (samtools depth)...")
    alltsv = outbam.replace("sorted.bam", "allcov.tsv")
    mintsv = outbam.replace("sorted.bam", "minuscov.tsv")
    plustsv = outbam.replace("sorted.bam", "pluscov.tsv")

    thr = max(int(threads), 1)
    depth_thr = ""
    view_thr = ""
    if thr > 1 and _samtools_supports_threads():
        depth_thr = f" -@ {thr}"
    if thr > 1:
        # samtools view supports -@ in practically all modern versions
        view_thr = f" -@ {thr}"

    # All coverage
    cmd_all = f"samtools depth{depth_thr} -g 0x700 -d 0 {outbam} | awk '$3 != 0' > {alltsv}"
    subprocess.run(cmd_all, shell=True, check=True)

    # Plus strand (forward reads): exclude FLAG 16
    cmd_plus = (
        f"samtools view{view_thr} -b -F 16 {out_stranded_bam} | "
        f"samtools depth{depth_thr} -g 0x700 -d 0 - | awk '$3 != 0' > {plustsv}"
    )
    subprocess.run(cmd_plus, shell=True, check=True)

    # Minus strand (reverse reads): include FLAG 16
    cmd_minus = (
        f"samtools view{view_thr} -b -f 16 {out_stranded_bam} | "
        f"samtools depth{depth_thr} -g 0x700 -d 0 - | awk '$3 != 0' > {mintsv}"
    )
    subprocess.run(cmd_minus, shell=True, check=True)


def output_creator(outsam, threads=1):
    """
    Uses pysam/samtools to sort, index and filter the output BAMs. ``threads``
    affects only samtools I/O/compression operations.
    """
    print("Sorting and indexing output...")
    outbam = outsam.replace(".sam", "_sorted.bam")
    out_stranded_bam = outsam.replace("out.sam", "stranded_only.bam")
    thr = max(int(threads), 1)
    pysam.sort("-@", str(thr), "-o", outbam, outsam)
    pysam.index("-@", str(thr), outbam)
    pysam.view("-@", str(thr),
               "-h",
               "-b",
               "-F", "512",
               "-o",
               out_stranded_bam,
               outbam,
               catch_stdout=False)
    pysam.index("-@", str(thr), out_stranded_bam)
    coverage_counter(outbam, out_stranded_bam, threads=thr)


###############################################################################
###                           Cigar functions                               ###
###############################################################################

QUERY_CONSUMING_CIGAR = {0, 1, 4, 7, 8}
TERMINAL_CLIP_CIGAR = {4, 5}


def _full_query_length_from_cigar(cigar):
    """Full molecule length implied by CIGAR, including hard-clipped bases."""
    if not cigar:
        return 0
    return sum(length for op, length in cigar
               if op in QUERY_CONSUMING_CIGAR or op == 5)


def _terminal_clip_length(cigar, left=True):
    """Sum consecutive terminal S/H operations on one CIGAR end."""
    if not cigar:
        return 0
    events = cigar if left else reversed(cigar)
    total = 0
    for op, length in events:
        if op in TERMINAL_CLIP_CIGAR:
            total += length
        else:
            break
    return total


def _sequence_is_complete_for_cigar(read_seq, cigar):
    """True if SEQ contains the full molecule represented by S/H-aware CIGAR."""
    if read_seq is None:
        return False
    expected = _full_query_length_from_cigar(cigar)
    if expected == 0:
        return False
    return len(read_seq) == expected


def get_left_end(read_seq, cigar, args):
    """Return an adapter-search window around the left alignment boundary.

    Unlike the legacy implementation, this function does not require a soft
    clip. If a complete molecule sequence is available, hard-clipped sequence
    can also be searched because the alignment boundary is reconstructed from
    the CIGAR. If the record itself is hard-clipped and no complete sequence
    was recovered from another record with the same QNAME, ``None`` is
    returned rather than indexing absent bases.
    """
    if read_seq is None or not cigar:
        return None, 0, ""
    if any(op == 5 for op, _ in cigar) and not _sequence_is_complete_for_cigar(read_seq, cigar):
        return None, 0, ""

    boundary = _terminal_clip_length(cigar, left=True)
    if boundary == 0:
        return None, 0, ""
    boundary = min(max(boundary, 0), len(read_seq))
    checkleft = max(0, boundary - args.check_in_soft)
    checkright = min(len(read_seq), boundary + args.check_in_match)
    cis = boundary - checkleft
    leftend_seq = read_seq[checkleft:checkright]
    left_match = read_seq[boundary:min(len(read_seq), boundary + args.match_in_first)]
    return leftend_seq, cis, left_match


def get_right_end(read_seq, cigar, args):
    """Return an adapter-search window around the right alignment boundary.

    The boundary is reconstructed from terminal S/H operations and therefore
    can use hard-clipped sequence when a complete same-QNAME molecule sequence
    has been recovered.
    """
    if read_seq is None or not cigar:
        return None, 0, ""
    if any(op == 5 for op, _ in cigar) and not _sequence_is_complete_for_cigar(read_seq, cigar):
        return None, 0, ""

    right_clip = _terminal_clip_length(cigar, left=False)
    if right_clip == 0:
        return None, 0, ""
    boundary = len(read_seq) - right_clip
    boundary = min(max(boundary, 0), len(read_seq))
    checkleft = max(0, boundary - args.check_in_match)
    checkright = min(len(read_seq), boundary + args.check_in_soft)
    cis = checkright - boundary
    rightend_seq = read_seq[checkleft:checkright]
    right_match = read_seq[max(0, boundary - args.match_in_first):boundary]
    return rightend_seq, cis, right_match


def _group_forward_sequence(reads):
    """Recover one canonical full molecule sequence for a same-QNAME group.

    ``AlignedSegment.get_forward_sequence()`` restores sequencing orientation.
    We prefer a record whose SEQ length equals the full CIGAR-implied molecule
    length (i.e. no unrecoverable hard-clipped bases), then primary over
    supplementary/secondary records, then longer/higher-MAPQ records.
    """
    candidates = []
    for read in reads:
        if read.query_sequence is None:
            continue
        try:
            seq = read.get_forward_sequence()
        except Exception:
            seq = None
        if not seq:
            continue
        expected = _full_query_length_from_cigar(read.cigartuples)
        complete = (expected > 0 and len(seq) == expected)
        primary = not read.is_secondary and not read.is_supplementary
        candidates.append((complete, primary, len(seq), int(read.mapping_quality), str(seq)))
    if not candidates:
        return None
    candidates.sort(reverse=True)
    best = candidates[0]
    return best[-1] if best[0] else None


def _sequence_for_alignment_orientation(read, full_forward_seq=None):
    """Return complete sequence in the same orientation as this SAM record."""
    if full_forward_seq is None:
        return read.query_sequence
    if read.is_reverse:
        return str(Seq(full_forward_seq).reverse_complement())
    return str(full_forward_seq)


def _query_interval_original_orientation(read):
    """Return (qstart, qend, full_length) in original sequencing orientation."""
    cigar = read.cigartuples or []
    full_len = _full_query_length_from_cigar(cigar)
    left_clip = _terminal_clip_length(cigar, left=True)
    right_clip = _terminal_clip_length(cigar, left=False)
    if read.is_reverse:
        qstart = right_clip
        qend = full_len - left_clip
    else:
        qstart = left_clip
        qend = full_len - right_clip
    qstart = max(0, min(qstart, full_len))
    qend = max(qstart, min(qend, full_len))
    return qstart, qend, full_len


def _query_overlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def _alignment_priority(read):
    """Rank alternative alignments occupying the same query interval."""
    try:
        as_score = float(read.get_tag("AS"))
    except Exception:
        as_score = float("-inf")
    qstart, qend, _ = _query_interval_original_orientation(read)
    return (
        1 if not read.is_secondary else 0,
        1 if not read.is_supplementary else 0,
        int(read.mapping_quality),
        as_score,
        qend - qstart,
    )


def _interval_union_length(intervals):
    """Length of the union of half-open query intervals."""
    cleaned = sorted((int(a), int(b)) for a, b in intervals if int(b) > int(a))
    if not cleaned:
        return 0
    total = 0
    cur_start, cur_end = cleaned[0]
    for start, end in cleaned[1:]:
        if start <= cur_end:
            cur_end = max(cur_end, end)
        else:
            total += cur_end - cur_start
            cur_start, cur_end = start, end
    total += cur_end - cur_start
    return int(total)


def _select_split_alignments(reads, max_query_overlap=20):
    """Keep distinct query segments and retain ambiguity/coverage metadata.

    A genuine chimeric/split read is defined in query space: different linear
    alignments cover different portions of the same molecule. Strongly
    overlapping alternatives are collapsed exactly as before, but their
    existence is now retained as QC metadata instead of being forgotten.

    Returns
    -------
    (kept, qc)
        ``kept`` is the unchanged selected alignment chain. ``qc`` contains
        query coverage of the selected chain plus counts describing mapping
        ambiguity in the original same-QNAME alignment group. These QC values
        do not alter alignment selection or fusion calling.
    """
    ranked = sorted(reads, key=_alignment_priority, reverse=True)
    kept = []
    kept_intervals = []
    alternative_mapping_count = 0
    ambiguous_kept_indices = set()

    for read in ranked:
        qi = _query_interval_original_orientation(read)
        interval = (qi[0], qi[1])
        if interval[1] <= interval[0]:
            continue

        conflicting = [
            i for i, old in enumerate(kept_intervals)
            if _query_overlap(interval, old) > max_query_overlap
        ]
        if conflicting:
            alternative_mapping_count += 1
            ambiguous_kept_indices.update(conflicting)
            continue

        kept.append(read)
        kept_intervals.append(interval)

    kept.sort(key=lambda r: _query_interval_original_orientation(r)[:2])

    selected_intervals = [
        _query_interval_original_orientation(r)[:2] for r in kept
    ]
    qlens = [
        _query_interval_original_orientation(r)[2] for r in kept
        if _query_interval_original_orientation(r)[2] > 0
    ]
    qlen = max(qlens, default=0)
    covered = _interval_union_length(selected_intervals)
    query_coverage = (float(covered) / float(qlen)) if qlen > 0 else 0.0

    secondary_alignment_count = sum(1 for r in reads if r.is_secondary)
    supplementary_alignment_count = sum(1 for r in reads if r.is_supplementary)
    mapping_ambiguous = int(
        alternative_mapping_count > 0 or secondary_alignment_count > 0
    )

    qc = {
        "query_covered_bases": int(covered),
        "query_length": int(qlen),
        "query_coverage": float(query_coverage),
        "alternative_mapping_count": int(alternative_mapping_count),
        "ambiguous_segment_count": int(len(ambiguous_kept_indices)),
        "secondary_alignment_count": int(secondary_alignment_count),
        "supplementary_alignment_count": int(supplementary_alignment_count),
        "mapping_ambiguous": int(mapping_ambiguous),
    }
    return kept, qc



def _adapter_state_from_summary(summary, tag):
    try:
        return summary[tag][3]
    except Exception:
        return "missing"


def _tag_on_original_query_side(read, side, kind):
    """Adapter tag corresponding to one end of the original sequenced molecule.

    adapter_checker() works in the SAM-record/alignment orientation. Reverse
    alignments therefore swap which alignment side corresponds to the original
    query left/right side.
    """
    if side not in {"left", "right"} or kind not in {"5", "3"}:
        raise ValueError("Invalid side/kind")
    if side == "left":
        lr = "r" if read.is_reverse else "l"
    else:
        lr = "l" if read.is_reverse else "r"
    return lr + kind


def _adapter_molecule_direction(prepared, args):
    """Infer molecule direction from the two external LoRTIA adapter calls."""
    if len(prepared) < 2:
        return "ambiguous"
    left_read, left_summary = prepared[0]
    right_read, right_summary = prepared[-1]

    left3 = _adapter_state_from_summary(
        left_summary, _tag_on_original_query_side(left_read, "left", "3")
    )
    right3 = _adapter_state_from_summary(
        right_summary, _tag_on_original_query_side(right_read, "right", "3")
    )

    if getattr(args, "dRNA", False):
        left_ok = left3 == "correct"
        right_ok = right3 == "correct"
        if right_ok and not left_ok:
            return "5to3"
        if left_ok and not right_ok:
            return "3to5"
        return "ambiguous"

    left5 = _adapter_state_from_summary(
        left_summary, _tag_on_original_query_side(left_read, "left", "5")
    )
    right5 = _adapter_state_from_summary(
        right_summary, _tag_on_original_query_side(right_read, "right", "5")
    )

    forward = (left5 == "correct" and right3 == "correct")
    reverse = (left3 == "correct" and right5 == "correct")
    if forward and not reverse:
        return "5to3"
    if reverse and not forward:
        return "3to5"
    return "ambiguous"


def _ts_consensus(prepared):
    """Return minimap2 transcript-strand consensus over retained split segments.

    ``ts:A:+`` means the original query/read orientation is the same as the
    transcript 5'->3' orientation, whereas ``ts:A:-`` means it is opposite.
    Only the query-space-selected alignment segments are considered here.
    """
    values = []
    for read, _summary in prepared:
        try:
            value = str(read.get_tag("ts"))
        except Exception:
            continue
        if value in {"+", "-"}:
            values.append(value)

    if not values:
        return "missing", "ambiguous"
    unique = set(values)
    if unique == {"+"}:
        return "+", "5to3"
    if unique == {"-"}:
        return "-", "3to5"
    return "conflict", "ambiguous"


def _molecule_orientation(prepared, args):
    """Resolve molecule direction while retaining adapter evidence separately.

    Priority:
      1. consistent minimap2 ``ts`` consensus across retained query segments;
      2. external LoRTIA adapter orientation as a fallback;
      3. ambiguous when neither resolves the molecule.

    Adapter calls are *not* removed from validation.  ``adapter_direction`` is
    stored independently and is used downstream as terminal-QC evidence.
    """
    ts_state, ts_direction = _ts_consensus(prepared)
    adapter_direction = _adapter_molecule_direction(prepared, args)

    if ts_direction in {"5to3", "3to5"}:
        return {
            "direction": ts_direction,
            "source": "minimap2_ts",
            "ts_consensus": ts_state,
            "adapter_direction": adapter_direction,
        }

    if adapter_direction in {"5to3", "3to5"}:
        source = "adapter_fallback_no_ts" if ts_state == "missing" else "adapter_fallback_conflicting_ts"
        return {
            "direction": adapter_direction,
            "source": source,
            "ts_consensus": ts_state,
            "adapter_direction": adapter_direction,
        }

    source = "conflicting_ts" if ts_state == "conflict" else "unresolved"
    return {
        "direction": "ambiguous",
        "source": source,
        "ts_consensus": ts_state,
        "adapter_direction": adapter_direction,
    }


def _internal_adapter_present(prepared, args):
    """Detect adapter signatures at internal split boundaries.

    A correct internal library adapter is evidence for a concatenated library
    molecule rather than a biological fusion. In dRNA mode only the 3' adapter
    is informative because 5' ends are intentionally marked correct without
    adapter search.
    """
    if len(prepared) < 2:
        return False
    kinds = ("3",) if getattr(args, "dRNA", False) else ("5", "3")
    for (left_read, left_summary), (right_read, right_summary) in zip(prepared, prepared[1:]):
        for kind in kinds:
            left_tag = _tag_on_original_query_side(left_read, "right", kind)
            right_tag = _tag_on_original_query_side(right_read, "left", kind)
            if _adapter_state_from_summary(left_summary, left_tag) == "correct":
                return True
            if _adapter_state_from_summary(right_summary, right_tag) == "correct":
                return True
    return False

###############################################################################
###                             Adapter functions                           ###
###############################################################################

_ADAPTER_ALIGNER_CACHE = {}

def _get_adapter_aligner(args):
    """
    Return a configured Bio.Align.PairwiseAligner for adapter detection.

    The scoring is intentionally identical to the previous
    local-alignment call. Aligners are cached by scoring parameters so
    that a new PairwiseAligner object is not constructed for every read end.
    """
    key = (float(args.match_score),
           float(args.mismatch_score),
           float(args.gap_open_score),
           float(args.gap_extend_score))

    aligner = _ADAPTER_ALIGNER_CACHE.get(key)
    if aligner is None:
        aligner = Align.PairwiseAligner()
        aligner.mode = "local"
        aligner.match_score = key[0]
        aligner.mismatch_score = key[1]
        aligner.open_gap_score = key[2]
        aligner.extend_gap_score = key[3]
        _ADAPTER_ALIGNER_CACHE[key] = aligner
    return aligner

def adapter_aligner(sequence, adapter, args):
    """
    Align the adapter to ``sequence`` using Bio.Align.PairwiseAligner.

    This is a compatibility-preserving local-alignment replacement. The
    downstream LoRTIA adapter-classification code is intentionally left
    unchanged.  Therefore this function returns five-element tuples with the
    same fields that the existing code uses:

        (adapter, query, score, query_start, query_end)

    ``query_start`` and ``query_end`` are already ungapped query coordinates.
    ``pos_wo_gap`` therefore returns them unchanged because ``query`` itself is
    ungapped.  PairwiseAligner returns all optimal local alignments lazily; all
    of them are retained here because LoRTIA checks their positions when
    deciding whether an adapter is correctly placed.
    """
    # PairwiseAligner raises ValueError on an empty sequence.  Empty/missing
    # terminal sequence simply means that no adapter can be demonstrated at
    # this alignment end, so preserve LoRTIA semantics by returning no hit.
    if sequence is None or adapter is None:
        return []
    sequence = str(sequence)
    adapter = str(adapter)
    if not sequence or not adapter:
        return []

    aligner = _get_adapter_aligner(args)
    pairwise_alignments = aligner.align(adapter, sequence)

    # Preserve the previous behavior: no positive-scoring local alignment
    # means no adapter alignment. PairwiseAligner likewise has no alignments at
    # score 0, but keep the explicit guard for compatibility and clarity.
    if pairwise_alignments.score <= 0:
        return []

    alignments = []
    for alignment in pairwise_alignments:
        query_start = int(alignment.coordinates[1, 0])
        query_end = int(alignment.coordinates[1, -1])
        alignments.append((adapter,
                           sequence,
                           float(alignment.score),
                           query_start,
                           query_end))
    return alignments


def pos_wo_gap(alignment, p):
    """
    Returns the position in the alignment, substracting the gaps.
    """
    # the alignment end is reported counting the gaps as well, we change that
    pos_wo_gap = alignment[p] - alignment[1][:alignment[p]].count("-")
    return pos_wo_gap

def in_place_checker(alignments,
                     cigar_info,
                     args,
                     check_in_soft,
                     matchseq,
                     adapter):
    """
    Checks whether the adapter is no further than nts away from the start of
    the alignment.
    """
    ts = check_in_soft + args.shs_for_ts
    soft_match = check_in_soft + args.check_in_match
    adapter_info = (alignments[0][2],
                    pos_wo_gap(alignments[0], 3),
                    pos_wo_gap(alignments[0], 4),
                    "out of place")
    for alignment in alignments:
        if pos_wo_gap(alignment, 4) in range(check_in_soft 
                                              - args.check_from_alignment,
                                              soft_match + 1):
            adapter_info = (alignment[2],
                            pos_wo_gap(alignment, 3),
                            pos_wo_gap(alignment, 4),
                            "correct")
            break
    if adapter_info[3] == "correct":
        exon = 0
        for event in cigar_info:
            event_type, event_count = event
            if event_type == 0:
                if exon + int(event_count) < args.first_exon + 1:
                    exon += int(event_count)
                else:
                    break
            elif event_type == 3:
                match_adapter = adapter_aligner(matchseq, adapter, args)
                if match_adapter:
                    score_limit = args.match_in_first * 0.5 * args.match_score
                    if match_adapter[0][2] >= score_limit:
                        adapter_info = (alignment[2],
                                        pos_wo_gap(alignment, 3),
                                        pos_wo_gap(alignment, 4),
                                        "false exon")
                break
            else:
                for alignment in alignments:
                    if pos_wo_gap(alignment, 4) in range(ts, soft_match + 1):
                        adapter_info = (alignment[2],
                                        pos_wo_gap(alignment, 3),
                                        pos_wo_gap(alignment, 4),
                                        "potential template switching")
                        break
    return adapter_info

def get_adapter_info(sequence,
                     adapter,
                     score_limit,
                     cigar_info,
                     args,
                     check_in_soft,
                     matchseq):
    """
    Sends sequence to the aligner, and gathers adapter information.
    """
    alignments = adapter_aligner(sequence, adapter, args)

    if alignments: # local alignment may be empty when there is no positive score.
        # Pairwise2 only lists alignments with the best score, so checking the
        # 1st one is enough.

        isadapter = float(alignments[0][2]) > score_limit
        if isadapter:
            adapter_info = in_place_checker(alignments,
                                            cigar_info,
                                            args,
                                            check_in_soft,
                                            matchseq,
                                            adapter)
        else:
            adapter_info = (alignments[0][2],
                            pos_wo_gap(alignments[0], 3),
                            pos_wo_gap(alignments[0], 4),
                            "missing")
    else:
        adapter_info = (0, 0, 0, "missing")
    return adapter_info

def adapter_checker(read, args, full_forward_seq=None):
    """
    Retrieves the end sequences and sorts adapter information for each end.

    dRNA mode (args.dRNA):
      - 5' end is ALWAYS treated as 'correct' (independent of softclip)
      - 5' adapter search is skipped
      - 3' end QC is performed as usual using three_adapter on the TRUE 3' end
        (+: r3,  -: l3)
    """
    cigar = read.cigartuples
    seq = _sequence_for_alignment_orientation(read, full_forward_seq)

    leftend, check_in_softl, left_match = get_left_end(seq, cigar, args)
    rightend, check_in_softr, right_match = get_right_end(seq, cigar, args)

    left_cigar_info = cigar
    right_cigar_info = list(cigar)[::-1]

    # --- dRNA: 5' always correct, 3' QC only on true 3' end ---
    if getattr(args, "dRNA", False):

        # defaults
        three_left  = (0, 0, 0, "missing")
        three_right = (0, 0, 0, "missing")
        five_left   = (0, 0, 0, "missing")
        five_right  = (0, 0, 0, "missing")

        if read.is_reverse:
            # (-) strand: 5' = RIGHT  -> r5 correct
            #            3' = LEFT   -> l3 QC
            five_right = (0, 0, 0, "correct")

            if leftend is not None:
                three_left = get_adapter_info(
                    Seq(leftend).complement(),
                    args.three_adapter,
                    args.three_score,
                    left_cigar_info,
                    args,
                    check_in_softl,
                    Seq(left_match).complement()
                )
            else:
                three_left = (0, 0, 0, "missing")

        else:
            # (+) strand: 5' = LEFT   -> l5 correct
            #            3' = RIGHT  -> r3 QC
            five_left = (0, 0, 0, "correct")

            if rightend is not None:
                three_right = get_adapter_info(
                    rightend[::-1],
                    args.three_adapter,
                    args.three_score,
                    right_cigar_info,
                    args,
                    check_in_softr,
                    right_match[::-1]
                )
            else:
                three_right = (0, 0, 0, "missing")

        return {"l3": three_left,
                "r3": three_right,
                "l5": five_left,
                "r5": five_right}

    # --- adapterless PacBio full-length cDNA mode ---
    # Search only for the 3' poly(A) signal on both alignment ends.  When one
    # end has a unique *correct* 3' call, the opposite physical end is the 5'
    # end and is marked correct without performing a 5' adapter search.  The
    # existing set_read_strand() logic then flips/retains transcript strand
    # from this end evidence while preserving mapper identity bits.
    if getattr(args, "pacbio", False):
        three_left  = (0, 0, 0, "missing")
        three_right = (0, 0, 0, "missing")
        five_left   = (0, 0, 0, "missing")
        five_right  = (0, 0, 0, "missing")

        if leftend is not None:
            three_left = get_adapter_info(
                Seq(leftend).complement(),
                args.three_adapter,
                args.three_score,
                left_cigar_info,
                args,
                check_in_softl,
                Seq(left_match).complement()
            )
        if rightend is not None:
            three_right = get_adapter_info(
                rightend[::-1],
                args.three_adapter,
                args.three_score,
                right_cigar_info,
                args,
                check_in_softr,
                right_match[::-1]
            )

        left_poly_a = (three_left[3] == "correct")
        right_poly_a = (three_right[3] == "correct")

        # Require an unambiguous terminal poly(A) call for full-length
        # orientation.  A right 3' end implies a left 5' end; a left 3' end
        # implies a right 5' end.
        if right_poly_a and not left_poly_a:
            five_left = (0, 0, 0, "correct")
        elif left_poly_a and not right_poly_a:
            five_right = (0, 0, 0, "correct")

        return {"l3": three_left,
                "r3": three_right,
                "l5": five_left,
                "r5": five_right}

    # --- default mode (cDNA): original logic ---
    if leftend is not None:
        three_left = get_adapter_info(Seq(leftend).complement(),
                                      args.three_adapter,
                                      args.three_score,
                                      left_cigar_info,
                                      args,
                                      check_in_softl,
                                      Seq(left_match).complement())
        five_left = get_adapter_info(leftend,
                                     args.five_adapter,
                                     args.five_score,
                                     left_cigar_info,
                                     args,
                                     check_in_softl,
                                     left_match)
    else:
        three_left = (0, 0, 0, "missing")
        five_left = (0, 0, 0, "missing")

    if rightend is not None:
        three_right = get_adapter_info(rightend[::-1],
                                       args.three_adapter,
                                       args.three_score,
                                       right_cigar_info,
                                       args,
                                       check_in_softr,
                                       right_match[::-1])
        five_right = get_adapter_info(Seq(rightend[::-1]).complement(),
                                      args.five_adapter,
                                      args.five_score,
                                      right_cigar_info,
                                      args,
                                      check_in_softr,
                                      Seq(right_match[::-1]).complement())
    else:
        three_right = (0, 0, 0, "missing")
        five_right = (0, 0, 0, "missing")

    return {"l3": three_left,
            "r3": three_right,
            "l5": five_left,
            "r5": five_right}

def out_appender(line, out_file):
    """
    Appends string to a file.
    """
    with open(out_file, "a") as out:
        out.write(line)

def out_writer(dictionary, prefix, suffix):
    """
    Writes a dictionary into a csv file.
    """
    dataframe = pd.DataFrame.from_dict(dictionary, orient="index")
    dataframe.to_csv("{}_{}.tsv".format(prefix, suffix), header=False,
                     sep="\t")

def set_split_transcript_strand(old_flag, direction):
    """Set transcript genomic strand for an oriented split segment.

    ``direction`` is the biological 5'->3' direction along the original query.
    For 5to3 the transcript genomic strand equals the mapper alignment strand;
    for 3to5 it is inverted.  All mapper identity bits (supplementary,
    secondary, pairing, etc.) are preserved and the legacy QC-fail bit is
    cleared because orientation has been resolved at molecule level.
    """
    old_flag = int(old_flag)
    original_reverse = bool(old_flag & 0x10)
    transcript_reverse = original_reverse if direction == "5to3" else not original_reverse
    new_flag = old_flag & ~0x10 & ~0x200
    if transcript_reverse:
        new_flag |= 0x10
    return new_flag


def set_pacbio_strand(summary, old_flag):
    """Orient an adapterless PacBio read from a unique correct 3' poly(A) end.

    Only a *correct* terminal poly(A) call is authoritative.  This avoids the
    legacy cDNA strand setter treating weaker states such as template-switching
    or out-of-place adapter hits as sufficient PacBio orientation evidence.
    Mapper identity bits are preserved; unresolved/bidirectional poly(A)
    evidence is marked QC-fail (0x200) and excluded from stranded_only.bam.
    """
    left3 = summary.get("l3", (0, 0, 0, "missing"))[3] == "correct"
    right3 = summary.get("r3", (0, 0, 0, "missing"))[3] == "correct"

    new_flag = int(old_flag) & ~0x10 & ~0x200
    if left3 and not right3:
        new_flag |= 0x10
    elif right3 and not left3:
        pass
    else:
        new_flag |= 0x200
    return new_flag


def set_read_strand(summary, old_flag):
    """Set only the reverse-strand bit while preserving alignment identity.

    The legacy implementation replaced the complete FLAG with 0/16, which
    erased supplementary (0x800), secondary (0x100), pairing and other flags.
    Ambiguous adapter orientation is still marked with 0x200 so the existing
    stranded-only filter can exclude it.
    """
    absent = {"missing", "no softclip", "hardclip unavailable", "sequence unavailable"}
    is_l3 = summary.get("l3")[3] not in absent
    is_r3 = summary.get("r3")[3] not in absent
    is_l5 = summary.get("l5")[3] not in absent
    is_r5 = summary.get("r5")[3] not in absent

    new_flag = int(old_flag) & ~0x10
    if (is_l3 or is_r5) and not (is_r3 or is_l5):
        new_flag |= 0x10
    elif (is_r3 or is_l5) and not (is_l3 or is_r5):
        pass
    else:
        new_flag |= 0x200
    return new_flag

def intron_finder(read, read_start, read_end, args, split_meta=None):
    """
    Finds introns based on the cigar code. Also filters introns based on the
    first_exon argument and removes introns which contain a long insertion
    preceding the intron (likely triple-chimeric reads).

    The reference position is tracked in a single pass through the CIGAR.
    This preserves the original LoRTIA coordinate logic while avoiding a
    repeated scan of all preceding CIGAR operations for every intron.

    dRNA fix:
      - only the biologically relevant end is used for is_left_false/is_right_false
        (+ strand: LEFT=5' (l5), RIGHT=3' (r3)
         - strand: LEFT=3' (l3), RIGHT=5' (r5))
      - prevents irrelevant "no softclip" tags from killing introns.
    """
    introns = ()
    previous_event = (0, 0)
    matches_so_far = 0
    matches_to_come = 0

    # Original intron coordinates are calculated as:
    #   read_start - 1 + sum(previous M/D/N lengths)
    # Keep that exact convention, but update it incrementally.
    reference_position = read_start - 1

    # --- false-end logic ---
    # For a split molecule, only the two *external biological molecule ends*
    # are allowed to invalidate a short terminal exon.  A fusion breakpoint is
    # an internal molecule boundary and must not be treated as a missing
    # adapter/TSS/TES artefact.
    split_oriented = bool(
        split_meta and split_meta.get("direction") in {"5to3", "3to5"}
    )
    if split_oriented:
        idx = int(split_meta["segment_index"])
        count = int(split_meta["segment_count"])
        if split_meta["direction"] == "5to3":
            is_bio_first = (idx == 1)
            is_bio_last = (idx == count)
        else:
            is_bio_first = (idx == count)
            is_bio_last = (idx == 1)

        if read.is_reverse:
            # On the transcript minus strand, genomic LEFT is biological 3'.
            left_external = is_bio_last
            left_tag = "l3"
        else:
            # On the transcript plus strand, genomic LEFT is biological 5'.
            left_external = is_bio_first
            left_tag = "l5"
        is_left_false = bool(
            left_external and read.get_tag(left_tag).split(",")[3] != "correct"
        )
    elif getattr(args, "dRNA", False):
        if read.is_reverse:
            is_left_false = (read.get_tag("l3").split(",")[3] != "correct")
        else:
            is_left_false = (read.get_tag("l5").split(",")[3] != "correct")
    elif getattr(args, "pacbio", False):
        # PacBio orientation has already been resolved from terminal poly(A)
        # evidence by set_read_strand(); only the biological terminal tag is
        # relevant for short-first-exon QC.
        if read.is_reverse:
            is_left_false = (read.get_tag("l3").split(",")[3] != "correct")
        else:
            is_left_false = (read.get_tag("l5").split(",")[3] != "correct")
    else:
        # Preserve legacy behavior for non-split cDNA reads.
        is_left_false = ((read.get_tag("l5").split(",")[3] != "correct") or
                         (read.get_tag("l3").split(",")[3] != "correct"))

    for event_type, event_count in read.cigartuples:
        # event coding: M:0, I:1, D:2, N:3, S:4, H:5
        previous_is_insert = (previous_event[0] == 1 and
                              previous_event[1] > args.insert_before_intron)

        if event_type == 3:  # N = intron
            intron_start = reference_position
            intron_end = intron_start + event_count + 1
            introns += (intron_start, intron_end),

            matches_to_come = 0
            if (previous_is_insert or
                (is_left_false and (matches_so_far < args.first_exon))):
                read.set_tag("ga", ",".join((str(i) for i in introns)), "Z")
                introns = ()

        elif event_type == 0:  # M
            matches_so_far += event_count
            matches_to_come += event_count

        # Preserve the original coordinate convention exactly: only M, D and N
        # advance the reference position here (the legacy code used [0, 2, 3]).
        if event_type in [0, 2, 3]:
            reference_position += event_count

        previous_event = (event_type, event_count)

    # --- right end logic ---
    if split_oriented:
        if read.is_reverse:
            # On the transcript minus strand, genomic RIGHT is biological 5'.
            right_external = is_bio_first
            right_tag = "r5"
        else:
            # On the transcript plus strand, genomic RIGHT is biological 3'.
            right_external = is_bio_last
            right_tag = "r3"
        is_right_false = bool(
            right_external and read.get_tag(right_tag).split(",")[3] != "correct"
        )
    elif getattr(args, "dRNA", False):
        if read.is_reverse:
            is_right_false = (read.get_tag("r5").split(",")[3] != "correct")
        else:
            is_right_false = (read.get_tag("r3").split(",")[3] != "correct")
    elif getattr(args, "pacbio", False):
        if read.is_reverse:
            is_right_false = (read.get_tag("r5").split(",")[3] != "correct")
        else:
            is_right_false = (read.get_tag("r3").split(",")[3] != "correct")
    else:
        is_right_false = ((read.get_tag("r5").split(",")[3] != "correct") or
                          (read.get_tag("r3").split(",")[3] != "correct"))

    if matches_to_come < args.first_exon and is_right_false and introns:
        read.set_tag("ga", str(introns[-1]), "Z")
        introns = introns[:-1]

    return introns

def put_in_dict(adapter_sum, end_type, string, dictionary, contig, end):
    """
    Increases adapter value in a dicrionary.
    """
    if adapter_sum.get(end_type)[3] == string:
        if (contig, end) in dictionary:
            dictionary[contig, end] += 1
        else:
            dictionary[contig, end] = 1
    return dictionary

def prepare_new_sam_line(read,
                         args,
                         out_handle,
                         l3_dict,
                         r3_dict,
                         l5_dict,
                         r5_dict,
                         introns_dict,
                         tsl3_dict,
                         tsr3_dict,
                         tsl5_dict,
                         tsr5_dict,
                         full_forward_seq=None,
                         adapter_sum=None,
                         original_flag=None,
                         split_meta=None):
    """Summarize one alignment and write it to the processed SAM.

    For split reads, adapter_sum is precomputed while every alignment still has
    its original mapper FLAG. original_flag is preserved in ZF so downstream
    fusion reconstruction can recover the alignment orientation even though the
    legacy LoRTIA cDNA path rewrites bit 0x10 to encode transcript strand.
    """
    read_start = read.reference_start + 1
    read_end = read.reference_end
    contig = read.reference_name

    if original_flag is None:
        original_flag = int(read.flag)
    if adapter_sum is None:
        adapter_sum = adapter_checker(read, args, full_forward_seq=full_forward_seq)

    # Preserve the mapper's original alignment FLAG before LoRTIA changes it.
    read.set_tag("ZF", int(original_flag), "i")

    if split_meta is not None:
        read.set_tag("ZN", int(split_meta["segment_count"]), "i")
        read.set_tag("ZI", int(split_meta["segment_index"]), "i")
        read.set_tag("ZS", int(split_meta["qstart"]), "i")
        read.set_tag("ZE", int(split_meta["qend"]), "i")
        read.set_tag("ZL", int(split_meta["qlen"]), "i")
        read.set_tag("ZO", str(split_meta["direction"]), "Z")
        read.set_tag("ZQ", str(split_meta.get("orientation_source", "unresolved")), "Z")
        read.set_tag("ZT", str(split_meta.get("ts_consensus", "missing")), "Z")
        read.set_tag("ZP", str(split_meta.get("adapter_direction", "ambiguous")), "Z")
        read.set_tag("ZA", int(split_meta["internal_adapter"]), "i")
        read.set_tag("ZU", int(split_meta["full_sequence_recovered"]), "i")
        # Split-read QC only; these tags never determine whether a fusion is called.
        read.set_tag("ZC", float(split_meta.get("query_coverage", 0.0)), "f")
        read.set_tag("ZX", int(split_meta.get("alternative_mapping_count", 0)), "i")
        read.set_tag("ZY", int(split_meta.get("secondary_alignment_count", 0)), "i")
        read.set_tag("ZV", int(split_meta.get("mapping_ambiguous", 0)), "i")

    if split_meta is not None and split_meta.get("direction") in {"5to3", "3to5"}:
        # For a resolved split molecule, molecule-level orientation (preferably
        # minimap2 ts consensus) is authoritative for transcript strand.  Do
        # not mark an internal segment as QC-failed merely because that segment
        # lacks a terminal adapter; adapter evidence remains stored separately.
        read.flag = set_split_transcript_strand(
            int(original_flag), split_meta["direction"]
        )
        read.mapping_quality = 9
    elif getattr(args, "dRNA", False):
        read.mapping_quality = 9
    elif getattr(args, "pacbio", False):
        read.flag = set_pacbio_strand(adapter_sum, int(original_flag))
        if read.flag & 0x200:
            read.mapping_quality = 0
        else:
            read.mapping_quality = 9
    else:
        read.flag = set_read_strand(adapter_sum, int(original_flag))
        if read.flag & 0x200:
            read.mapping_quality = 0
        else:
            read.mapping_quality = 9

    l3_dict = put_in_dict(adapter_sum, "l3", "correct", l3_dict, contig, read_start)
    r3_dict = put_in_dict(adapter_sum, "r3", "correct", r3_dict, contig, read_end)
    l5_dict = put_in_dict(adapter_sum, "l5", "correct", l5_dict, contig, read_start)
    r5_dict = put_in_dict(adapter_sum, "r5", "correct", r5_dict, contig, read_end)
    tsl3_dict = put_in_dict(adapter_sum, "l3", "potential template switching", tsl3_dict, contig, read_start)
    tsr3_dict = put_in_dict(adapter_sum, "r3", "potential template switching", tsr3_dict, contig, read_end)
    tsl5_dict = put_in_dict(adapter_sum, "l5", "potential template switching", tsl5_dict, contig, read_start)
    tsr5_dict = put_in_dict(adapter_sum, "r5", "potential template switching", tsr5_dict, contig, read_end)

    read.set_tag("re", read_end, "i")
    read.set_tag("l3", ",".join((str(i) for i in adapter_sum.get("l3"))), "Z")
    read.set_tag("r3", ",".join((str(i) for i in adapter_sum.get("r3"))), "Z")
    read.set_tag("l5", ",".join((str(i) for i in adapter_sum.get("l5"))), "Z")
    read.set_tag("r5", ",".join((str(i) for i in adapter_sum.get("r5"))), "Z")

    introns = intron_finder(read, read_start, read_end, args, split_meta=split_meta)
    if introns:
        for intron in introns:
            if (contig, intron) in introns_dict:
                introns_dict[contig, intron] += 1
            else:
                introns_dict[contig, intron] = 1
    read.set_tag("in", ",".join((str(i) for i in introns)), "Z")
    out_handle.write(read.to_string() + "\n")
    return (l3_dict, r3_dict, l5_dict, r5_dict, introns_dict,
            tsl3_dict, tsr3_dict, tsl5_dict, tsr5_dict)


###############################################################################
###           Legacy-compatible conventional transcript processing           ###
###############################################################################
# Compatibility baseline: Samprocessor(20260906-121723).py.
# Preserve its terminal S/H slicing exactly for conventional cDNA/dRNA output.
# The separate split/fusion path retains its full-sequence-aware clip handling.
# These functions intentionally preserve the pre-fusion LoRTIA linear path.
# Fusion/split-molecule reconstruction is performed in a separate second pass
# and must never alter conventional l5/r5/l3/r3/in counts.

def legacy_get_left_end(read_seq, cigar, args):
    """
    Gets the left end of the alignment, taking 'check_in_soft' number of 
    nucleotides from the softclip and check_in_match number of nucleotides
    from the mapped part.
    """
    if cigar[0][0] in [4, 5]:
        softleftpos = cigar[0][1]
        if softleftpos > args.check_in_soft - 1:
            checkleft = softleftpos - args.check_in_soft
            cis = args.check_in_soft
        else:
            checkleft = 0
            cis = softleftpos
        leftend_seq = read_seq[checkleft:softleftpos + args.check_in_match]
        left_match = read_seq[softleftpos:softleftpos + args.match_in_first]
    else:
        leftend_seq = "no softclip"
        cis = 0
        left_match = 0
    return leftend_seq, cis, left_match


def legacy_get_right_end(read_seq, cigar, args):
    """
    Gets the right end of the alignment, taking 'check_in_soft' number of 
    nucleotides from the softclip and check_in_match number of nucleotides
    from the mapped part.
    """
# here I first calculate the length of the read and I
# substract the length of the closing softclip
    if cigar[-1][0] in [4, 5]:
        softrightpos = (len(read_seq) - cigar[-1][1])
        if args.check_in_soft < cigar[-1][1] - 1:
            checkright = softrightpos + args.check_in_soft
            cis = args.check_in_soft
        else:
            checkright = len(read_seq)
            cis = cigar[-1][1]
        rightend_seq = read_seq[softrightpos - args.check_in_match:checkright]
        right_match = read_seq[softrightpos - args.match_in_first:softrightpos]
    else:
        rightend_seq = "no softclip"
        cis = 0
        right_match = 0
    return rightend_seq, cis, right_match


def legacy_adapter_checker(read, args):
    """Pre-fusion adapter calling, using only each alignment record's own SEQ."""
    cigar = read.cigartuples
    seq = read.query_sequence

    leftend, check_in_softl, left_match = legacy_get_left_end(seq, cigar, args)
    rightend, check_in_softr, right_match = legacy_get_right_end(seq, cigar, args)

    left_cigar_info = cigar
    right_cigar_info = list(cigar)[::-1]

    if getattr(args, "dRNA", False):
        three_left  = (0, 0, 0, "no softclip")
        three_right = (0, 0, 0, "no softclip")
        five_left   = (0, 0, 0, "no softclip")
        five_right  = (0, 0, 0, "no softclip")

        if read.is_reverse:
            five_right = (0, 0, 0, "correct")
            if leftend != "no softclip":
                three_left = get_adapter_info(
                    Seq(leftend).complement(), args.three_adapter,
                    args.three_score, left_cigar_info, args,
                    check_in_softl, Seq(left_match).complement()
                )
        else:
            five_left = (0, 0, 0, "correct")
            if rightend != "no softclip":
                three_right = get_adapter_info(
                    rightend[::-1], args.three_adapter, args.three_score,
                    right_cigar_info, args, check_in_softr,
                    right_match[::-1]
                )

        return {"l3": three_left, "r3": three_right,
                "l5": five_left, "r5": five_right}

    if leftend != "no softclip":
        three_left = get_adapter_info(
            Seq(leftend).complement(), args.three_adapter, args.three_score,
            left_cigar_info, args, check_in_softl,
            Seq(left_match).complement()
        )
        five_left = get_adapter_info(
            leftend, args.five_adapter, args.five_score,
            left_cigar_info, args, check_in_softl, left_match
        )
    else:
        three_left = (0, 0, 0, "no softclip")
        five_left = (0, 0, 0, "no softclip")

    if rightend != "no softclip":
        three_right = get_adapter_info(
            rightend[::-1], args.three_adapter, args.three_score,
            right_cigar_info, args, check_in_softr, right_match[::-1]
        )
        five_right = get_adapter_info(
            Seq(rightend[::-1]).complement(), args.five_adapter,
            args.five_score, right_cigar_info, args, check_in_softr,
            Seq(right_match[::-1]).complement()
        )
    else:
        three_right = (0, 0, 0, "no softclip")
        five_right = (0, 0, 0, "no softclip")

    return {"l3": three_left, "r3": three_right,
            "l5": five_left, "r5": five_right}


def legacy_set_read_strand(summary, old_flag):
    """Exact pre-fusion LoRTIA strand/QC flag logic."""
    is_l3 = summary.get("l3")[3] not in ["missing", "no softclip"]
    is_r3 = summary.get("r3")[3] not in ["missing", "no softclip"]
    is_l5 = summary.get("l5")[3] not in ["missing", "no softclip"]
    is_r5 = summary.get("r5")[3] not in ["missing", "no softclip"]
    if (is_l3 or is_r5) and not (is_r3 or is_l5):
        new_flag = 16
    elif (is_r3 or is_l5) and not (is_l3 or is_r5):
        new_flag = 0
    else:
        new_flag = old_flag + 512
    return new_flag


def legacy_intron_finder(read, read_start, read_end, args):
    """Exact pre-fusion intron acceptance logic for cDNA/dRNA."""
    introns = ()
    previous_event = (0, 0)
    matches_so_far = 0
    matches_to_come = 0
    reference_position = read_start - 1

    if getattr(args, "dRNA", False):
        if read.is_reverse:
            is_left_false = (read.get_tag("l3").split(",")[3] != "correct")
        else:
            is_left_false = (read.get_tag("l5").split(",")[3] != "correct")
    else:
        is_left_false = ((read.get_tag("l5").split(",")[3] != "correct") or
                         (read.get_tag("l3").split(",")[3] != "correct"))

    for event_type, event_count in read.cigartuples:
        previous_is_insert = (
            previous_event[0] == 1 and
            previous_event[1] > args.insert_before_intron
        )
        if event_type == 3:
            intron_start = reference_position
            intron_end = intron_start + event_count + 1
            introns += (intron_start, intron_end),
            matches_to_come = 0
            if (previous_is_insert or
                    (is_left_false and (matches_so_far < args.first_exon))):
                read.set_tag("ga", ",".join((str(i) for i in introns)), "Z")
                introns = ()
        elif event_type == 0:
            matches_so_far += event_count
            matches_to_come += event_count

        if event_type in [0, 2, 3]:
            reference_position += event_count
        previous_event = (event_type, event_count)

    if getattr(args, "dRNA", False):
        if read.is_reverse:
            is_right_false = (read.get_tag("r5").split(",")[3] != "correct")
        else:
            is_right_false = (read.get_tag("r3").split(",")[3] != "correct")
    else:
        is_right_false = ((read.get_tag("r5").split(",")[3] != "correct") or
                          (read.get_tag("r3").split(",")[3] != "correct"))

    if matches_to_come < args.first_exon and is_right_false and introns:
        read.set_tag("ga", str(introns[-1]), "Z")
        introns = introns[:-1]
    return introns


def linear_prepare_new_sam_line(read, args, out_handle,
                                l3_dict, r3_dict, l5_dict, r5_dict,
                                introns_dict, tsl3_dict, tsr3_dict,
                                tsl5_dict, tsr5_dict):
    """Conventional LoRTIA processing isolated from split/fusion reconstruction."""
    read_start = read.reference_start + 1
    read_end = read.reference_end
    contig = read.reference_name

    if getattr(args, "pacbio", False):
        # PacBio is a new linear library mode. It remains isolated from the
        # fusion pass but uses the PacBio-specific terminal poly(A) logic.
        adapter_sum = adapter_checker(read, args, full_forward_seq=None)
        read.flag = set_pacbio_strand(adapter_sum, read.flag)
        read.mapping_quality = 0 if (read.flag & 0x200) else 9
    else:
        # Standard cDNA and dRNA reproduce the tested pre-fusion baseline.
        adapter_sum = legacy_adapter_checker(read, args)
        if getattr(args, "dRNA", False):
            read.mapping_quality = 9
        else:
            read.flag = legacy_set_read_strand(adapter_sum, read.flag)
            read.mapping_quality = 0 if read.flag > 511 else 9

    l3_dict = put_in_dict(adapter_sum, "l3", "correct", l3_dict, contig, read_start)
    r3_dict = put_in_dict(adapter_sum, "r3", "correct", r3_dict, contig, read_end)
    l5_dict = put_in_dict(adapter_sum, "l5", "correct", l5_dict, contig, read_start)
    r5_dict = put_in_dict(adapter_sum, "r5", "correct", r5_dict, contig, read_end)
    tsl3_dict = put_in_dict(adapter_sum, "l3", "potential template switching", tsl3_dict, contig, read_start)
    tsr3_dict = put_in_dict(adapter_sum, "r3", "potential template switching", tsr3_dict, contig, read_end)
    tsl5_dict = put_in_dict(adapter_sum, "l5", "potential template switching", tsl5_dict, contig, read_start)
    tsr5_dict = put_in_dict(adapter_sum, "r5", "potential template switching", tsr5_dict, contig, read_end)

    read.set_tag("re", read_end, "i")
    read.set_tag("l3", ",".join((str(i) for i in adapter_sum.get("l3"))), "Z")
    read.set_tag("r3", ",".join((str(i) for i in adapter_sum.get("r3"))), "Z")
    read.set_tag("l5", ",".join((str(i) for i in adapter_sum.get("l5"))), "Z")
    read.set_tag("r5", ",".join((str(i) for i in adapter_sum.get("r5"))), "Z")

    if getattr(args, "pacbio", False):
        introns = intron_finder(read, read_start, read_end, args, split_meta=None)
    else:
        introns = legacy_intron_finder(read, read_start, read_end, args)

    if introns:
        for intron in introns:
            if (contig, intron) in introns_dict:
                introns_dict[contig, intron] += 1
            else:
                introns_dict[contig, intron] = 1
    read.set_tag("in", ",".join((str(i) for i in introns)), "Z")
    out_handle.write(read.to_string() + "\n")

    return (l3_dict, r3_dict, l5_dict, r5_dict, introns_dict,
            tsl3_dict, tsr3_dict, tsl5_dict, tsr5_dict)


def linear_deal_with_same_name(previous_lines, args, out_handle,
                               l3_dict, r3_dict, l5_dict, r5_dict,
                               introns_dict, tsl3_dict, tsr3_dict,
                               tsl5_dict, tsr5_dict):
    """Exact pre-fusion same-QNAME selection for the conventional caller."""
    counter = 1
    previous_lines.sort(key=operator.itemgetter(3), reverse=True)
    ranges_by_far = ()
    for read, read_start, read_end, read_span in previous_lines:
        is_in_previous = False
        for ran in ranges_by_far:
            if read_start in ran or read_end in ran:
                is_in_previous = True
                break
        if not is_in_previous:
            if len(previous_lines) > 1:
                read.query_name = "_".join((read.query_name, str(counter)))
                counter += 1
            (l3_dict, r3_dict, l5_dict, r5_dict, introns_dict,
             tsl3_dict, tsr3_dict, tsl5_dict, tsr5_dict) = linear_prepare_new_sam_line(
                read, args, out_handle,
                l3_dict, r3_dict, l5_dict, r5_dict, introns_dict,
                tsl3_dict, tsr3_dict, tsl5_dict, tsr5_dict
            )
        ranges_by_far += range(read_start, read_end + 1),

    return (l3_dict, r3_dict, l5_dict, r5_dict, introns_dict,
            tsl3_dict, tsr3_dict, tsl5_dict, tsr5_dict)


def _safe_read_tag(read, tag, default=None):
    try:
        return read.get_tag(tag)
    except Exception:
        return default


def _write_chimeric_group(chim_handle, reads, full_forward_seq):
    """Write retained split-alignment structure for transparent fusion QC."""
    if chim_handle is None or len(reads) < 2:
        return
    for read in reads:
        tags = {}
        for tag in ("l3", "r3", "l5", "r5"):
            try:
                tags[tag] = read.get_tag(tag).split(",")[3]
            except Exception:
                tags[tag] = "missing"
        chim_handle.write("\t".join(map(str, [
            read.query_name,
            read.get_tag("ZI"), read.get_tag("ZN"),
            read.get_tag("ZS"), read.get_tag("ZE"), read.get_tag("ZL"),
            read.reference_name, read.reference_start + 1, read.reference_end,
            "-" if (read.get_tag("ZF") & 0x10) else "+",
            int(read.is_supplementary), int(read.is_secondary),
            read.get_tag("ZF"), read.flag, read.mapping_quality,
            read.get_tag("ZU"), read.get_tag("ZO"),
            _safe_read_tag(read, "ZQ", "unresolved"),
            _safe_read_tag(read, "ZT", "missing"),
            _safe_read_tag(read, "ZP", "ambiguous"),
            _safe_read_tag(read, "ts", "."),
            read.get_tag("ZA"),
            _safe_read_tag(read, "ZC", 0.0),
            _safe_read_tag(read, "ZX", 0),
            _safe_read_tag(read, "ZY", 0),
            _safe_read_tag(read, "ZV", 0),
            tags["l5"], tags["r5"], tags["l3"], tags["r3"],
        ])) + "\n")


def deal_with_same_name(previous_lines,
                        args,
                        out_handle,
                        l3_dict,
                        r3_dict,
                        l5_dict,
                        r5_dict,
                        introns_dict,
                        tsl3_dict,
                        tsr3_dict,
                        tsl5_dict,
                        tsr5_dict,
                        chim_handle=None):
    """Process all alignments belonging to one physical read jointly.

    Distinct query segments are retained, alternative overlapping mappings are
    collapsed, and one complete same-QNAME sequence is reused when available.
    For split molecules, adapter calls are first computed for every segment
    while mapper flags are intact; the two physical molecule termini are then
    evaluated jointly to determine 5'->3' molecule direction. This permits the
    downstream Fusion_Annotator to validate TSS on one sub-alignment and TES on
    another without requiring both adapters on a single BAM record.
    """
    if not previous_lines:
        return (l3_dict, r3_dict, l5_dict, r5_dict, introns_dict,
                tsl3_dict, tsr3_dict, tsl5_dict, tsr5_dict)

    reads = [item[0] for item in previous_lines]
    full_forward_seq = _group_forward_sequence(reads)
    max_qov = int(getattr(args, "max_split_query_overlap", 15))
    selected, split_qc = _select_split_alignments(reads, max_query_overlap=max_qov)

    # Critical: compute every adapter summary before legacy strand rewriting can
    # change FLAG 0x10 on any record in the group.
    prepared = [
        (read, adapter_checker(read, args, full_forward_seq=full_forward_seq))
        for read in selected
    ]

    if len(prepared) > 1:
        orientation = _molecule_orientation(prepared, args)
        direction = orientation["direction"]
        internal_adapter = int(_internal_adapter_present(prepared, args))
        full_recovered = int(full_forward_seq is not None)
    else:
        orientation = {
            "direction": "ambiguous",
            "source": "single_segment",
            "ts_consensus": "missing",
            "adapter_direction": "ambiguous",
        }
        direction = "ambiguous"
        internal_adapter = 0
        full_recovered = int(full_forward_seq is not None)

    for idx, (read, adapter_sum) in enumerate(prepared, start=1):
        original_flag = int(read.flag)
        split_meta = None
        if len(prepared) > 1:
            qstart, qend, qlen = _query_interval_original_orientation(read)
            split_meta = {
                "segment_count": len(prepared),
                "segment_index": idx,
                "qstart": qstart,
                "qend": qend,
                "qlen": qlen,
                "direction": direction,
                "orientation_source": orientation["source"],
                "ts_consensus": orientation["ts_consensus"],
                "adapter_direction": orientation["adapter_direction"],
                "internal_adapter": internal_adapter,
                "full_sequence_recovered": full_recovered,
                "query_coverage": split_qc["query_coverage"],
                "alternative_mapping_count": split_qc["alternative_mapping_count"],
                "secondary_alignment_count": split_qc["secondary_alignment_count"],
                "mapping_ambiguous": split_qc["mapping_ambiguous"],
            }

        (l3_dict,
         r3_dict,
         l5_dict,
         r5_dict,
         introns_dict,
         tsl3_dict,
         tsr3_dict,
         tsl5_dict,
         tsr5_dict) = prepare_new_sam_line(
             read, args, out_handle,
             l3_dict, r3_dict, l5_dict, r5_dict, introns_dict,
             tsl3_dict, tsr3_dict, tsl5_dict, tsr5_dict,
             full_forward_seq=full_forward_seq,
             adapter_sum=adapter_sum,
             original_flag=original_flag,
             split_meta=split_meta,
         )

    _write_chimeric_group(chim_handle, selected, full_forward_seq)
    return (l3_dict, r3_dict, l5_dict, r5_dict, introns_dict,
            tsl3_dict, tsr3_dict, tsl5_dict, tsr5_dict)



def sam_iterator(args):
    """Run the conventional LoRTIA pass using the tested pre-fusion logic."""
    sam = pysam.AlignmentFile("{}_temp.sam".format(args.prefix), "r")
    previous_lines = []
    prog_comment = "@PG\tID:LoRTIA\tPN:LoRTIA\tVN:0.1\tCL:Samprocessor.py {}\n".format(args)
    outsam = pysam.AlignmentFile(args.out_file, "w", template=sam)
    outsam.close()
    out_appender(prog_comment, args.out_file)
    out_handle = open(args.out_file, "a")

    l3_dict = {}
    r3_dict = {}
    l5_dict = {}
    r5_dict = {}
    introns_dict = {}
    tsl3_dict = {}
    tsr3_dict = {}
    tsl5_dict = {}
    tsr5_dict = {}

    print("Iterating over {} (legacy-compatible linear pass)...".format(args.in_file))
    for read in sam:
        if not read.is_unmapped:
            read_start = read.reference_start + 1
            read_end = read.reference_end
            read_span = read_end - read_start
            if (previous_lines and
                    read.query_name != previous_lines[0][0].query_name):
                (l3_dict, r3_dict, l5_dict, r5_dict, introns_dict,
                 tsl3_dict, tsr3_dict, tsl5_dict, tsr5_dict) = linear_deal_with_same_name(
                    previous_lines, args, out_handle,
                    l3_dict, r3_dict, l5_dict, r5_dict, introns_dict,
                    tsl3_dict, tsr3_dict, tsl5_dict, tsr5_dict
                )
                previous_lines = []
            previous_lines.append((read, read_start, read_end, read_span))

    (l3_dict, r3_dict, l5_dict, r5_dict, introns_dict,
     tsl3_dict, tsr3_dict, tsl5_dict, tsr5_dict) = linear_deal_with_same_name(
        previous_lines, args, out_handle,
        l3_dict, r3_dict, l5_dict, r5_dict, introns_dict,
        tsl3_dict, tsr3_dict, tsl5_dict, tsr5_dict
    )

    out_handle.close()
    sam.close()

    print("Generating linear feature statistics...")
    out_writer(r3_dict, args.prefix, "r3")
    out_writer(l3_dict, args.prefix, "l3")
    out_writer(l5_dict, args.prefix, "l5")
    out_writer(r5_dict, args.prefix, "r5")
    out_writer(tsr3_dict, args.prefix, "ts_r3")
    out_writer(tsl3_dict, args.prefix, "ts_l3")
    out_writer(tsl5_dict, args.prefix, "ts_l5")
    out_writer(tsr5_dict, args.prefix, "ts_r5")
    out_writer(introns_dict, args.prefix, "in")



def fusion_sam_iterator(args):
    """Build a fusion-only processed SAM from the untouched name-sorted input.

    This second pass may use full same-QNAME sequence recovery, query-space
    split selection and molecule-level orientation, but none of its counts are
    written to the conventional l5/r5/l3/r3/in TSV files.
    """
    sam = pysam.AlignmentFile("{}_temp.sam".format(args.prefix), "r")
    fusion_sam = args.prefix + "_fusion_out.sam"
    outsam = pysam.AlignmentFile(fusion_sam, "w", template=sam)
    outsam.close()
    out_handle = open(fusion_sam, "a")

    previous_lines = []
    dummy = [{} for _ in range(9)]
    chim_path = args.prefix + "_chimeric_segments.tsv"
    chim_handle = open(chim_path, "w")
    chim_handle.write(
        "read_name\tsegment_index\tsegment_count\tquery_start\tquery_end\tquery_length\t"
        "contig\tref_start\tref_end\toriginal_alignment_strand\tis_supplementary\t"
        "is_secondary\toriginal_flag\tprocessed_flag\tmapq\tfull_sequence_recovered\t"
        "molecule_direction\torientation_source\tts_consensus\tadapter_direction\tts_tag\t"
        "internal_adapter\tquery_coverage\talternative_mapping_count\t"
        "secondary_alignment_count\tmapping_ambiguous\tl5\tr5\tl3\tr3\n"
    )

    print("Building isolated fusion/split-alignment stream...")
    for read in sam:
        if read.is_unmapped:
            continue
        read_start = read.reference_start + 1
        read_end = read.reference_end
        read_span = read_end - read_start
        if previous_lines and read.query_name != previous_lines[0][0].query_name:
            dummy = list(deal_with_same_name(
                previous_lines, args, out_handle,
                *dummy, chim_handle=chim_handle
            ))
            previous_lines = []
        previous_lines.append((read, read_start, read_end, read_span))

    if previous_lines:
        dummy = list(deal_with_same_name(
            previous_lines, args, out_handle,
            *dummy, chim_handle=chim_handle
        ))

    out_handle.close()
    chim_handle.close()
    sam.close()
    return fusion_sam


def fusion_output_creator(fusion_sam, threads=1):
    """Sort/index the fusion-only SAM without producing coverage files."""
    thr = max(int(threads), 1)
    fusion_bam = fusion_sam.replace(".sam", "_sorted.bam")
    pysam.sort("-@", str(thr), "-o", fusion_bam, fusion_sam)
    pysam.index("-@", str(thr), fusion_bam)
    return fusion_bam



def Samprocessor(args):
    """Run isolated conventional and (optionally) fusion processing passes."""
    if getattr(args, "dRNA", False) and getattr(args, "pacbio", False):
        raise ValueError("--dRNA and --pacbio are mutually exclusive library modes.")
    if int(getattr(args, "threads", 8)) < 1:
        raise ValueError("--threads must be >= 1.")
    if not os.path.isdir(args.out_path):
        os.mkdir(args.out_path)

    in_prefix = os.path.basename(args.in_file)
    args.prefix = "{}/{}".format(args.out_path, in_prefix[:len(in_prefix) - 4])
    args.out_file = "{}_out.sam".format(args.prefix)
    thr = max(int(getattr(args, "threads", 8)), 1)

    # Keep name sorting single-threaded to preserve the tested pre-fusion
    # same-QNAME ordering. Parallelism is retained for downstream BAM I/O.
    input_sorter(args.in_file, args.prefix)

    # Pass 1: legacy-compatible conventional feature evidence. Split/fusion
    # selection and molecule orientation are confined to the second pass.
    # Only this pass writes standard feature TSVs and coverage files.
    sam_iterator(args)
    output_creator(args.out_file, threads=thr)

    # Pass 2: fusion-specific split reconstruction. It reads the untouched
    # name-sorted SAM again and writes a separate tagged BAM.
    args.fusion_bam = None
    if not getattr(args, "no_fusion_detection", False):
        fusion_sam = fusion_sam_iterator(args)
        args.fusion_bam = fusion_output_creator(fusion_sam, threads=thr)

    print("Processed files saved to {}\n".format(args.out_path))


###############################################################################
###                             Main function                               ###
###############################################################################

def main():
    args = parsing()
    Samprocessor(args)

def parsing():
    """
    This part handles the commandline arguments
    """
    parser = ArgumentParser(description="This is the first module of LoRTIA:\
                            a Long-read RNA-Seq Transcript Isofom Annotator")
    parser.add_argument("in_file",
                        help="Input file. Both .sam and .bam files are\
                        accepted.",
                        metavar="input_file")
    parser.add_argument("out_path",
                        help="Output folder. Multiple output files are going\
                        to be created using the input file's prefix (ie. the\
                        part that precedes '.bam' or '.sam')",
                        metavar="output_path")
    parser.add_argument("--match_score", 
                        dest="match_score",
                        help="The alignment scores for each match when \
                        searching for adapters. Penalty scores should be \
                        supplied as negative vaules. The default is: 2",
                        type=float,
                        metavar="[float]", 
                        default=2.0)
    parser.add_argument("--mismatch_score", 
                        dest="mismatch_score",
                        help="The alignment scores for each mismatch when \
                        searching for adapters. Penalty scores should be \
                        supplied as negative vaules. The default is: -3",
                        type=float,
                        metavar="[float]", 
                        default=-3.0)
    parser.add_argument("--gap_open_score", 
                        dest="gap_open_score",
                        help="The alignment scores for each gap opening when\
                        searching for adapters. Penalty scores should be \
                        supplied as negative vaules. The default is: -3",
                        type=float,
                        metavar="[float]", 
                        default=-3.0)
    parser.add_argument("--gap_extend_score", 
                        dest="gap_extend_score",
                        help="The alignment scores for each gap extension \
                        when searching for adapters. Penalty scores should be \
                        supplied as negative vaules. The default is: -3",
                        type=float,
                        metavar="[float]", 
                        default=-3.0)
    parser.add_argument("-3", "--three_adapter", 
                        dest="three_adapter",
                        help="The 3' adapter to look for, the default is a \
                        polyA tail of 30 adenines",
                        metavar="[string]",
                        default=30*"A")
    parser.add_argument("--dRNA",
                        dest="dRNA",
                        help="Enable direct RNA (dRNA) mode: skip 5' adapter search; mark the true 5' end as 'correct' based on mapper FLAG strand; keep 3' poly(A)/three_adapter QC as usual; intron logic unchanged.",
                        action="store_true",
                        default=False)
    parser.add_argument("--pacbio",
                        dest="pacbio",
                        help="Enable adapterless PacBio full-length cDNA mode: search only for terminal 3' poly(A), orient from a unique correct poly(A) end, and mark the opposite end as 5' correct.",
                        action="store_true",
                        default=False)
    parser.add_argument("--threads",
                        dest="threads",
                        help="Threads used only for pysam/samtools I/O. Default: 8.",
                        type=int,
                        default=8,
                        metavar="[integer]")
    parser.add_argument("-5", "--five_adapter", 
                        dest="five_adapter",
                        help="The 5' adapter to look for. The default is the\
                        TeloPrime cap adapter: TGGATTGATATGTAATACGACTCACTATAG",
                        metavar="[string]",
                        default="TGGATTGATATGTAATACGACTCACTATAG")
    parser.add_argument("--five_score",
                        dest="five_score",
                        help="5' local-alignment score threshold. The score must be "
                             "strictly greater than this value; equality fails. "
                             "Generic default: 16.0.",
                        type=float,
                        metavar="[float]",
                        default=16.0)
    parser.add_argument("--three_score",
                        dest="three_score",
                        help="3' local-alignment score threshold. The score must be "
                             "strictly greater than this value; equality fails. "
                             "Generic default: 16.0.",
                        type=float,
                        metavar="[float]",
                        default=16.0)
    parser.add_argument("--check_in_soft",
                        dest="check_in_soft",
                        help="Maximum number of terminal soft-clipped query "
                             "nucleotides included in the adapter-search window "
                             "immediately adjacent to the mapped boundary. "
                             "Default: 30 nt.",
                        type=int,
                        metavar="[integer]",
                        default=30)
    parser.add_argument("--check_in_match",
                        dest="check_in_match",
                        help="Number of immediately adjacent aligned read "
                             "nucleotides included from inside the mapped boundary "
                             "in the adapter-search window. These bases need not "
                             "match the reference perfectly; this is not a required "
                             "adapter-match count. Default: 10 nt.",
                        type=int,
                        metavar="[integer]",
                        default=10)
    parser.add_argument("--check_from_alignment",
                        dest="check_from_alignment",
                        help="Maximum number of read bases by which the boundary-"
                             "facing end of a local adapter hit may stop before "
                             "the mapped boundary on the clipped side. This is "
                             "not a genomic distance or a distance from the read "
                             "tip. Default: 3 nt.",
                        type=int,
                        metavar="[integer]",
                        default=3)
    parser.add_argument("--shs_for_ts",
                        dest="shs_for_ts",
                        help="Minimum adapter-associated short-homology length "
                             "used to flag a terminal adapter call as potential "
                             "template switching. The threshold should be smaller "
                             "than --check_in_match; values greater than "
                             "--check_in_match make this read-end test inactive. "
                             "Default: 3 nt.",
                        type=int,
                        metavar="[integer]",
                        default=3)
    parser.add_argument("--first_exon", 
                        dest="first_exon",
                        help="Alignment ends are often placed far away from \
                        from the rest of the read if the adapter maps to a \
                        nearby part of the genome. This option sets the length\
                        of the first exon, under which the matching part of\
                        of the alignment should be checked for the presence of\
                        the adapters. The default is 30.",
                        type=int,
                        metavar="[integer]",
                        default=30)
    parser.add_argument("--match_in_first", 
                        dest="match_in_first",
                        help="Alignment ends are often placed far away from \
                        from the rest of the read if the adapter maps to a \
                        nearby part of the genome. With this option, the user\
                        can set how many nucleotides from the matching part \
                        of the alignment should be aligned to the adapter \
                        if at least half of the nucleotides match to the \
                        adapter, the exon will be considered false. The \
                        default is 15.",
                        type=int,
                        metavar="[integer]",
                        default=15)
    parser.add_argument("--max_split_query_overlap",
                        dest="max_split_query_overlap",
                        help="Maximum query-coordinate overlap (nt) allowed "
                             "between retained split-alignment segments from the "
                             "same read. More-overlapping records are treated as "
                             "alternative mappings and collapsed. Default: 15.",
                        type=int,
                        metavar="[integer]",
                        default=15)
    parser.add_argument("--insert_before_intron", 
                        dest="insert_before_intron",
                        help="The maximum allowed insert length immediately \
                        before an intron. Triple-chimeric reads are often \
                        as an exon, a long insert and another exon. In these \
                        cases the inserts are usually several hundred nts \
                        long and are unmapped because they stem from another \
                        contig or would be mapped to the complementary strand.\
                        The default value is 20 nts.",
                        type=int,
                        metavar="[integer]",
                        default=20)

    return parser.parse_args()



if __name__== "__main__":
    main()

#Explanation of constants for adapter search: 

#S: softclip      M: Match (alignment)

# SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
#    |____________ADAPTER________|
#    -----------------------------: CHECK_IN_SOFT
#    how many nucleotides are checked from the softclip

# SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
#                                 ADAPTER\
#                  CHECK_IN_MATCH:-----
# how many nucleotides are checked from the beginning of the alignment

# SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
#        |______ADAPTER______|
#        CHECK_FROM_ALIGNMENT:---
# how many bases are allowed between the start of the alignment and the adapter

# SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
#                      FIRST_EXON:--------------------------------------
# if an intron starts in this region, the adapter will be classified as
# 'out of place',

