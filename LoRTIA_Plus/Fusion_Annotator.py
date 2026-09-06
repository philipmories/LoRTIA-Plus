#!/usr/bin/env python3
"""LoRTIA Plus chimeric/fusion transcript annotator.

This module consumes the full LoRTIA-processed BAM (including supplementary
segments), groups retained split alignments by QNAME, validates the complete
molecule using accepted LoRTIA TSS/TES/intron features, and reports validated
multi-locus/chimeric transcript chains and their junctions.

It intentionally does not force multi-contig transcripts into the conventional
single-seqid GFF3 representation. The primary TSV outputs are layered:
  <output_prefix>_fusion_candidates.tsv  all structural breakpoint candidates
  <output_prefix>_fusion_molecules.tsv   every split molecule with orientation,
                                         adapter QC and TSS/TES annotation/counts
  <output_prefix>_fusion_transcripts.tsv reconstructed fusion isoforms from
                                         accepted LoRTIA TSS/TES/introns; adapters
                                         are not required
  <output_prefix>_fusion_transcripts_adapter.tsv
                                         the same reconstructed isoform space,
                                         restricted to isoforms with reads carrying
                                         correct biological 5' AND 3' adapters
  <output_prefix>_fusions.tsv            legacy strict LoRTIA-validated chains
BEDPE companions are also written for structural and legacy strict junctions.
"""

from argparse import ArgumentParser
from ast import literal_eval
from bisect import bisect_left
from collections import defaultdict, Counter
from statistics import median
import os
import csv
import pysam

RECONSTRUCTION_QCFIX_VERSION = "2026-08-19-v3-candidate-qc"
CANDIDATE_QC_VERSION = "2026-08-19-v4-query-coverage-ambiguity"


ACCEPTED_ADAPTER_STATES = {"correct"}


def _load_end_index(path):
    """Load accepted TSS/TES positions and their LoRTIA feature counts.

    Returns ``(index, support)`` where ``index`` maps ``(contig, strand)`` to
    sorted accepted positions and ``support`` maps
    ``(contig, strand, position)`` to the GFF score/count.
    """
    index = defaultdict(list)
    support = {}
    with open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 7:
                continue
            try:
                pos = int(fields[3])
            except ValueError:
                continue
            strand = fields[6]
            key = (fields[0], strand)
            index[key].append(pos)
            try:
                raw = float(fields[5])
                count = int(raw) if raw.is_integer() else raw
            except Exception:
                count = None
            support[(fields[0], strand, pos)] = count
    for key in index:
        index[key] = sorted(set(index[key]))
    return dict(index), support


def _load_intron_index(path):
    """Load accepted introns using the same coordinate convention as LoRTIA."""
    introns = set()
    with open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 7:
                continue
            try:
                # GFF intron = SAM/LoRTIA intron boundaries shifted inward by 1.
                left = int(fields[3]) - 1
                right = int(fields[4]) + 1
            except ValueError:
                continue
            introns.add((fields[0], fields[6], left, right))
    return introns


def _nearest_feature(index, contig, strand, position, wobble):
    arr = index.get((contig, strand))
    if not arr:
        return None
    i = bisect_left(arr, position)
    candidates = []
    if i < len(arr):
        candidates.append(arr[i])
    if i > 0:
        candidates.append(arr[i - 1])
    if not candidates:
        return None
    best = min(candidates, key=lambda x: (abs(x - position), x))
    return best if abs(best - position) <= wobble else None


def _parse_introns(value):
    if value in (None, "", "None"):
        return ()
    try:
        parsed = literal_eval(value)
    except Exception:
        return None
    if parsed == ():
        return ()
    if isinstance(parsed, tuple) and len(parsed) == 2 and all(isinstance(x, int) for x in parsed):
        return (tuple(parsed),)
    if isinstance(parsed, (tuple, list)):
        out = []
        for item in parsed:
            if isinstance(item, (tuple, list)) and len(item) == 2:
                try:
                    out.append((int(item[0]), int(item[1])))
                except Exception:
                    return None
        return tuple(out)
    return None


def _get_tag(read, tag, default=None):
    try:
        return read.get_tag(tag)
    except Exception:
        return default


def _adapter_state(read, tag):
    value = _get_tag(read, tag, None)
    if value is None:
        return "missing"
    try:
        return str(value).split(",")[3]
    except Exception:
        return "missing"


def _segment_from_read(read):
    count = int(_get_tag(read, "ZN", 1))
    index = int(_get_tag(read, "ZI", 1))
    direction = str(_get_tag(read, "ZO", "ambiguous"))
    orientation_source = str(_get_tag(read, "ZQ", "legacy_or_unresolved"))
    ts_consensus = str(_get_tag(read, "ZT", "missing"))
    adapter_direction = str(_get_tag(read, "ZP", "ambiguous"))
    ts_tag = str(_get_tag(read, "ts", "."))
    original_flag = int(_get_tag(read, "ZF", read.flag))
    original_reverse = bool(original_flag & 0x10)
    original_strand = "-" if original_reverse else "+"

    ref_start = int(read.reference_start) + 1
    ref_end = int(read.reference_end)
    if original_reverse:
        query_left_ref = ref_end
        query_right_ref = ref_start
    else:
        query_left_ref = ref_start
        query_right_ref = ref_end

    introns = _parse_introns(_get_tag(read, "in", ""))
    ga = _get_tag(read, "ga", None)
    gap_flag = ga not in (None, "", "None")

    return {
        "read_name": read.query_name,
        "segment_index": index,
        "segment_count": count,
        "direction": direction,
        "orientation_source": orientation_source,
        "ts_consensus": ts_consensus,
        "adapter_direction": adapter_direction,
        "ts_tag": ts_tag,
        "internal_adapter": int(_get_tag(read, "ZA", 0)),
        "full_sequence_recovered": int(_get_tag(read, "ZU", 0)),
        "qstart": int(_get_tag(read, "ZS", 0)),
        "qend": int(_get_tag(read, "ZE", 0)),
        "qlen": int(_get_tag(read, "ZL", 0)),
        # QC metadata from Samprocessor; backward-compatible defaults keep older BAMs usable.
        "query_coverage": float(_get_tag(read, "ZC", 0.0)),
        "alternative_mapping_count": int(_get_tag(read, "ZX", 0)),
        "secondary_alignment_count": int(_get_tag(read, "ZY", 0)),
        "mapping_ambiguous": int(_get_tag(read, "ZV", 0)),
        "contig": read.reference_name,
        "ref_start": ref_start,
        "ref_end": ref_end,
        "strand": ".",
        "original_strand": original_strand,
        "query_left_ref": query_left_ref,
        "query_right_ref": query_right_ref,
        "original_flag": original_flag,
        "mapq": int(read.mapping_quality),
        "introns": introns,
        "gap_flag": gap_flag,
        "l5": _adapter_state(read, "l5"),
        "r5": _adapter_state(read, "r5"),
        "l3": _adapter_state(read, "l3"),
        "r3": _adapter_state(read, "r3"),
        "bio5": None,
        "bio3": None,
    }


def _resolve_group_orientation(segments):
    """Resolve split-molecule orientation, preferring minimap2 ``ts`` consensus."""
    ts_values = {s["ts_tag"] for s in segments if s.get("ts_tag") in {"+", "-"}}
    if ts_values == {"+"}:
        return "5to3", "minimap2_ts", "+"
    if ts_values == {"-"}:
        return "3to5", "minimap2_ts", "-"

    adapter_values = {s.get("adapter_direction") for s in segments
                      if s.get("adapter_direction") in {"5to3", "3to5"}}
    # Backward-compatible fallback: derive external adapter orientation directly
    # when ZP was not written by an older Samprocessor.
    if len(adapter_values) != 1 and segments:
        ordered = sorted(segments, key=lambda s: s["segment_index"])
        left, right = ordered[0], ordered[-1]
        left5 = _state_on_query_side(left, "left", "5")
        left3 = _state_on_query_side(left, "left", "3")
        right5 = _state_on_query_side(right, "right", "5")
        right3 = _state_on_query_side(right, "right", "3")
        forward = left5 in ACCEPTED_ADAPTER_STATES and right3 in ACCEPTED_ADAPTER_STATES
        reverse = left3 in ACCEPTED_ADAPTER_STATES and right5 in ACCEPTED_ADAPTER_STATES
        if forward and not reverse:
            adapter_values = {"5to3"}
        elif reverse and not forward:
            adapter_values = {"3to5"}

    if len(adapter_values) == 1:
        direction = next(iter(adapter_values))
        source = "adapter_fallback_conflicting_ts" if ts_values == {"+", "-"} else "adapter_fallback_no_ts"
        ts_state = "conflict" if ts_values == {"+", "-"} else "missing"
        return direction, source, ts_state

    if ts_values == {"+", "-"}:
        return "ambiguous", "conflicting_ts", "conflict"
    return "ambiguous", "unresolved", "missing"


def _apply_direction(segments, direction):
    """Set transcript genomic strand and biological 5'/3' coordinates."""
    for s in segments:
        original_reverse = s["original_strand"] == "-"
        if direction == "5to3":
            strand = "-" if original_reverse else "+"
        elif direction == "3to5":
            strand = "+" if original_reverse else "-"
        else:
            strand = "."
        s["strand"] = strand
        if strand == "+":
            s["bio5"], s["bio3"] = s["ref_start"], s["ref_end"]
        elif strand == "-":
            s["bio5"], s["bio3"] = s["ref_end"], s["ref_start"]
        else:
            s["bio5"] = s["bio3"] = None
    return segments


def _state_on_query_side(segment, side, kind):
    if side == "left":
        lr = "r" if segment["original_strand"] == "-" else "l"
    else:
        lr = "l" if segment["original_strand"] == "-" else "r"
    return segment.get(lr + kind, "missing")


def _terminal_adapter_evidence(segments, direction):
    ordered = sorted(segments, key=lambda s: s["segment_index"])
    left, right = ordered[0], ordered[-1]
    evidence = {
        "left5": _state_on_query_side(left, "left", "5"),
        "left3": _state_on_query_side(left, "left", "3"),
        "right5": _state_on_query_side(right, "right", "5"),
        "right3": _state_on_query_side(right, "right", "3"),
    }
    forward = evidence["left5"] in ACCEPTED_ADAPTER_STATES and evidence["right3"] in ACCEPTED_ADAPTER_STATES
    reverse = evidence["left3"] in ACCEPTED_ADAPTER_STATES and evidence["right5"] in ACCEPTED_ADAPTER_STATES
    if forward and not reverse:
        adapter_direction = "5to3"
    elif reverse and not forward:
        adapter_direction = "3to5"
    else:
        adapter_direction = "ambiguous"
    if direction == "5to3":
        five_valid = evidence["left5"] in ACCEPTED_ADAPTER_STATES
        three_valid = evidence["right3"] in ACCEPTED_ADAPTER_STATES
        valid = five_valid and three_valid
    elif direction == "3to5":
        five_valid = evidence["right5"] in ACCEPTED_ADAPTER_STATES
        three_valid = evidence["left3"] in ACCEPTED_ADAPTER_STATES
        valid = five_valid and three_valid
    else:
        five_valid = False
        three_valid = False
        valid = False
    evidence["adapter_direction_observed"] = adapter_direction
    evidence["five_prime_adapter_valid"] = int(five_valid)
    evidence["three_prime_adapter_valid"] = int(three_valid)
    evidence["terminal_adapter_valid"] = int(valid)
    return evidence


def _validate_introns(segment, intron_index):
    if segment["introns"] is None:
        return False, "malformed_intron_tag"
    for left, right in segment["introns"]:
        direct = (segment["contig"], segment["strand"], left, right)
        unstranded = (segment["contig"], ".", left, right)
        if direct not in intron_index and unstranded not in intron_index:
            return False, "unaccepted_intron"
    return True, ""


def _junction_class(a, b, local_distance):
    if a["contig"] != b["contig"]:
        return "interchromosomal"
    if a["strand"] != b["strand"]:
        return "intrachromosomal_strand_switch"
    donor = a["bio3"]
    acceptor = b["bio5"]
    if a["strand"] == "+":
        collinear = acceptor >= donor
    else:
        collinear = acceptor <= donor
    if not collinear:
        return "intrachromosomal_rearranged"
    if abs(acceptor - donor) <= local_distance:
        return "intrachromosomal_local"
    return "intrachromosomal_distant"


def _annotate_group(read_name, segments, tss_index, tes_index, tss_support, tes_support, intron_index, args):
    """Annotate one split molecule without discarding failed validation states."""
    reasons = []
    segments = sorted((dict(s) for s in segments), key=lambda s: s["segment_index"])
    if len(segments) < 2:
        reasons.append("single_segment")

    expected = max((s["segment_count"] for s in segments), default=len(segments))
    if len(segments) != expected:
        reasons.append("incomplete_split_chain")

    query_gaps = [right["qstart"] - left["qend"] for left, right in zip(segments, segments[1:])]
    if any(gap > int(args.fusion_max_query_gap) for gap in query_gaps):
        reasons.append("large_unaligned_query_gap")

    direction, orientation_source, ts_consensus = _resolve_group_orientation(segments)
    _apply_direction(segments, direction)
    adapter_evidence = _terminal_adapter_evidence(segments, direction)

    if direction not in {"5to3", "3to5"}:
        reasons.append("unresolved_molecule_orientation")

    if not adapter_evidence["terminal_adapter_valid"]:
        reasons.append("terminal_adapter_validation_failed")

    internal_adapter = int(any(s["internal_adapter"] for s in segments))
    if internal_adapter and not args.fusion_allow_internal_adapters:
        reasons.append("internal_adapter_detected")

    if any(s["gap_flag"] for s in segments):
        reasons.append("gap_or_chimera_artifact_tag")

    biological = list(segments)
    if direction == "3to5":
        biological.reverse()

    introns_valid = True
    intron_reason = "."
    if direction in {"5to3", "3to5"}:
        for segment in biological:
            ok, reason = _validate_introns(segment, intron_index)
            if not ok:
                introns_valid = False
                intron_reason = reason
                reasons.append(reason)
                break
    else:
        introns_valid = False
        intron_reason = "orientation_required"

    tss_raw = tes_raw = None
    tss_match = tes_match = None
    tss_count = tes_count = None
    tss_contig = tes_contig = "."
    tss_strand = tes_strand = "."
    if direction in {"5to3", "3to5"} and biological:
        first, last = biological[0], biological[-1]
        tss_contig, tss_strand, tss_raw = first["contig"], first["strand"], first["bio5"]
        tes_contig, tes_strand, tes_raw = last["contig"], last["strand"], last["bio3"]
        tss_match = _nearest_feature(tss_index, tss_contig, tss_strand, tss_raw, args.wobble)
        tes_match = _nearest_feature(tes_index, tes_contig, tes_strand, tes_raw, args.wobble)
        if tss_match is not None:
            tss_count = tss_support.get((tss_contig, tss_strand, int(tss_match)))
        else:
            reasons.append("no_accepted_tss")
        if tes_match is not None:
            tes_count = tes_support.get((tes_contig, tes_strand, int(tes_match)))
        else:
            reasons.append("no_accepted_tes")

    classes = []
    if direction in {"5to3", "3to5"}:
        for a, b in zip(biological, biological[1:]):
            classes.append(_junction_class(a, b, args.fusion_local_distance))
    else:
        for a, b in zip(segments, segments[1:]):
            classes.append("interchromosomal" if a["contig"] != b["contig"] else "intrachromosomal")

    validation_reasons = sorted(set(reasons))
    validated = len(validation_reasons) == 0
    annotation = {
        "read_name": read_name,
        "direction": direction,
        "orientation_source": orientation_source,
        "ts_consensus": ts_consensus,
        "adapter_direction": adapter_evidence["adapter_direction_observed"],
        "five_prime_adapter_valid": adapter_evidence["five_prime_adapter_valid"],
        "three_prime_adapter_valid": adapter_evidence["three_prime_adapter_valid"],
        "terminal_adapter_valid": adapter_evidence["terminal_adapter_valid"],
        "left5": adapter_evidence["left5"],
        "left3": adapter_evidence["left3"],
        "right5": adapter_evidence["right5"],
        "right3": adapter_evidence["right3"],
        "internal_adapter": internal_adapter,
        "segments": biological if direction in {"5to3", "3to5"} else segments,
        "query_gaps": tuple(query_gaps),
        "tss_contig": tss_contig,
        "tss_strand": tss_strand,
        "tss_raw": tss_raw,
        "tss": int(tss_match) if tss_match is not None else None,
        "tss_feature_count": tss_count,
        "tes_contig": tes_contig,
        "tes_strand": tes_strand,
        "tes_raw": tes_raw,
        "tes": int(tes_match) if tes_match is not None else None,
        "tes_feature_count": tes_count,
        "introns_valid": int(introns_valid),
        "intron_reason": intron_reason,
        "classes": tuple(classes),
        "class": ",".join(sorted(set(classes))) if classes else ".",
        "full_sequence_recovered": int(any(s["full_sequence_recovered"] for s in segments)),
        "validated": int(validated),
        "validation_reasons": validation_reasons,
    }
    return annotation


def _candidate_from_group(read_name, segments, tss_index, tes_index, tss_support, tes_support, intron_index, args):
    annotated = _annotate_group(
        read_name, segments, tss_index, tes_index, tss_support, tes_support, intron_index, args
    )
    if annotated["validated"]:
        return annotated, []
    return None, annotated["validation_reasons"]


def _same_chain(a, b, wobble):
    if a["direction"] != b["direction"]:
        return False
    if a["tss"] != b["tss"] or a["tes"] != b["tes"]:
        return False
    if len(a["segments"]) != len(b["segments"]):
        return False
    for sa, sb in zip(a["segments"], b["segments"]):
        if sa["contig"] != sb["contig"] or sa["strand"] != sb["strand"]:
            return False
        if sa["introns"] != sb["introns"]:
            return False
        if abs(sa["bio5"] - sb["bio5"]) > wobble:
            return False
        if abs(sa["bio3"] - sb["bio3"]) > wobble:
            return False
    return True


def _cluster_candidates(candidates, wobble):
    clusters = []
    for cand in candidates:
        target = None
        for cluster in clusters:
            if _same_chain(cand, cluster["representative"], wobble):
                target = cluster
                break
        if target is None:
            clusters.append({"representative": cand, "members": [cand]})
        else:
            target["members"].append(cand)
    return clusters


def _median_int(values):
    return int(round(median(values)))


def _cluster_segments(cluster):
    members = cluster["members"]
    n = len(members[0]["segments"])
    out = []
    for i in range(n):
        base = dict(members[0]["segments"][i])
        base["bio5"] = _median_int([m["segments"][i]["bio5"] for m in members])
        base["bio3"] = _median_int([m["segments"][i]["bio3"] for m in members])
        out.append(base)
    # Accepted outer features are stronger representatives than raw alignment ends.
    out[0]["bio5"] = members[0]["tss"]
    out[-1]["bio3"] = members[0]["tes"]
    return out


def _segment_text(segment):
    intron_text = ",".join(f"{a}-{b}" for a, b in segment["introns"]) if segment["introns"] else "."
    return (
        f'{segment["contig"]}:{segment["strand"]}:'
        f'{segment["bio5"]}>{segment["bio3"]}:introns={intron_text}'
    )


def _junction_text(a, b):
    return f'{a["contig"]}:{a["bio3"]}:{a["strand"]}>{b["contig"]}:{b["bio5"]}:{b["strand"]}'




def _reconstruction_info(annotation):
    """Return the adapter-independent LoRTIA fusion reconstruction status.

    This is intentionally independent of the legacy strict fusion validation.
    A read is reconstructable when the biological orientation is resolved, the
    retained split chain is complete, an accepted LoRTIA TSS and TES are
    available, and every intron in every retained genomic segment is accepted.

    Adapter absence, internal-adapter evidence, large query gaps and generic
    gap/chimera QC tags are preserved as QC annotations but do NOT block the
    adapter-independent transcript reconstruction.  The separate
    *_fusion_transcripts_adapter.tsv output supplies the stricter subset with
    both external terminal adapters.
    """
    five_adapter = bool(annotation.get("five_prime_adapter_valid", 0))
    three_adapter = bool(annotation.get("three_prime_adapter_valid", 0))

    tss_pos = int(annotation["tss"]) if annotation.get("tss") is not None else None
    tes_pos = int(annotation["tes"]) if annotation.get("tes") is not None else None
    tss_source = (
        "accepted_tss+adapter" if tss_pos is not None and five_adapter
        else "accepted_tss" if tss_pos is not None
        else "missing"
    )
    tes_source = (
        "accepted_tes+adapter" if tes_pos is not None and three_adapter
        else "accepted_tes" if tes_pos is not None
        else "missing"
    )

    # IMPORTANT: do not inherit the complete legacy validation_reasons list.
    # Reconstruction has its own biological requirements.
    reasons = set()

    direction = annotation.get("direction")
    if direction not in {"5to3", "3to5"}:
        reasons.add("unresolved_molecule_orientation")

    segments = annotation.get("segments") or ()
    if len(segments) < 2:
        reasons.add("single_segment")
    else:
        # A missing retained subalignment means the full fusion chain itself is
        # incomplete; unlike generic gap/QC flags, this remains a blocker.
        expected = max((int(s.get("segment_count", len(segments))) for s in segments),
                       default=len(segments))
        if len(segments) != expected:
            reasons.add("incomplete_split_chain")

    if tss_pos is None:
        reasons.add("no_accepted_tss")
    if tes_pos is None:
        reasons.add("no_accepted_tes")

    if not bool(annotation.get("introns_valid", 0)):
        intron_reason = annotation.get("intron_reason")
        if intron_reason and intron_reason not in {".", ""}:
            reasons.add(str(intron_reason))
        else:
            reasons.add("invalid_intron_chain")

    forbidden_legacy_qc = {
        "terminal_adapter_validation_failed",
        "internal_adapter_detected",
        "large_unaligned_query_gap",
        "gap_or_chimera_artifact_tag",
    }
    leaked = reasons & forbidden_legacy_qc
    if leaked:
        raise RuntimeError(
            "Legacy QC reason leaked into adapter-independent reconstruction: "
            + ";".join(sorted(leaked))
        )

    return {
        "reconstructable": int(len(reasons) == 0),
        "reconstruction_reasons": tuple(sorted(reasons)),
        "tss_effective": tss_pos,
        "tes_effective": tes_pos,
        "tss_source": tss_source,
        "tes_source": tes_source,
    }


def _same_transcript_chain(a, b, end_wobble, junction_wobble):
    """Compare biological fusion transcript isoforms, ignoring read direction."""
    ra, rb = a["reconstruction"], b["reconstruction"]
    if not ra["reconstructable"] or not rb["reconstructable"]:
        return False
    if a["tss_contig"] != b["tss_contig"] or a["tss_strand"] != b["tss_strand"]:
        return False
    if a["tes_contig"] != b["tes_contig"] or a["tes_strand"] != b["tes_strand"]:
        return False
    if abs(int(ra["tss_effective"]) - int(rb["tss_effective"])) > end_wobble:
        return False
    if abs(int(ra["tes_effective"]) - int(rb["tes_effective"])) > end_wobble:
        return False
    if len(a["segments"]) != len(b["segments"]):
        return False

    n = len(a["segments"])
    for i, (sa, sb) in enumerate(zip(a["segments"], b["segments"])):
        if sa["contig"] != sb["contig"] or sa["strand"] != sb["strand"]:
            return False
        if sa["introns"] != sb["introns"]:
            return False
        # Outer transcript ends are represented by TSS/TES above.  Only the
        # internal segment boundaries define fusion junction identity here.
        if i > 0 and abs(int(sa["bio5"]) - int(sb["bio5"])) > junction_wobble:
            return False
        if i < n - 1 and abs(int(sa["bio3"]) - int(sb["bio3"])) > junction_wobble:
            return False
    return True


def _cluster_reconstructed_transcripts(candidates, end_wobble, junction_wobble):
    clusters = []
    for cand in candidates:
        target = None
        for cluster in clusters:
            if _same_transcript_chain(cand, cluster["representative"], end_wobble, junction_wobble):
                target = cluster
                break
        if target is None:
            clusters.append({"representative": cand, "members": [cand]})
        else:
            target["members"].append(cand)
    return clusters


def _canonical_outer_end(members, kind):
    """Prefer an accepted statistical LoRTIA end; otherwise use adapter end median."""
    accepted_key = "tss" if kind == "tss" else "tes"
    effective_key = "tss_effective" if kind == "tss" else "tes_effective"
    accepted = [int(m[accepted_key]) for m in members if m.get(accepted_key) is not None]
    if accepted:
        counts = Counter(accepted)
        max_count = max(counts.values())
        tied = sorted(pos for pos, n in counts.items() if n == max_count)
        pos = tied[len(tied) // 2]
        source = "accepted_feature"
    else:
        pos = _median_int([int(m["reconstruction"][effective_key]) for m in members])
        source = "adapter"
    return pos, source


def _segment_exons(segment):
    """Reconstruct exon intervals for one validated genomic segment."""
    low = min(int(segment["bio5"]), int(segment["bio3"]))
    high = max(int(segment["bio5"]), int(segment["bio3"]))
    introns = sorted(tuple(map(int, x)) for x in (segment.get("introns") or ()))

    exons = []
    cursor = low
    for left, right in introns:
        if right < low or left > high:
            continue
        left = max(left, low)
        right = min(right, high)
        if cursor <= left:
            exons.append((cursor, left))
        cursor = max(cursor, right)
    if cursor <= high:
        exons.append((cursor, high))

    if segment["strand"] == "-":
        exons.reverse()
    return tuple(exons)


def _canonical_reconstructed_segments(cluster, tss_pos, tes_pos):
    members = cluster["members"]
    n = len(members[0]["segments"])
    out = []
    for i in range(n):
        base = dict(members[0]["segments"][i])
        base["bio5"] = _median_int([m["segments"][i]["bio5"] for m in members])
        base["bio3"] = _median_int([m["segments"][i]["bio3"] for m in members])
        out.append(base)
    out[0]["bio5"] = int(tss_pos)
    out[-1]["bio3"] = int(tes_pos)
    return out


def _source_summary(values):
    c = Counter(values)
    return ",".join(f"{k}:{c[k]}" for k in sorted(c))


def _write_reconstructed_transcripts(annotations, args):
    """Write adapter-optional and both-terminal-adapter fusion transcript isoforms."""
    candidates = []
    for ann in annotations.values():
        ann = dict(ann)
        ann["reconstruction"] = _reconstruction_info(ann)
        if ann["reconstruction"]["reconstructable"]:
            candidates.append(ann)

    clusters = _cluster_reconstructed_transcripts(
        candidates, int(args.wobble), int(args.fusion_wobble)
    )

    all_path = args.output_prefix + "_fusion_transcripts.tsv"
    adapter_path = args.output_prefix + "_fusion_transcripts_adapter.tsv"
    fields = [
        "fusion_transcript_id", "support_total", "support_5_adapter",
        "support_3_adapter", "support_both_adapters", "segment_count", "class",
        "tss_contig", "tss_strand", "tss_pos", "tss_source",
        "tss_feature_count", "tss_statistical_read_support",
        "tes_contig", "tes_strand", "tes_pos", "tes_source",
        "tes_feature_count", "tes_statistical_read_support",
        "orientation_sources", "ts_consensus",
        "support_terminal_adapter_failed", "support_internal_adapter",
        "support_large_unaligned_query_gap", "support_gap_or_chimera_artifact",
        "legacy_qc_reason_summary",
        "read_names", "adapter_read_names",
        "segments", "exon_chain", "intron_chain", "fusion_junctions",
    ]

    all_rows = []
    adapter_rows = []
    tr_no = 0
    min_reads = int(args.fusion_min_reads)

    for cluster in clusters:
        members = cluster["members"]
        read_names = sorted({m["read_name"] for m in members})
        support_total = len(read_names)
        if support_total < min_reads:
            continue

        adapter_members = [m for m in members if m.get("terminal_adapter_valid")]
        adapter_read_names = sorted({m["read_name"] for m in adapter_members})
        support_both = len(adapter_read_names)
        support_5 = len({m["read_name"] for m in members if m.get("five_prime_adapter_valid")})
        support_3 = len({m["read_name"] for m in members if m.get("three_prime_adapter_valid")})

        tss_pos, tss_source = _canonical_outer_end(members, "tss")
        tes_pos, tes_source = _canonical_outer_end(members, "tes")
        segments = _canonical_reconstructed_segments(cluster, tss_pos, tes_pos)
        exons = [_segment_exons(s) for s in segments]
        classes = [_junction_class(a, b, args.fusion_local_distance)
                   for a, b in zip(segments, segments[1:])]
        chain_class = ",".join(sorted(set(classes))) if classes else "."

        tss_feature_count = "."
        for m in members:
            if m.get("tss") == tss_pos and m.get("tss_feature_count") is not None:
                tss_feature_count = m["tss_feature_count"]
                break
        tes_feature_count = "."
        for m in members:
            if m.get("tes") == tes_pos and m.get("tes_feature_count") is not None:
                tes_feature_count = m["tes_feature_count"]
                break

        # Legacy strict-QC reasons are retained for auditability, but the
        # adapter-independent reconstruction above does not inherit them as
        # exclusion criteria (except incomplete split chain, handled earlier).
        legacy_reason_counts = Counter(
            reason
            for m in members
            for reason in (m.get("validation_reasons") or ())
        )
        legacy_qc_reason_summary = (
            ";".join(f"{reason}:{legacy_reason_counts[reason]}"
                     for reason in sorted(legacy_reason_counts))
            if legacy_reason_counts else "."
        )

        tr_no += 1
        row = {
            "fusion_transcript_id": f"fusion_tr{tr_no}",
            "support_total": support_total,
            "support_5_adapter": support_5,
            "support_3_adapter": support_3,
            "support_both_adapters": support_both,
            "segment_count": len(segments),
            "class": chain_class,
            "tss_contig": segments[0]["contig"],
            "tss_strand": segments[0]["strand"],
            "tss_pos": tss_pos,
            "tss_source": tss_source,
            "tss_feature_count": tss_feature_count,
            "tss_statistical_read_support": len({m["read_name"] for m in members if m.get("tss") is not None}),
            "tes_contig": segments[-1]["contig"],
            "tes_strand": segments[-1]["strand"],
            "tes_pos": tes_pos,
            "tes_source": tes_source,
            "tes_feature_count": tes_feature_count,
            "tes_statistical_read_support": len({m["read_name"] for m in members if m.get("tes") is not None}),
            "orientation_sources": _source_summary([m["orientation_source"] for m in members]),
            "ts_consensus": _source_summary([m["ts_consensus"] for m in members]),
            "support_terminal_adapter_failed": sum(
                "terminal_adapter_validation_failed" in (m.get("validation_reasons") or ())
                for m in members
            ),
            "support_internal_adapter": sum(
                "internal_adapter_detected" in (m.get("validation_reasons") or ())
                for m in members
            ),
            "support_large_unaligned_query_gap": sum(
                "large_unaligned_query_gap" in (m.get("validation_reasons") or ())
                for m in members
            ),
            "support_gap_or_chimera_artifact": sum(
                "gap_or_chimera_artifact_tag" in (m.get("validation_reasons") or ())
                for m in members
            ),
            "legacy_qc_reason_summary": legacy_qc_reason_summary,
            "read_names": ",".join(read_names),
            "adapter_read_names": ",".join(adapter_read_names) if adapter_read_names else ".",
            "segments": ";".join(_segment_text(s) for s in segments),
            "exon_chain": ";".join(
                f'{s["contig"]}:{s["strand"]}:' +
                ",".join(f"{start}-{end}" for start, end in seg_exons)
                for s, seg_exons in zip(segments, exons)
            ),
            "intron_chain": ";".join(
                f'{s["contig"]}:{s["strand"]}:' +
                (",".join(f"{a}-{b}" for a, b in (s.get("introns") or ())) or ".")
                for s in segments
            ),
            "fusion_junctions": ";".join(
                _junction_text(a, b) for a, b in zip(segments, segments[1:])
            ),
        }
        all_rows.append(row)
        # The adapter file is reconstructed from the same biological isoform
        # definition but requires >= fusion_min_reads reads with BOTH expected
        # external adapters.  Soft clipping itself is never a requirement.
        if support_both >= min_reads:
            adapter_rows.append(row)

    for path, rows in ((all_path, all_rows), (adapter_path, adapter_rows)):
        with open(path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    return len(candidates), len(all_rows), len(adapter_rows)


def _canonical_breakpoint_pair(a, b):
    """Return an orientation-independent breakpoint pair for event clustering.

    The structural candidate layer deliberately does not infer transcript
    5'->3' direction.  That inference belongs to the LoRTIA-validated layer and
    may be unavailable in external benchmark libraries whose terminal adapters
    are not LoRTIA-compatible.
    """
    e1 = (str(a["contig"]), int(a["query_right_ref"]))
    e2 = (str(b["contig"]), int(b["query_left_ref"]))
    if e2 < e1:
        e1, e2 = e2, e1
    return e1, e2


def _structural_junction_class(left, right, local_distance):
    """Classify a split junction without requiring biological transcript orientation.

    Query order and the original mapper strand are sufficient to distinguish
    same-contig collinear splits from strand switches and rearranged splits.
    """
    if left["contig"] != right["contig"]:
        return "interchromosomal"

    left_bp = int(left["query_right_ref"])
    right_bp = int(right["query_left_ref"])
    genomic_distance = abs(right_bp - left_bp)

    # Opposite mapper strands on the same contig are not, by themselves,
    # evidence of a biological strand-switch fusion. Very local opposite-strand
    # split alignments are treated as mapping/split-structure audit events and
    # excluded from the structural fusion-candidate layer by default.
    if left["original_strand"] != right["original_strand"]:
        if genomic_distance <= int(local_distance):
            return "intrachromosomal_local_strand_switch"
        return "intrachromosomal_strand_switch"

    if left["original_strand"] == "+":
        collinear = right_bp >= left_bp
    else:
        collinear = right_bp <= left_bp

    if not collinear:
        return "intrachromosomal_rearranged"
    if abs(right_bp - left_bp) <= int(local_distance):
        return "intrachromosomal_local_collinear"
    return "intrachromosomal_distant_collinear"


def _structural_junctions(read_name, segments, args):
    """Build adapter-independent split-alignment junctions for one molecule.

    Only properties of the split alignment chain are used:
      * at least two retained query segments,
      * complete selected chain,
      * adjacent query segments separated by no more than
        ``fusion_max_query_gap``.

    TSS/TES, adapter orientation, intron validation and LoRTIA artifact tags
    are intentionally NOT used here.  They remain part of the stricter
    validated-fusion output.
    """
    if len(segments) < 2:
        return []

    segments = sorted(segments, key=lambda s: s["segment_index"])
    expected = max(s["segment_count"] for s in segments)
    if len(segments) != expected:
        return []

    for left, right in zip(segments, segments[1:]):
        if (right["qstart"] - left["qend"]) > int(args.fusion_max_query_gap):
            return []

    out = []
    for j, (left, right) in enumerate(zip(segments, segments[1:]), start=1):
        jclass = _structural_junction_class(
            left, right, getattr(args, "fusion_local_filter_distance", 1000)
        )
        if (
            jclass in {
                "intrachromosomal_local_collinear",
                "intrachromosomal_local_strand_switch",
            }
            and getattr(args, "fusion_exclude_local_collinear", True)
        ):
            continue
        (chrom1, pos1), (chrom2, pos2) = _canonical_breakpoint_pair(left, right)
        out.append({
            "read_name": read_name,
            "junction_index": j,
            "chrom1": chrom1,
            "pos1": int(pos1),
            "chrom2": chrom2,
            "pos2": int(pos2),
            "class": jclass,
            "segment_count": len(segments),
            "query_gap": int(right["qstart"] - left["qend"]),
            # Molecule-level QC is repeated on each junction observation and later
            # collapsed to one value per QNAME before candidate-level summaries.
            "query_coverage": float(segments[0].get("query_coverage", 0.0)),
            "alternative_mapping_count": int(segments[0].get("alternative_mapping_count", 0)),
            "secondary_alignment_count": int(segments[0].get("secondary_alignment_count", 0)),
            "mapping_ambiguous": int(segments[0].get("mapping_ambiguous", 0)),
            "original_query_pair": (
                f'{left["contig"]}:{left["query_right_ref"]}:{left["original_strand"]}>'
                f'{right["contig"]}:{right["query_left_ref"]}:{right["original_strand"]}'
            ),
        })
    return out


def _median_abs_deviation(values):
    """Median absolute deviation around the median, returned as a float."""
    vals = [float(v) for v in values]
    if not vals:
        return 0.0
    centre = median(vals)
    return float(median([abs(v - centre) for v in vals]))


def _cluster_read_representatives(cluster):
    """Collapse a structural cluster to one robust observation per supporting read.

    A multi-junction molecule can occasionally contribute more than one nearby
    observation to the same breakpoint cluster. Candidate support is defined by
    distinct QNAMEs, so QC spread statistics should weight each physical read once.
    """
    by_read = defaultdict(list)
    for member in cluster["members"]:
        by_read[member["read_name"]].append(member)

    out = []
    for read_name, members in sorted(by_read.items()):
        out.append({
            "read_name": read_name,
            "pos1": _median_int([m["pos1"] for m in members]),
            "pos2": _median_int([m["pos2"] for m in members]),
            "query_gap": _median_int([m["query_gap"] for m in members]),
            "query_coverage": float(median([float(m.get("query_coverage", 0.0)) for m in members])),
            "alternative_mapping_count": max(int(m.get("alternative_mapping_count", 0)) for m in members),
            "secondary_alignment_count": max(int(m.get("secondary_alignment_count", 0)) for m in members),
            "mapping_ambiguous": max(int(m.get("mapping_ambiguous", 0)) for m in members),
        })
    return out


def _candidate_qc_stats(cluster):
    """Return robust, read-weighted breakpoint and query-gap QC statistics."""
    reps = _cluster_read_representatives(cluster)
    p1 = [r["pos1"] for r in reps]
    p2 = [r["pos2"] for r in reps]
    gaps = [r["query_gap"] for r in reps]
    coverages = [float(r.get("query_coverage", 0.0)) for r in reps]
    alt_counts = [int(r.get("alternative_mapping_count", 0)) for r in reps]
    secondary_counts = [int(r.get("secondary_alignment_count", 0)) for r in reps]
    ambiguous = [int(r.get("mapping_ambiguous", 0)) for r in reps]

    def _span(vals):
        return int(max(vals) - min(vals)) if vals else 0

    return {
        "junction_observations": len(cluster["members"]),
        "breakpoint1_mad": _median_abs_deviation(p1),
        "breakpoint2_mad": _median_abs_deviation(p2),
        "breakpoint1_span": _span(p1),
        "breakpoint2_span": _span(p2),
        "query_gap_median": _median_int(gaps) if gaps else 0,
        "query_gap_mad": _median_abs_deviation(gaps),
        "query_gap_min": int(min(gaps)) if gaps else 0,
        "query_gap_max": int(max(gaps)) if gaps else 0,
        "query_coverage_median": float(median(coverages)) if coverages else 0.0,
        "query_coverage_min": float(min(coverages)) if coverages else 0.0,
        "alternative_mapping_count_median": float(median(alt_counts)) if alt_counts else 0.0,
        "alternative_mapping_count_max": int(max(alt_counts)) if alt_counts else 0,
        "secondary_alignment_read_count": sum(1 for x in secondary_counts if x > 0),
        "ambiguous_read_count": sum(1 for x in ambiguous if x > 0),
        "ambiguous_read_fraction": (
            float(sum(1 for x in ambiguous if x > 0)) / float(len(ambiguous))
            if ambiguous else 0.0
        ),
    }


def _cluster_structural_junctions(junctions, wobble):
    """Greedily cluster structural breakpoint observations across reads."""
    clusters = []
    for junc in sorted(
        junctions,
        key=lambda x: (x["chrom1"], x["chrom2"], x["pos1"], x["pos2"], x["read_name"])
    ):
        target = None
        for cluster in clusters:
            rep = cluster["representative"]
            if rep["chrom1"] != junc["chrom1"] or rep["chrom2"] != junc["chrom2"]:
                continue
            if rep.get("class") != junc.get("class"):
                continue
            if abs(rep["pos1"] - junc["pos1"]) <= wobble and abs(rep["pos2"] - junc["pos2"]) <= wobble:
                target = cluster
                break
        if target is None:
            clusters.append({"representative": dict(junc), "members": [junc]})
        else:
            target["members"].append(junc)
            # Keep the representative near the centre as the cluster grows.
            target["representative"]["pos1"] = _median_int([m["pos1"] for m in target["members"]])
            target["representative"]["pos2"] = _median_int([m["pos2"] for m in target["members"]])
    return clusters


def _write_structural_candidates(groups, args):
    """Write adapter-independent read- and event-level fusion breakpoint output."""
    observations = []
    for read_name, segments in groups.items():
        observations.extend(_structural_junctions(read_name, segments, args))

    read_bedpe = args.output_prefix + "_fusion_candidate_junction_reads.bedpe"
    with open(read_bedpe, "w", newline="", encoding="utf-8") as handle:
        w = csv.writer(handle, delimiter="\t", lineterminator="\n")
        w.writerow([
            "chrom1", "start1", "end1", "chrom2", "start2", "end2",
            "read_name", "score", "strand1", "strand2", "class",
            "junction_index", "segment_count", "query_gap",
            "query_coverage", "alternative_mapping_count",
            "secondary_alignment_count", "mapping_ambiguous", "query_order"
        ])
        for j in observations:
            w.writerow([
                j["chrom1"], max(0, j["pos1"] - 1), j["pos1"],
                j["chrom2"], max(0, j["pos2"] - 1), j["pos2"],
                j["read_name"], 1, ".", ".", j["class"],
                j["junction_index"], j["segment_count"], j["query_gap"],
                f'{j.get("query_coverage", 0.0):.6f}',
                j.get("alternative_mapping_count", 0),
                j.get("secondary_alignment_count", 0),
                j.get("mapping_ambiguous", 0),
                j["original_query_pair"],
            ])

    clusters = _cluster_structural_junctions(observations, int(args.fusion_wobble))

    event_tsv = args.output_prefix + "_fusion_candidates.tsv"
    event_bedpe = args.output_prefix + "_fusion_candidate_junctions.bedpe"
    with open(event_tsv, "w", newline="", encoding="utf-8") as tout, \
         open(event_bedpe, "w", newline="", encoding="utf-8") as bout:
        fields = [
            "candidate_id", "support", "class",
            "chrom1", "breakpoint1", "chrom2", "breakpoint2",
            "junction_observations",
            "breakpoint1_mad", "breakpoint2_mad",
            "breakpoint1_span", "breakpoint2_span",
            "query_gap_median", "query_gap_mad",
            "query_gap_min", "query_gap_max",
            "query_coverage_median", "query_coverage_min",
            "alternative_mapping_count_median", "alternative_mapping_count_max",
            "secondary_alignment_read_count",
            "ambiguous_read_count", "ambiguous_read_fraction",
            "read_names",
        ]
        tw = csv.DictWriter(tout, fieldnames=fields, delimiter="\t")
        tw.writeheader()
        bw = csv.writer(bout, delimiter="\t", lineterminator="\n")
        bw.writerow([
            "chrom1", "start1", "end1", "chrom2", "start2", "end2",
            "name", "score", "strand1", "strand2", "class"
        ])

        event_no = 0
        for cluster in clusters:
            read_names = sorted({m["read_name"] for m in cluster["members"]})
            support = len(read_names)
            if support < int(args.fusion_min_reads):
                continue
            event_no += 1
            cid = f"candidate{event_no}"
            # Preserve the established consensus breakpoint definition.
            p1 = _median_int([m["pos1"] for m in cluster["members"]])
            p2 = _median_int([m["pos2"] for m in cluster["members"]])
            chrom1 = cluster["representative"]["chrom1"]
            chrom2 = cluster["representative"]["chrom2"]
            jclass = cluster["representative"].get(
                "class", "interchromosomal" if chrom1 != chrom2 else "intrachromosomal"
            )
            qc = _candidate_qc_stats(cluster)
            tw.writerow({
                "candidate_id": cid,
                "support": support,
                "class": jclass,
                "chrom1": chrom1,
                "breakpoint1": p1,
                "chrom2": chrom2,
                "breakpoint2": p2,
                "junction_observations": qc["junction_observations"],
                "breakpoint1_mad": f'{qc["breakpoint1_mad"]:.3f}',
                "breakpoint2_mad": f'{qc["breakpoint2_mad"]:.3f}',
                "breakpoint1_span": qc["breakpoint1_span"],
                "breakpoint2_span": qc["breakpoint2_span"],
                "query_gap_median": qc["query_gap_median"],
                "query_gap_mad": f'{qc["query_gap_mad"]:.3f}',
                "query_gap_min": qc["query_gap_min"],
                "query_gap_max": qc["query_gap_max"],
                "query_coverage_median": f'{qc["query_coverage_median"]:.6f}',
                "query_coverage_min": f'{qc["query_coverage_min"]:.6f}',
                "alternative_mapping_count_median": f'{qc["alternative_mapping_count_median"]:.3f}',
                "alternative_mapping_count_max": qc["alternative_mapping_count_max"],
                "secondary_alignment_read_count": qc["secondary_alignment_read_count"],
                "ambiguous_read_count": qc["ambiguous_read_count"],
                "ambiguous_read_fraction": f'{qc["ambiguous_read_fraction"]:.6f}',
                "read_names": ",".join(read_names),
            })
            bw.writerow([
                chrom1, max(0, p1 - 1), p1,
                chrom2, max(0, p2 - 1), p2,
                cid, support, ".", ".", jclass,
            ])

    return len(observations), sum(
        1 for c in clusters
        if len({m["read_name"] for m in c["members"]}) >= int(args.fusion_min_reads)
    )


def annotate_fusions(args):
    # Standalone mode may target a new directory; create it automatically.
    out_dir = os.path.dirname(os.path.abspath(args.output_prefix))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    print("Annotating chimeric/fusion transcripts...")
    print(f"Fusion transcript reconstruction QC mode: {RECONSTRUCTION_QCFIX_VERSION} (legacy adapter/gap QC is non-blocking)")
    print(f"Fusion candidate QC mode: {CANDIDATE_QC_VERSION} (query coverage + mapping ambiguity; QC-only)")
    tss_index, tss_support = _load_end_index(args.gff_prefix + "_tss.gff3")
    tes_index, tes_support = _load_end_index(args.gff_prefix + "_tes.gff3")
    intron_index = _load_intron_index(args.gff_prefix + "_intron.gff3")

    groups = defaultdict(list)
    bam = pysam.AlignmentFile(args.in_bam, "rb")
    for read in bam:
        if read.is_unmapped:
            continue
        if int(_get_tag(read, "ZN", 1)) < 2:
            continue
        groups[read.query_name].append(_segment_from_read(read))
    bam.close()

    # LEVEL 1: every structurally reconstructed breakpoint event.
    structural_obs, structural_events = _write_structural_candidates(groups, args)

    # LEVEL 2: every split molecule is retained with orientation, adapter QC,
    # TSS/TES annotation/counts and validation reasons. Nothing is silently lost.
    molecule_rows = []
    validated_candidates = []
    annotations = {}
    for read_name, segments in groups.items():
        ann = _annotate_group(
            read_name, segments, tss_index, tes_index, tss_support, tes_support,
            intron_index, args
        )
        ann["reconstruction"] = _reconstruction_info(ann)
        annotations[read_name] = ann
        if ann["validated"]:
            validated_candidates.append(ann)

        seg_sorted = sorted(segments, key=lambda x: x["segment_index"])
        molecule_rows.append({
            "read_name": read_name,
            "segment_count": len(seg_sorted),
            "direction": ann["direction"],
            "orientation_source": ann["orientation_source"],
            "ts_consensus": ann["ts_consensus"],
            "adapter_direction": ann["adapter_direction"],
            "five_prime_adapter_valid": ann["five_prime_adapter_valid"],
            "three_prime_adapter_valid": ann["three_prime_adapter_valid"],
            "terminal_adapter_valid": ann["terminal_adapter_valid"],
            "left5": ann["left5"], "left3": ann["left3"],
            "right5": ann["right5"], "right3": ann["right3"],
            "internal_adapter": ann["internal_adapter"],
            "tss_contig": ann["tss_contig"],
            "tss_strand": ann["tss_strand"],
            "tss_raw_pos": ann["tss_raw"] if ann["tss_raw"] is not None else ".",
            "tss_pos": ann["tss"] if ann["tss"] is not None else ".",
            "tss_feature_count": ann["tss_feature_count"] if ann["tss_feature_count"] is not None else ".",
            "tes_contig": ann["tes_contig"],
            "tes_strand": ann["tes_strand"],
            "tes_raw_pos": ann["tes_raw"] if ann["tes_raw"] is not None else ".",
            "tes_pos": ann["tes"] if ann["tes"] is not None else ".",
            "tes_feature_count": ann["tes_feature_count"] if ann["tes_feature_count"] is not None else ".",
            "introns_valid": ann["introns_valid"],
            "full_sequence_recovered": ann["full_sequence_recovered"],
            "class": ann["class"],
            "validation_status": "validated" if ann["validated"] else "not_validated",
            "validation_reasons": "." if ann["validated"] else ";".join(ann["validation_reasons"]),
            "reconstruction_status": "reconstructable" if ann["reconstruction"]["reconstructable"] else "not_reconstructable",
            "reconstruction_reasons": "." if ann["reconstruction"]["reconstructable"] else ";".join(ann["reconstruction"]["reconstruction_reasons"]),
            "reconstruction_tss_pos": ann["reconstruction"]["tss_effective"] if ann["reconstruction"]["tss_effective"] is not None else ".",
            "reconstruction_tss_source": ann["reconstruction"]["tss_source"],
            "reconstruction_tes_pos": ann["reconstruction"]["tes_effective"] if ann["reconstruction"]["tes_effective"] is not None else ".",
            "reconstruction_tes_source": ann["reconstruction"]["tes_source"],
            "query_intervals": ";".join(f'{s["qstart"]}-{s["qend"]}' for s in seg_sorted),
            "segments": ";".join(
                f'{s["contig"]}:{s["ref_start"]}-{s["ref_end"]}:{s["original_strand"]}:ts={s.get("ts_tag", ".")}'
                for s in seg_sorted
            ),
        })

    # Sanity check: adapter-independent reconstruction must never be blocked by
    # legacy terminal-adapter/gap/internal-adapter QC flags.
    forbidden = {
        "terminal_adapter_validation_failed",
        "internal_adapter_detected",
        "large_unaligned_query_gap",
        "gap_or_chimera_artifact_tag",
    }
    leaked = []
    for rn, ann in annotations.items():
        rr = set(ann["reconstruction"]["reconstruction_reasons"])
        bad = rr & forbidden
        if bad:
            leaked.append((rn, sorted(bad)))
    if leaked:
        raise RuntimeError(
            f"QC-fix invariant failed for {len(leaked)} reads; example: {leaked[0]}"
        )

    molecule_path = args.output_prefix + "_fusion_molecules.tsv"
    molecule_fields = [
        "read_name", "segment_count", "direction", "orientation_source", "ts_consensus",
        "adapter_direction", "five_prime_adapter_valid", "three_prime_adapter_valid",
        "terminal_adapter_valid", "left5", "left3", "right5", "right3",
        "internal_adapter", "tss_contig", "tss_strand", "tss_raw_pos", "tss_pos",
        "tss_feature_count", "tes_contig", "tes_strand", "tes_raw_pos", "tes_pos",
        "tes_feature_count", "introns_valid", "full_sequence_recovered", "class",
        "validation_status", "validation_reasons", "reconstruction_status",
        "reconstruction_reasons", "reconstruction_tss_pos", "reconstruction_tss_source",
        "reconstruction_tes_pos", "reconstruction_tes_source", "query_intervals", "segments",
    ]
    with open(molecule_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=molecule_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(molecule_rows)

    # Backward-compatible read-level filename: same comprehensive no-discard table.
    read_path = args.output_prefix + "_fusion_reads.tsv"
    with open(read_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=molecule_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(molecule_rows)

    # TRANSCRIPT RECONSTRUCTION: LoRTIA accepted TSS/TES/introns define the
    # biological isoform; external adapters are optional positive evidence.
    reconstructed_reads, reconstructed_chains, adapter_chains = _write_reconstructed_transcripts(annotations, args)

    # LEVEL 3: validated full fusion transcript chains.
    clusters = _cluster_candidates(validated_candidates, int(args.fusion_wobble))
    fusion_path = args.output_prefix + "_fusions.tsv"
    bedpe_path = args.output_prefix + "_fusion_junctions.bedpe"

    with open(fusion_path, "w", newline="", encoding="utf-8") as fout, \
         open(bedpe_path, "w", newline="", encoding="utf-8") as bout:
        ffields = [
            "fusion_id", "support", "segment_count", "class", "direction",
            "orientation_source", "ts_consensus", "adapter_valid_support",
            "tss_contig", "tss_pos", "tss_feature_count", "tss_fusion_read_support",
            "tes_contig", "tes_pos", "tes_feature_count", "tes_fusion_read_support",
            "full_sequence_support", "read_names", "segments", "junctions",
        ]
        fw = csv.DictWriter(fout, fieldnames=ffields, delimiter="\t")
        fw.writeheader()
        bw = csv.writer(bout, delimiter="\t", lineterminator="\n")
        bw.writerow([
            "chrom1", "start1", "end1", "chrom2", "start2", "end2",
            "name", "score", "strand1", "strand2", "class"
        ])

        fusion_no = 0
        for cluster in clusters:
            members = cluster["members"]
            support = len({m["read_name"] for m in members})
            if support < int(args.fusion_min_reads):
                continue
            fusion_no += 1
            fusion_id = f"fusion{fusion_no}"
            segments = _cluster_segments(cluster)
            classes = [_junction_class(a, b, args.fusion_local_distance)
                       for a, b in zip(segments, segments[1:])]
            chain_class = ",".join(sorted(set(classes)))
            sources = defaultdict(int)
            ts_states = defaultdict(int)
            for m in members:
                sources[m["orientation_source"]] += 1
                ts_states[m["ts_consensus"]] += 1
            fw.writerow({
                "fusion_id": fusion_id,
                "support": support,
                "segment_count": len(segments),
                "class": chain_class,
                "direction": members[0]["direction"],
                "orientation_source": ",".join(f"{k}:{v}" for k, v in sorted(sources.items())),
                "ts_consensus": ",".join(f"{k}:{v}" for k, v in sorted(ts_states.items())),
                "adapter_valid_support": sum(m["terminal_adapter_valid"] for m in members),
                "tss_contig": segments[0]["contig"],
                "tss_pos": members[0]["tss"],
                "tss_feature_count": members[0]["tss_feature_count"],
                "tss_fusion_read_support": sum(1 for m in members if m["tss"] == members[0]["tss"]),
                "tes_contig": segments[-1]["contig"],
                "tes_pos": members[0]["tes"],
                "tes_feature_count": members[0]["tes_feature_count"],
                "tes_fusion_read_support": sum(1 for m in members if m["tes"] == members[0]["tes"]),
                "full_sequence_support": sum(m["full_sequence_recovered"] for m in members),
                "read_names": ",".join(sorted({m["read_name"] for m in members})),
                "segments": ";".join(_segment_text(s) for s in segments),
                "junctions": ";".join(_junction_text(a, b) for a, b in zip(segments, segments[1:])),
            })

            for j, (a, b) in enumerate(zip(segments, segments[1:]), start=1):
                p1, p2 = int(a["bio3"]), int(b["bio5"])
                jclass = _junction_class(a, b, args.fusion_local_distance)
                bw.writerow([
                    a["contig"], max(0, p1 - 1), p1,
                    b["contig"], max(0, p2 - 1), p2,
                    f"{fusion_id}:J{j}", support,
                    a["strand"], b["strand"], jclass,
                ])

    validated_chains = sum(
        1 for c in clusters
        if len({m["read_name"] for m in c["members"]}) >= int(args.fusion_min_reads)
    )
    oriented = sum(1 for a in annotations.values() if a["direction"] in {"5to3", "3to5"})
    ts_oriented = sum(1 for a in annotations.values() if a["orientation_source"] == "minimap2_ts")
    print(
        f"Fusion/chimeric output: {structural_obs} structural junction observations, "
        f"{structural_events} structural candidate events; {len(groups)} split molecules retained, "
        f"{oriented} oriented ({ts_oriented} from minimap2 ts consensus); "
        f"{reconstructed_reads} reconstructable reads, {reconstructed_chains} reconstructed fusion transcript isoforms "
        f"({adapter_chains} with >= {int(args.fusion_min_reads)} both-terminal-adapter reads); "
        f"{len(validated_candidates)} legacy strict LoRTIA-validated reads, {validated_chains} legacy strict transcript chains."
    )


def parsing():
    parser = ArgumentParser(description="Annotate LoRTIA split alignments as validated chimeric/fusion transcript chains.")
    parser.add_argument("in_bam", help="LoRTIA full processed BAM containing supplementary segments and Z* split tags.")
    parser.add_argument("gff_prefix", help="Prefix of accepted *_tss.gff3, *_tes.gff3 and *_intron.gff3 files.")
    parser.add_argument("output_prefix", help="Output prefix.")
    parser.add_argument("--wobble", type=int, default=10, help="TSS/TES matching window (+/- nt).")
    parser.add_argument("--fusion-wobble", dest="fusion_wobble", type=int, default=50,
                        help="Maximum breakpoint difference for clustering fusion reads. Default 50 nt.")
    parser.add_argument("--fusion-min-reads", dest="fusion_min_reads", type=int, default=1,
                        help="Minimum distinct reads supporting a reported chimeric transcript. Default 1.")
    parser.add_argument("--fusion-max-query-gap", dest="fusion_max_query_gap", type=int, default=20,
                        help="Maximum unaligned query gap between adjacent retained segments. Default 20 nt.")
    parser.add_argument("--fusion-local-distance", dest="fusion_local_distance", type=int, default=100000,
                        help="Same-contig distance below which a collinear split is labelled local. This is classification only, not a rejection threshold. Default 100000 nt.")
    parser.add_argument("--fusion-local-filter-distance", dest="fusion_local_filter_distance", type=int, default=1000,
                        help="Maximum genomic distance for classifying a same-contig same-strand collinear structural split as local. Default 1000 nt.")
    parser.set_defaults(fusion_exclude_local_collinear=True)
    parser.add_argument("--fusion-exclude-local-collinear", dest="fusion_exclude_local_collinear",
                        action="store_true",
                        help="Exclude local same-contig same-strand collinear split events from structural fusion-candidate output (default: enabled).")
    parser.add_argument("--fusion-keep-local-collinear", dest="fusion_exclude_local_collinear",
                        action="store_false",
                        help="Retain local same-contig same-strand collinear split events in structural candidate output.")
    parser.add_argument("--fusion-allow-internal-adapters", dest="fusion_allow_internal_adapters",
                        action="store_true", default=False,
                        help="Allow validated chimeric calls even when a correct adapter is detected at an internal split boundary. Disabled by default because this pattern is consistent with concatenated library molecules.")
    return parser.parse_args()


def main():
    annotate_fusions(parsing())


if __name__ == "__main__":
    main()
