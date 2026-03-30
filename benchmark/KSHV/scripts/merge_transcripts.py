#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Merge transcript annotations from multiple GFF/GFF3/GTF files by exon structure.

The script:
- scans an input directory recursively for GFF/GFF3/GTF files
- loads transcript and exon annotations from each file
- merges transcripts by exon composition
- optionally applies a terminal window for fuzzy transcript merging
- writes merged transcript tables and GFF3 outputs
- can also collapse identical TSS/TES/intron sites across transcripts

Outputs:
- <outdir_name>.tsv
- <outdir_name>.gff3
- <outdir_name>.TSS.gff3
- <outdir_name>.TES.gff3
- <outdir_name>.introns.gff3
"""

from __future__ import annotations

import argparse
import os
import sys
from collections import Counter, defaultdict
from typing import Any, Dict, Iterable, List, Optional, Tuple


def parse_gff_attributes(attr_field: str) -> Dict[str, str]:
    """
    Parse attributes from GFF3 or GTF format into a dictionary.

    GFF3 example:
        key=value;key2=value2

    GTF example:
        key "value"; key2 "value2";
    """
    result: Dict[str, str] = {}
    if not attr_field or attr_field == ".":
        return result

    text = attr_field.strip()

    # GFF3-style attributes
    if "=" in text:
        parts = [p.strip() for p in text.split(";") if p.strip()]
        for part in parts:
            if "=" in part:
                key, value = part.split("=", 1)
                result[key.strip()] = value.strip()
            else:
                tokens = part.split()
                if len(tokens) >= 2:
                    result[tokens[0].strip()] = " ".join(tokens[1:]).strip().strip('"')
        return result

    # GTF-style attributes
    parts = [p.strip() for p in text.split(";") if p.strip()]
    for part in parts:
        tokens = part.split(None, 1)
        if not tokens:
            continue

        key = tokens[0].strip()
        if not key:
            continue

        value = ""
        if len(tokens) > 1:
            value = tokens[1].strip()
            if value.startswith('"') and value.endswith('"') and len(value) >= 2:
                value = value[1:-1]
            else:
                value = value.strip().strip('"')

        result[key] = value

    return result


def iter_gff_records(path: str) -> Iterable[Dict[str, object]]:
    """
    Yield parsed feature records from a GFF/GFF3/GTF file.
    Supports plain text and .gz files.
    """
    opener = open
    if path.endswith(".gz"):
        import gzip
        opener = lambda p, mode: gzip.open(p, mode)  # type: ignore

    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue

            seqid, source, feature_type, start, end, score, strand, phase, attrs = cols[:9]

            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                continue

            if score in (".", ""):
                score_val = None
            else:
                try:
                    score_val = float(score)
                except ValueError:
                    score_val = None

            yield {
                "seqid": seqid,
                "source": source,
                "type": feature_type,
                "start": start_i,
                "end": end_i,
                "score": score_val,
                "strand": strand,
                "phase": phase,
                "attrs": parse_gff_attributes(attrs),
            }


class SampleData:
    """
    Container for transcript and exon information from a single input file.
    """

    def __init__(self, name: str):
        self.name = name
        self.transcript_exons: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
        self.transcript_loci: Dict[str, Tuple[str, str]] = {}

        self.mrna_score: Dict[str, float] = {}
        self.exon_score_sum: Dict[str, float] = defaultdict(float)
        self.exon_score_max: Dict[str, float] = defaultdict(lambda: float("-inf"))
        self.exon_score_n: Dict[str, int] = defaultdict(int)
        self.exon_score_first: Dict[str, float] = {}

    @staticmethod
    def get_transcript_id(attrs: Dict[str, str]) -> Optional[str]:
        return attrs.get("ID") or attrs.get("transcript_id") or attrs.get("Name")

    def add_transcript_record(self, record: Dict[str, object]) -> None:
        attrs: Dict[str, str] = record["attrs"]  # type: ignore
        transcript_id = self.get_transcript_id(attrs)
        if not transcript_id:
            return

        self.transcript_loci[transcript_id] = (record["seqid"], record["strand"])  # type: ignore

        score: Optional[float] = record["score"]  # type: ignore
        if score is not None:
            self.mrna_score[transcript_id] = float(score)

    def add_exon_record(self, record: Dict[str, object]) -> None:
        attrs: Dict[str, str] = record["attrs"]  # type: ignore

        parents: List[str] = []
        if "Parent" in attrs:
            parents = [p.strip() for p in attrs["Parent"].split(",") if p.strip()]
        elif "transcript_id" in attrs:
            parents = [attrs["transcript_id"]]

        if not parents:
            return

        start = int(record["start"])  # type: ignore
        end = int(record["end"])  # type: ignore
        seqid = record["seqid"]  # type: ignore
        strand = record["strand"]  # type: ignore

        score: Optional[float] = record["score"]  # type: ignore
        score_val = None if score is None else float(score)

        for transcript_id in parents:
            self.transcript_exons[transcript_id].append((start, end))

            if transcript_id not in self.transcript_loci:
                self.transcript_loci[transcript_id] = (seqid, strand)

            if score_val is not None:
                self.exon_score_sum[transcript_id] += score_val
                self.exon_score_n[transcript_id] += 1

                if transcript_id not in self.exon_score_first:
                    self.exon_score_first[transcript_id] = score_val

                if score_val > self.exon_score_max[transcript_id]:
                    self.exon_score_max[transcript_id] = score_val


def load_sample(path: str, sample_name: Optional[str] = None) -> SampleData:
    """
    Load one annotation file into a SampleData object.
    """
    if sample_name is None:
        base = os.path.basename(path)
        if base.endswith(".gz"):
            base = base[:-3]

        for ext in (".gff3", ".gff", ".gtf", ".GFF3", ".GFF", ".GTF"):
            if base.endswith(ext):
                base = base[:-len(ext)]
                break

        sample_name = base

    sample = SampleData(sample_name)

    for record in iter_gff_records(path):
        feature_type = record["type"]  # type: ignore
        if feature_type in ("mRNA", "transcript"):
            sample.add_transcript_record(record)
        elif feature_type == "exon":
            sample.add_exon_record(record)

    return sample


def normalized_exons(exons: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """
    Return a sorted, unique exon list.
    """
    return sorted(set((int(s), int(e)) for s, e in exons), key=lambda x: (x[0], x[1]))


def choose_count(sample: SampleData, transcript_id: str, exon_score_mode: str) -> float:
    """
    Choose a transcript count value.

    Priority:
    1. transcript/mRNA score
    2. exon-derived score using the requested mode
    3. default presence value = 1.0
    """
    if transcript_id in sample.mrna_score:
        return sample.mrna_score[transcript_id]

    n = sample.exon_score_n.get(transcript_id, 0)
    if n > 0:
        mode = exon_score_mode.lower()

        if mode == "sum":
            return float(sample.exon_score_sum[transcript_id])
        if mode == "mean":
            return float(sample.exon_score_sum[transcript_id]) / float(n)
        if mode == "first":
            return float(sample.exon_score_first[transcript_id])

        max_score = sample.exon_score_max[transcript_id]
        if max_score != float("-inf"):
            return float(max_score)

    return 1.0


def internal_signature(exons: List[Tuple[int, int]]) -> Tuple[int, ...]:
    """
    Build a signature for multi-exon transcripts that ignores only:
    - the first exon start
    - the last exon end

    This allows fuzzy merging where only transcript ends can move.
    """
    if len(exons) <= 1:
        return tuple()

    signature: List[int] = []
    signature.append(exons[0][1])  # end of first exon

    for start, end in exons[1:-1]:
        signature.extend([start, end])

    signature.append(exons[-1][0])  # start of last exon
    return tuple(signature)


def cluster_by_window(instances: List[Dict[str, Any]], window: int) -> List[List[Dict[str, Any]]]:
    """
    Greedy clustering using transcript start/end coordinates.

    Two transcripts are placed into the same cluster if both:
    - abs(start - ref_start) <= window
    - abs(end - ref_end) <= window

    The cluster reference is the first transcript in the cluster after sorting.
    """
    if window <= 0 or len(instances) <= 1:
        return [instances]

    sorted_instances = sorted(instances, key=lambda d: (d["start"], d["end"]))
    clusters: List[List[Dict[str, Any]]] = []

    for instance in sorted_instances:
        if not clusters:
            clusters.append([instance])
            continue

        ref = clusters[-1][0]
        if (
            abs(int(instance["start"]) - int(ref["start"])) <= window
            and abs(int(instance["end"]) - int(ref["end"])) <= window
        ):
            clusters[-1].append(instance)
        else:
            clusters.append([instance])

    return clusters


def merge_samples(
    sample_datas: List[SampleData],
    exon_score_mode: str = "max",
    window: int = 0,
):
    """
    Merge transcripts across samples by exon structure.

    Exact mode (window=0):
        transcripts must have identical exon composition

    Window mode (window>0):
        - single-exon transcripts: start and end may differ within the window
        - multi-exon transcripts: internal exon structure must match exactly,
          while the first exon start and last exon end may vary within the window
    """
    sample_names = [sample.name for sample in sample_datas]

    base_groups: Dict[Tuple[str, str, int, Tuple[int, ...]], List[Dict[str, Any]]] = defaultdict(list)

    for sample in sample_datas:
        for transcript_id, exons in sample.transcript_exons.items():
            if not exons:
                continue

            seqid, strand = sample.transcript_loci.get(transcript_id, (None, None))
            if not seqid or not strand:
                continue

            exons_norm = normalized_exons(exons)
            if not exons_norm:
                continue

            start = min(s for s, _ in exons_norm)
            end = max(e for _, e in exons_norm)
            exon_count = len(exons_norm)
            signature = internal_signature(exons_norm)

            count_value = float(choose_count(sample, transcript_id, exon_score_mode))

            base_key = (seqid, strand, exon_count, signature)
            base_groups[base_key].append(
                {
                    "seqid": seqid,
                    "strand": strand,
                    "exons": exons_norm,
                    "start": start,
                    "end": end,
                    "sample": sample.name,
                    "count": count_value,
                }
            )

    key_to_counts: Dict[Tuple[str, str, Tuple[Tuple[int, int], ...]], Counter] = defaultdict(Counter)
    key_to_locus: Dict[Tuple[str, str, Tuple[Tuple[int, int], ...]], Tuple[int, int]] = {}

    for (seqid, strand, exon_count, _signature), instances in base_groups.items():
        if window <= 0:
            for instance in instances:
                exons_norm = instance["exons"]
                key = (seqid, strand, tuple(exons_norm))
                key_to_locus.setdefault(key, (int(instance["start"]), int(instance["end"])))
                key_to_counts[key][instance["sample"]] += float(instance["count"])
            continue

        clusters = cluster_by_window(instances, window=window)

        for cluster in clusters:
            cluster_start = min(int(item["start"]) for item in cluster)
            cluster_end = max(int(item["end"]) for item in cluster)

            representative_exons = [tuple(map(int, exon)) for exon in cluster[0]["exons"]]

            if exon_count == 1:
                merged_exons = [(cluster_start, cluster_end)]
            else:
                merged_exons = list(representative_exons)
                merged_exons[0] = (cluster_start, merged_exons[0][1])
                merged_exons[-1] = (merged_exons[-1][0], cluster_end)

            merged_exons = normalized_exons(merged_exons)
            key = (seqid, strand, tuple(merged_exons))
            key_to_locus[key] = (cluster_start, cluster_end)

            for instance in cluster:
                key_to_counts[key][instance["sample"]] += float(instance["count"])

    sorted_keys = sorted(
        key_to_counts.keys(),
        key=lambda key: (key[0], key[1], key_to_locus[key][0], key_to_locus[key][1], len(key[2])),
    )

    transcript_meta: Dict[Tuple[str, str, Tuple[Tuple[int, int], ...]], Dict[str, object]] = {}
    counts: Dict[str, Dict[str, float]] = {}

    for idx, key in enumerate(sorted_keys, start=1):
        seqid, strand, exons_norm = key
        start, end = key_to_locus[key]
        exon_comp = ";".join(f"{s}-{e}" for s, e in exons_norm)
        transcript_id = f"TR_{idx:06d}"

        transcript_meta[key] = {
            "tr_id": transcript_id,
            "seqid": seqid,
            "strand": strand,
            "start": start,
            "end": end,
            "exon_comp": exon_comp,
            "exon_count": len(exons_norm),
        }

        counts[transcript_id] = {
            sample_name: float(key_to_counts[key].get(sample_name, 0.0))
            for sample_name in sample_names
        }

    return transcript_meta, counts, sample_names


def int_if_integer(value: float):
    """
    Convert a float to int if it has no fractional part.
    """
    return int(value) if float(value).is_integer() else float(value)


def format_plain_number(value: float, decimals: int = 6) -> str:
    """
    Format numeric values without unnecessary trailing zeros.
    """
    if float(value).is_integer():
        return str(int(value))

    text = f"{value:.{decimals}f}".rstrip("0").rstrip(".")
    return text if text else "0"


def write_tsv(
    path: str,
    transcript_meta: Dict[Tuple[str, str, Tuple[Tuple[int, int], ...]], Dict[str, object]],
    counts: Dict[str, Dict[str, float]],
    sample_names: List[str],
) -> None:
    """
    Write merged transcript table as TSV.
    """
    import csv

    headers = ["isoform", "chromosome", "start", "end", "strand", "exon_composition"] + sample_names

    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(headers)

        for meta in sorted(transcript_meta.values(), key=lambda d: d["tr_id"]):  # type: ignore
            transcript_id = meta["tr_id"]  # type: ignore
            row = [
                transcript_id,
                meta["seqid"],
                meta["start"],
                meta["end"],
                meta["strand"],
                meta["exon_comp"],
            ] + [int_if_integer(counts[transcript_id].get(sample_name, 0.0)) for sample_name in sample_names]
            writer.writerow(row)


def write_gff3(
    path: str,
    transcript_meta: Dict[Tuple[str, str, Tuple[Tuple[int, int], ...]], Dict[str, object]],
    counts: Dict[str, Dict[str, float]],
    sample_names: List[str],
    write_sample_attrs: bool = False,
) -> None:
    """
    Write merged transcripts as a GFF3 file.
    """
    with open(path, "w", encoding="utf-8", newline="") as handle:
        handle.write("##gff-version 3\n")

        metas = sorted(transcript_meta.values(), key=lambda d: d["tr_id"])  # type: ignore
        for meta in metas:
            transcript_id: str = meta["tr_id"]  # type: ignore
            seqid: str = meta["seqid"]  # type: ignore
            strand: str = meta["strand"]  # type: ignore
            start: int = int(meta["start"])  # type: ignore
            end: int = int(meta["end"])  # type: ignore
            exon_comp: str = meta["exon_comp"]  # type: ignore

            exons = [tuple(map(int, part.split("-"))) for part in exon_comp.split(";") if part]
            total_count = sum(float(v) for v in counts.get(transcript_id, {}).values())
            score_str = "." if total_count == 0 else format_plain_number(total_count)

            attrs = [
                f"ID={transcript_id}",
                f"Name={transcript_id}",
                f"exon_count={len(exons)}",
                f"exon_comp={exon_comp}",
            ]

            if write_sample_attrs:
                sample_attr = ";".join(
                    f"{sample_name}={format_plain_number(value)}"
                    for sample_name, value in counts.get(transcript_id, {}).items()
                )
                if sample_attr:
                    attrs.append(sample_attr)

            attr_str = ";".join(attrs)

            handle.write(
                "\t".join(
                    [seqid, "merge", "mRNA", str(start), str(end), score_str, strand, ".", attr_str]
                )
                + "\n"
            )

            for exon_idx, (exon_start, exon_end) in enumerate(exons, start=1):
                handle.write(
                    "\t".join(
                        [
                            seqid,
                            "merge",
                            "exon",
                            str(exon_start),
                            str(exon_end),
                            ".",
                            strand,
                            ".",
                            f"ID={transcript_id}.exon{exon_idx};Parent={transcript_id}",
                        ]
                    )
                    + "\n"
                )


def write_tss_gff3(
    path: str,
    transcript_meta: Dict[Tuple[str, str, Tuple[Tuple[int, int], ...]], Dict[str, object]],
    counts: Dict[str, Dict[str, float]],
    collapse_sites: bool = False,
) -> None:
    """
    Write TSS sites as GFF3.

    If collapse_sites is False:
        one TSS line per transcript

    If collapse_sites is True:
        identical (seqid, strand, position) sites are merged and their scores are summed
    """
    with open(path, "w", encoding="utf-8", newline="") as handle:
        handle.write("##gff-version 3\n")

        if not collapse_sites:
            metas = sorted(transcript_meta.values(), key=lambda d: d["tr_id"])  # type: ignore
            for meta in metas:
                transcript_id: str = meta["tr_id"]  # type: ignore
                seqid: str = meta["seqid"]  # type: ignore
                strand: str = meta["strand"]  # type: ignore
                start: int = int(meta["start"])  # type: ignore
                end: int = int(meta["end"])  # type: ignore

                tss_pos = start if strand == "+" else end
                total_count = sum(float(v) for v in counts.get(transcript_id, {}).values())
                score_str = "." if total_count == 0 else format_plain_number(total_count)

                attrs = f"ID={transcript_id}.TSS;Parent={transcript_id};Name={transcript_id}.TSS"
                handle.write(
                    "\t".join(
                        [seqid, "merge", "TSS", str(tss_pos), str(tss_pos), score_str, strand, ".", attrs]
                    )
                    + "\n"
                )
            return

        site_sum: Dict[Tuple[str, str, int], float] = defaultdict(float)

        for meta in transcript_meta.values():
            transcript_id: str = meta["tr_id"]  # type: ignore
            seqid: str = meta["seqid"]  # type: ignore
            strand: str = meta["strand"]  # type: ignore
            start: int = int(meta["start"])  # type: ignore
            end: int = int(meta["end"])  # type: ignore

            pos = start if strand == "+" else end
            total_count = sum(float(v) for v in counts.get(transcript_id, {}).values())
            site_sum[(seqid, strand, pos)] += float(total_count)

        items = sorted(site_sum.items(), key=lambda x: (x[0][0], x[0][1], x[0][2]))
        for idx, ((seqid, strand, pos), total_count) in enumerate(items, start=1):
            score_str = "." if total_count == 0 else format_plain_number(total_count)
            attrs = f"ID=TSS_{idx:06d};Name=TSS_{idx:06d}"
            handle.write(
                "\t".join([seqid, "merge", "TSS", str(pos), str(pos), score_str, strand, ".", attrs]) + "\n"
            )


def write_tes_gff3(
    path: str,
    transcript_meta: Dict[Tuple[str, str, Tuple[Tuple[int, int], ...]], Dict[str, object]],
    counts: Dict[str, Dict[str, float]],
    collapse_sites: bool = False,
) -> None:
    """
    Write TES sites as GFF3.

    If collapse_sites is False:
        one TES line per transcript

    If collapse_sites is True:
        identical (seqid, strand, position) sites are merged and their scores are summed
    """
    with open(path, "w", encoding="utf-8", newline="") as handle:
        handle.write("##gff-version 3\n")

        if not collapse_sites:
            metas = sorted(transcript_meta.values(), key=lambda d: d["tr_id"])  # type: ignore
            for meta in metas:
                transcript_id: str = meta["tr_id"]  # type: ignore
                seqid: str = meta["seqid"]  # type: ignore
                strand: str = meta["strand"]  # type: ignore
                start: int = int(meta["start"])  # type: ignore
                end: int = int(meta["end"])  # type: ignore

                tes_pos = end if strand == "+" else start
                total_count = sum(float(v) for v in counts.get(transcript_id, {}).values())
                score_str = "." if total_count == 0 else format_plain_number(total_count)

                attrs = f"ID={transcript_id}.TES;Parent={transcript_id};Name={transcript_id}.TES"
                handle.write(
                    "\t".join(
                        [seqid, "merge", "TES", str(tes_pos), str(tes_pos), score_str, strand, ".", attrs]
                    )
                    + "\n"
                )
            return

        site_sum: Dict[Tuple[str, str, int], float] = defaultdict(float)

        for meta in transcript_meta.values():
            transcript_id: str = meta["tr_id"]  # type: ignore
            seqid: str = meta["seqid"]  # type: ignore
            strand: str = meta["strand"]  # type: ignore
            start: int = int(meta["start"])  # type: ignore
            end: int = int(meta["end"])  # type: ignore

            pos = end if strand == "+" else start
            total_count = sum(float(v) for v in counts.get(transcript_id, {}).values())
            site_sum[(seqid, strand, pos)] += float(total_count)

        items = sorted(site_sum.items(), key=lambda x: (x[0][0], x[0][1], x[0][2]))
        for idx, ((seqid, strand, pos), total_count) in enumerate(items, start=1):
            score_str = "." if total_count == 0 else format_plain_number(total_count)
            attrs = f"ID=TES_{idx:06d};Name=TES_{idx:06d}"
            handle.write(
                "\t".join([seqid, "merge", "TES", str(pos), str(pos), score_str, strand, ".", attrs]) + "\n"
            )


def write_introns_gff3(
    path: str,
    transcript_meta: Dict[Tuple[str, str, Tuple[Tuple[int, int], ...]], Dict[str, object]],
    counts: Dict[str, Dict[str, float]],
    collapse_sites: bool = False,
) -> None:
    """
    Write introns as GFF3.

    If collapse_sites is False:
        write introns for each transcript separately

    If collapse_sites is True:
        identical intron intervals are merged and their scores are summed
    """
    with open(path, "w", encoding="utf-8", newline="") as handle:
        handle.write("##gff-version 3\n")

        if not collapse_sites:
            metas = sorted(transcript_meta.values(), key=lambda d: d["tr_id"])  # type: ignore
            for meta in metas:
                transcript_id: str = meta["tr_id"]  # type: ignore
                seqid: str = meta["seqid"]  # type: ignore
                strand: str = meta["strand"]  # type: ignore
                exon_comp: str = meta["exon_comp"]  # type: ignore

                exons = [tuple(map(int, part.split("-"))) for part in exon_comp.split(";") if part]
                exons = sorted(exons, key=lambda x: (x[0], x[1]))

                if len(exons) < 2:
                    continue

                intron_idx = 0
                for (s1, e1), (s2, _e2) in zip(exons, exons[1:]):
                    intron_start = e1 + 1
                    intron_end = s2 - 1

                    if intron_start <= intron_end:
                        intron_idx += 1
                        attrs = (
                            f"ID={transcript_id}.intron{intron_idx};"
                            f"Parent={transcript_id};"
                            f"Name={transcript_id}.intron{intron_idx}"
                        )
                        handle.write(
                            "\t".join(
                                [
                                    seqid,
                                    "merge",
                                    "intron",
                                    str(intron_start),
                                    str(intron_end),
                                    ".",
                                    strand,
                                    ".",
                                    attrs,
                                ]
                            )
                            + "\n"
                        )
            return

        intron_sum: Dict[Tuple[str, str, int, int], float] = defaultdict(float)

        for meta in transcript_meta.values():
            transcript_id: str = meta["tr_id"]  # type: ignore
            seqid: str = meta["seqid"]  # type: ignore
            strand: str = meta["strand"]  # type: ignore
            exon_comp: str = meta["exon_comp"]  # type: ignore

            transcript_total = sum(float(v) for v in counts.get(transcript_id, {}).values())

            exons = [tuple(map(int, part.split("-"))) for part in exon_comp.split(";") if part]
            exons = sorted(exons, key=lambda x: (x[0], x[1]))

            if len(exons) < 2:
                continue

            for (s1, e1), (s2, _e2) in zip(exons, exons[1:]):
                intron_start = e1 + 1
                intron_end = s2 - 1

                if intron_start <= intron_end:
                    intron_sum[(seqid, strand, intron_start, intron_end)] += float(transcript_total)

        items = sorted(intron_sum.items(), key=lambda x: (x[0][0], x[0][1], x[0][2], x[0][3]))
        for idx, ((seqid, strand, start, end), total_count) in enumerate(items, start=1):
            score_str = "." if total_count == 0 else format_plain_number(total_count)
            attrs = f"ID=INTRON_{idx:06d};Name=INTRON_{idx:06d}"
            handle.write(
                "\t".join(
                    [seqid, "merge", "intron", str(start), str(end), score_str, strand, ".", attrs]
                )
                + "\n"
            )


def collect_gff_files(indir: str, pattern: str = "auto") -> List[str]:
    """
    Recursively collect input annotation files from a directory.
    """
    patt = pattern.lower()

    if patt == "auto":
        extensions = (".gff3", ".gff", ".gtf", ".gff3.gz", ".gff.gz", ".gtf.gz")
    elif patt in (".gff3", "gff3"):
        extensions = (".gff3", ".gff3.gz")
    elif patt in (".gff", "gff"):
        extensions = (".gff", ".gff.gz")
    elif patt in (".gtf", "gtf"):
        extensions = (".gtf", ".gtf.gz")
    else:
        gz_pattern = pattern + ".gz" if not pattern.endswith(".gz") else pattern
        extensions = (pattern, gz_pattern)

    files: List[str] = []

    for root, _, filenames in os.walk(indir):
        for filename in filenames:
            lower_name = filename.lower()
            if any(lower_name.endswith(ext) for ext in extensions):
                path = os.path.join(root, filename)
                try:
                    if os.path.getsize(path) > 0:
                        files.append(path)
                except OSError:
                    continue

    files.sort()
    return files


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Merge multiple GFF/GFF3/GTF files by transcript exon structure."
    )
    parser.add_argument(
        "-i",
        "--indir",
        required=True,
        help="Input directory containing GFF/GFF3/GTF files (searched recursively).",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        required=True,
        help="Output directory for TSV and GFF3 files.",
    )
    parser.add_argument(
        "--pattern",
        default="auto",
        help="File extension filter: auto | .gff3 | .gff | .gtf (default: auto).",
    )
    parser.add_argument(
        "--write-sample-attrs",
        action="store_true",
        help="Include per-sample counts in merged mRNA attributes.",
    )
    parser.add_argument(
        "--exon-score-mode",
        default="max",
        choices=["max", "mean", "sum", "first"],
        help="How to derive transcript count from exon scores if transcript score is missing (default: max).",
    )
    parser.add_argument(
        "--window",
        type=int,
        default=0,
        help=(
            "Window size in bp for fuzzy transcript merging. "
            "For multi-exon transcripts only transcript ends may vary; "
            "internal exon structure must remain identical."
        ),
    )
    parser.add_argument(
        "--collapse-sites",
        action="store_true",
        help="Collapse identical TSS/TES/intron coordinates and sum their scores.",
    )

    args = parser.parse_args()

    indir = os.path.abspath(args.indir)
    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)

    base_name = os.path.basename(outdir.rstrip("/"))
    tsv_path = os.path.join(outdir, f"{base_name}.tsv")
    gff_path = os.path.join(outdir, f"{base_name}.gff3")
    tss_path = os.path.join(outdir, f"{base_name}.TSS.gff3")
    tes_path = os.path.join(outdir, f"{base_name}.TES.gff3")
    intron_path = os.path.join(outdir, f"{base_name}.introns.gff3")

    files = collect_gff_files(indir, args.pattern)
    if not files:
        sys.exit("No readable GFF/GFF3/GTF files were found in the input directory.")

    print(f"Loaded files: {len(files)}")
    samples = [load_sample(path) for path in files]

    transcript_meta, counts, sample_names = merge_samples(
        samples,
        exon_score_mode=args.exon_score_mode,
        window=max(0, int(args.window)),
    )

    write_tsv(tsv_path, transcript_meta, counts, sample_names)
    write_gff3(gff_path, transcript_meta, counts, sample_names, write_sample_attrs=args.write_sample_attrs)
    write_tss_gff3(tss_path, transcript_meta, counts, collapse_sites=args.collapse_sites)
    write_tes_gff3(tes_path, transcript_meta, counts, collapse_sites=args.collapse_sites)
    write_introns_gff3(intron_path, transcript_meta, counts, collapse_sites=args.collapse_sites)

    print(f"Done: {tsv_path}")
    print(f"Done: {gff_path}")
    print(f"Done: {tss_path}")
    print(f"Done: {tes_path}")
    print(f"Done: {intron_path}")
    print(f"Merged transcripts: {len(transcript_meta)}")
    print(f"Samples: {', '.join(sample_names)}")

    if args.window > 0:
        print(f"Window-based merging enabled: {args.window} bp")

    if args.collapse_sites:
        print("Site collapsing enabled for TSS/TES/introns")


if __name__ == "__main__":
    main()