#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
merge_tsv_transcripts_stats_by_windows_clean.py

Merge transcript rows from multiple TSV files using separate TSS, TES, and intron windows.

Main features:
- generates new shared Cluster_ID values
- supports multiple input TSV files
- prefixes *_present columns using the input file name
- can disable merging within the same input file
- can apply optional 1:1 greedy matching between input files
- keeps only the requested minimal output columns
- missing *_present values are written as 0
- outputs one row per cluster
"""

from __future__ import annotations

import argparse
import csv
import os
import re
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional


def natural_basename_noext(path: str) -> str:
    """Return a sanitized file stem without extension."""
    base = os.path.basename(path)
    stem, _ = os.path.splitext(base)
    stem = re.sub(r"[^\w.\-]+", "_", stem)
    return stem


def detect_column(columns: List[str], candidates: List[str], required: bool = True) -> Optional[str]:
    """Find a matching column name case-insensitively."""
    lower_map = {c.lower(): c for c in columns}
    for candidate in candidates:
        if candidate.lower() in lower_map:
            return lower_map[candidate.lower()]
    if required:
        raise ValueError(f"Missing required column. None of these were found: {candidates}")
    return None


def parse_int(value: str, default: Optional[int] = None) -> Optional[int]:
    """Parse integer values while tolerating commas and spaces."""
    if value is None:
        return default
    text = str(value).strip().replace(",", "").replace("\u00A0", "").replace(" ", "")
    if text == "":
        return default
    try:
        return int(text)
    except ValueError:
        return default


def parse_exon_comp(exon_comp: str) -> List[Tuple[int, int]]:
    """Parse exon composition string into sorted exon tuples."""
    if exon_comp is None:
        return []

    exon_comp = exon_comp.strip()
    if not exon_comp:
        return []

    exons = []
    for part in exon_comp.split(";"):
        part = part.strip()
        if not part:
            continue

        match = re.match(r"^\s*(\d+)\s*-\s*(\d+)\s*$", part)
        if not match:
            continue

        a, b = int(match.group(1)), int(match.group(2))
        if a <= b:
            exons.append((a, b))
        else:
            exons.append((b, a))

    exons.sort()
    return exons


def introns_from_exons(exons: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Compute introns from exon coordinates."""
    if len(exons) < 2:
        return []

    introns = []
    for i in range(len(exons) - 1):
        intron_start = exons[i][1] + 1
        intron_end = exons[i + 1][0] - 1
        if intron_start <= intron_end:
            introns.append((intron_start, intron_end))
    return introns


def introns_match(
    introns_a: List[Tuple[int, int]],
    introns_b: List[Tuple[int, int]],
    window: int,
) -> bool:
    """Check if two intron chains match within the given window."""
    if len(introns_a) != len(introns_b):
        return False

    for (a1, a2), (b1, b2) in zip(introns_a, introns_b):
        if abs(a1 - b1) > window or abs(a2 - b2) > window:
            return False

    return True


def tss_tes_from_row(start: int, end: int, strand: str) -> Tuple[int, int]:
    """Return TSS and TES based on strand."""
    strand = (strand or "").strip()
    if strand == "-":
        return end, start
    return start, end


@dataclass
class TranscriptRow:
    idx: int
    source_prefix: str
    original_id: str
    chrom: str
    start: int
    end: int
    strand: str
    exon_comp: str
    exons: List[Tuple[int, int]]
    introns: List[Tuple[int, int]]
    tss: int
    tes: int
    raw_present: Dict[str, str]


class UnionFind:
    """Disjoint-set data structure for clustering matches."""

    def __init__(self, n: int):
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, x: int) -> int:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: int, b: int) -> None:
        root_a = self.find(a)
        root_b = self.find(b)

        if root_a == root_b:
            return

        if self.rank[root_a] < self.rank[root_b]:
            self.parent[root_a] = root_b
        elif self.rank[root_a] > self.rank[root_b]:
            self.parent[root_b] = root_a
        else:
            self.parent[root_b] = root_a
            self.rank[root_a] += 1


def compute_match_score(a: TranscriptRow, b: TranscriptRow) -> Tuple[int, int, int]:
    """Return a score tuple used for greedy pairwise matching."""
    tss_diff = abs(a.tss - b.tss)
    tes_diff = abs(a.tes - b.tes)

    intron_sum = 0
    if len(a.introns) == len(b.introns):
        for (a1, a2), (b1, b2) in zip(a.introns, b.introns):
            intron_sum += abs(a1 - b1) + abs(a2 - b2)
    else:
        intron_sum = 10**9

    return (tss_diff, tes_diff, intron_sum)


def rows_match(
    a: TranscriptRow,
    b: TranscriptRow,
    tss_window: int,
    tes_window: int,
    intron_window: int,
) -> bool:
    """Check if two transcript rows match within the requested windows."""
    if a.chrom != b.chrom:
        return False
    if a.strand != b.strand:
        return False
    if abs(a.tss - b.tss) > tss_window:
        return False
    if abs(a.tes - b.tes) > tes_window:
        return False
    if not introns_match(a.introns, b.introns, intron_window):
        return False
    return True


def build_pairwise_edges(
    rows: List[TranscriptRow],
    tss_window: int,
    tes_window: int,
    intron_window: int,
    no_within_source_merge: bool,
    pairwise_exclusive: bool,
) -> List[Tuple[int, int]]:
    """Build matching edges between transcript rows."""
    by_source: Dict[str, List[TranscriptRow]] = {}
    for row in rows:
        by_source.setdefault(row.source_prefix, []).append(row)

    sources = sorted(by_source.keys())
    edges: List[Tuple[int, int]] = []

    if not no_within_source_merge:
        for source in sources:
            group = by_source[source]
            n = len(group)
            for i in range(n):
                for j in range(i + 1, n):
                    a = group[i]
                    b = group[j]
                    if rows_match(a, b, tss_window, tes_window, intron_window):
                        edges.append((a.idx, b.idx))

    for i in range(len(sources)):
        for j in range(i + 1, len(sources)):
            src_a = sources[i]
            src_b = sources[j]
            group_a = by_source[src_a]
            group_b = by_source[src_b]

            candidates: List[Tuple[Tuple[int, int, int], int, int]] = []
            for a in group_a:
                for b in group_b:
                    if rows_match(a, b, tss_window, tes_window, intron_window):
                        score = compute_match_score(a, b)
                        candidates.append((score, a.idx, b.idx))

            candidates.sort(key=lambda x: x[0])

            if pairwise_exclusive:
                used_a = set()
                used_b = set()
                for _, aidx, bidx in candidates:
                    if aidx in used_a or bidx in used_b:
                        continue
                    edges.append((aidx, bidx))
                    used_a.add(aidx)
                    used_b.add(bidx)
            else:
                for _, aidx, bidx in candidates:
                    edges.append((aidx, bidx))

    return edges


def load_tsv(path: str, start_idx: int = 0) -> List[TranscriptRow]:
    """Load one TSV file into TranscriptRow records."""
    rows: List[TranscriptRow] = []
    source_prefix = natural_basename_noext(path)

    with open(path, "r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"No header found in file: {path}")

        columns = reader.fieldnames

        id_col = detect_column(columns, ["Transcript_ID", "transcript_id", "ID", "Id"])
        chrom_col = detect_column(columns, ["chromosome", "Chromosome", "chrom", "seqid"])
        start_col = detect_column(columns, ["Start", "start"])
        end_col = detect_column(columns, ["End", "end", "Stop", "stop"])
        strand_col = detect_column(columns, ["Strand", "strand"])
        exon_col = detect_column(columns, ["Exon_Composition", "exon_composition", "Exons"], required=False)

        idx = start_idx
        for record in reader:
            start = parse_int(record.get(start_col, ""))
            end = parse_int(record.get(end_col, ""))
            if start is None or end is None:
                continue

            chrom = (record.get(chrom_col) or "").strip()
            strand = (record.get(strand_col) or "").strip()
            original_id = (record.get(id_col) or "").strip()

            exon_comp = (record.get(exon_col) or "").strip() if exon_col else ""
            if start > end:
                start, end = end, start

            if not exon_comp:
                exon_comp = f"{start}-{end}"

            exons = parse_exon_comp(exon_comp)
            if not exons:
                exons = [(start, end)]
                exon_comp = f"{start}-{end}"

            introns = introns_from_exons(exons)
            tss, tes = tss_tes_from_row(start, end, strand)

            raw_present: Dict[str, str] = {}
            for key, value in record.items():
                if key.endswith("_present"):
                    val = str(value).strip()
                    raw_present[f"{source_prefix}__{key}"] = "0" if val == "" else val

            rows.append(
                TranscriptRow(
                    idx=idx,
                    source_prefix=source_prefix,
                    original_id=original_id,
                    chrom=chrom,
                    start=start,
                    end=end,
                    strand=strand,
                    exon_comp=exon_comp,
                    exons=exons,
                    introns=introns,
                    tss=tss,
                    tes=tes,
                    raw_present=raw_present,
                )
            )
            idx += 1

    return rows


def build_cluster_table(rows: List[TranscriptRow], uf: UnionFind) -> List[Dict[str, str]]:
    """Build the final cluster-level output table."""
    root_to_members: Dict[int, List[TranscriptRow]] = {}
    for row in rows:
        root = uf.find(row.idx)
        root_to_members.setdefault(root, []).append(row)

    cluster_roots = sorted(
        root_to_members.keys(),
        key=lambda root: min(member.idx for member in root_to_members[root]),
    )

    root_to_cluster_id = {
        root: f"CLUST_{i:06d}"
        for i, root in enumerate(cluster_roots, start=1)
    }

    present_cols = set()
    for row in rows:
        for col in row.raw_present.keys():
            if col.endswith("_present"):
                present_cols.add(col)
    present_cols = sorted(present_cols)

    out_rows: List[Dict[str, str]] = []

    for root in cluster_roots:
        members = root_to_members[root]
        cluster_id = root_to_cluster_id[root]

        cluster_start = min(row.start for row in members)
        cluster_end = max(row.end for row in members)

        cluster_chroms = sorted(set(row.chrom for row in members if row.chrom))
        cluster_chrom = cluster_chroms[0] if len(cluster_chroms) == 1 else ";".join(cluster_chroms)

        cluster_strands = sorted(set(row.strand for row in members if row.strand))
        cluster_strand = cluster_strands[0] if len(cluster_strands) == 1 else ";".join(cluster_strands)

        representative = members[0]
        cluster_exon_comp = representative.exon_comp
        cluster_introns = ";".join(f"{a}-{b}" for a, b in representative.introns)

        merged_present = {col: "0" for col in present_cols}
        for row in members:
            for col in present_cols:
                value = str(row.raw_present.get(col, "")).strip()
                if value == "1":
                    merged_present[col] = "1"
                elif value in {"", "0"}:
                    pass
                else:
                    if merged_present[col] == "0":
                        merged_present[col] = value

        merged_ids = ";".join(sorted(set(row.original_id for row in members if row.original_id)))

        record: Dict[str, str] = {
            "Cluster_ID": cluster_id,
            "Original_Transcript_ID": merged_ids,
            "chrom": cluster_chrom,
            "Start": str(cluster_start),
            "End": str(cluster_end),
            "Strand": cluster_strand,
            "Exon_Composition": cluster_exon_comp,
            "Introns": cluster_introns,
        }

        for col in present_cols:
            record[col] = merged_present.get(col, "0") or "0"

        out_rows.append(record)

    return out_rows


def write_tsv(path: str, rows: List[Dict[str, str]]) -> None:
    """Write rows to TSV."""
    if not rows:
        with open(path, "w", encoding="utf-8", newline="") as handle:
            handle.write("")
        return

    columns = list(rows[0].keys())
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Merge TSV transcript rows using TSS/TES/intron windows."
    )

    parser.add_argument(
        "-i", "--inputs",
        nargs="+",
        required=True,
        help="Input TSV files",
    )
    parser.add_argument(
        "-o", "--out",
        required=True,
        help="Output TSV file",
    )
    parser.add_argument(
        "-t", "--tss-window",
        type=int,
        required=True,
        help="TSS window size",
    )
    parser.add_argument(
        "-e", "--tes-window",
        type=int,
        required=True,
        help="TES window size",
    )
    parser.add_argument(
        "-n", "--intron-window",
        type=int,
        required=True,
        help="Intron window size",
    )
    parser.add_argument(
        "-x", "--no-within-source-merge",
        action="store_true",
        help="Do not merge records within the same input TSV",
    )
    parser.add_argument(
        "-1", "--pairwise-exclusive",
        action="store_true",
        help="Apply greedy 1:1 matching between input files",
    )

    args = parser.parse_args()

    all_rows: List[TranscriptRow] = []
    idx = 0

    for path in args.inputs:
        rows = load_tsv(path, start_idx=idx)
        all_rows.extend(rows)
        idx += len(rows)

    if not all_rows:
        raise ValueError("No valid transcript rows could be loaded.")

    uf = UnionFind(len(all_rows))

    edges = build_pairwise_edges(
        rows=all_rows,
        tss_window=args.tss_window,
        tes_window=args.tes_window,
        intron_window=args.intron_window,
        no_within_source_merge=args.no_within_source_merge,
        pairwise_exclusive=args.pairwise_exclusive,
    )

    for a, b in edges:
        uf.union(a, b)

    out_rows = build_cluster_table(all_rows, uf)
    write_tsv(args.out, out_rows)

    print(f"Done: {args.out}")
    print(f"Loaded records: {len(all_rows)}")
    print(f"Created edges: {len(edges)}")
    print(f"Output rows: {len(out_rows)}")


if __name__ == "__main__":
    main()