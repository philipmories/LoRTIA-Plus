#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
evaluate_transcripts_tsv_vs_gff3_clean.py

Benchmark transcript predictions from GFF3 files against a reference TSV transcript set.

Reference TSV expected columns:
- isoform
- chromosome
- start
- end
- strand
- exon_composition
- reads

Prediction input:
- a directory containing GFF3 files with mRNA and exon features

Comparison criteria:
- seqid / chromosome
- strand
- TSS window
- TES window
- intron chain window

Outputs:
- <out>_matches.tsv
- <out>_stats.tsv
- <out>_false_positives.tsv
- <out>_reference_not_found.tsv

Additional notes:
- one reference transcript may match multiple predictions, but still counts as one TP
- predicted read counts are taken primarily from the mRNA score field
- if score is not numeric, the script tries these attributes:
  read_count, reads, count, support, coverage
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import sys
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
from collections import defaultdict


# ----------------------------
# Data classes
# ----------------------------

@dataclass(frozen=True)
class TranscriptRef:
    ref_id: str
    isoform: str
    seqid: str
    strand: str
    start: int
    end: int
    exon_comp: str
    reads: int
    exons: Tuple[Tuple[int, int], ...]
    introns: Tuple[Tuple[int, int], ...]


@dataclass(frozen=True)
class TranscriptPred:
    pred_id: str
    seqid: str
    strand: str
    start: int
    end: int
    reads: float
    exons: Tuple[Tuple[int, int], ...]
    introns: Tuple[Tuple[int, int], ...]


@dataclass(frozen=True)
class FalsePositivePred:
    gff3_file: str
    pred_id: str
    seqid: str
    strand: str
    start: int
    end: int
    reads: float
    exons: Tuple[Tuple[int, int], ...]
    introns: Tuple[Tuple[int, int], ...]


@dataclass
class FalsePositiveCluster:
    seqid: str
    strand: str
    start: int
    end: int
    exons: Tuple[Tuple[int, int], ...]
    introns: Tuple[Tuple[int, int], ...]
    members: List[FalsePositivePred]


# ----------------------------
# Helpers
# ----------------------------

def norm_strand(value: str) -> str:
    value = (value or ".").strip()
    return value if value in {"+", "-", "."} else value


def safe_float_div(numerator: float, denominator: float) -> float:
    return (numerator / denominator) if denominator != 0 else 0.0


def f1_score(precision: float, recall: float) -> float:
    return safe_float_div(2 * precision * recall, precision + recall)


def resolve_column(possible_names: List[str], fieldnames: List[str]) -> Optional[str]:
    lower_map = {name.lower(): name for name in fieldnames}
    for name in possible_names:
        if name in fieldnames:
            return name
        if name.lower() in lower_map:
            return lower_map[name.lower()]
    return None


def parse_int(value) -> Optional[int]:
    if value is None:
        return None
    text = str(value).strip().replace(",", "")
    if not text or text.lower() == "nan":
        return None
    try:
        return int(float(text))
    except Exception:
        return None


def parse_float(value) -> Optional[float]:
    if value is None:
        return None
    text = str(value).strip().replace(",", "")
    if not text or text.lower() == "nan" or text == ".":
        return None
    try:
        return float(text)
    except Exception:
        return None


def parse_reads(value) -> int:
    parsed = parse_int(value)
    return 0 if parsed is None else parsed


def parse_exon_composition(exon_comp: str) -> List[Tuple[int, int]]:
    """Parse exon composition string into exon tuples."""
    if exon_comp is None:
        return []

    text = str(exon_comp).strip()
    if not text or text.lower() == "nan":
        return []

    exons: List[Tuple[int, int]] = []
    for part in text.split(";"):
        part = part.strip()
        if not part or "-" not in part:
            continue

        a, b = part.split("-", 1)
        start = parse_int(a)
        end = parse_int(b)
        if start is None or end is None:
            continue

        exon_start, exon_end = (start, end) if start <= end else (end, start)
        exons.append((exon_start, exon_end))

    exons.sort()
    return exons


def exons_to_string(exons: Tuple[Tuple[int, int], ...]) -> str:
    if not exons:
        return "."
    return ";".join(f"{start}-{end}" for start, end in exons)


def exons_to_introns(exons: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Convert exon coordinates into intron intervals."""
    if len(exons) <= 1:
        return []

    introns: List[Tuple[int, int]] = []
    exons_sorted = sorted(exons)

    for i in range(len(exons_sorted) - 1):
        left = exons_sorted[i]
        right = exons_sorted[i + 1]
        intron_start = left[1] + 1
        intron_end = right[0] - 1
        introns.append((intron_start, intron_end))

    return introns


def tss_tes(start: int, end: int, strand: str) -> Tuple[int, int]:
    if strand == "-":
        return end, start
    return start, end


def parse_gff3_attrs(attr: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    if not attr or attr == ".":
        return attrs

    for part in attr.split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
            attrs[key.strip()] = value.strip()

    return attrs


def get_pred_reads(score: str, attrs: Dict[str, str]) -> float:
    """
    Extract predicted read support:
    1) mRNA score field
    2) attributes: read_count, reads, count, support, coverage
    """
    score_val = parse_float(score)
    if score_val is not None:
        return score_val

    for key in ["read_count", "reads", "count", "support", "coverage"]:
        if key in attrs:
            value = parse_float(attrs[key])
            if value is not None:
                return value

    return 0.0


def fmt_reads(value: float) -> str:
    if abs(value - round(value)) < 1e-9:
        return str(int(round(value)))
    return f"{value:.6f}".rstrip("0").rstrip(".")


def transcript_match(
    ref: TranscriptRef,
    pred: TranscriptPred,
    tss_window: int,
    tes_window: int,
    intron_window: int,
) -> bool:
    """Check whether a prediction matches a reference transcript."""
    if ref.seqid != pred.seqid:
        return False
    if ref.strand != pred.strand:
        return False

    ref_tss, ref_tes = tss_tes(ref.start, ref.end, ref.strand)
    pred_tss, pred_tes = tss_tes(pred.start, pred.end, pred.strand)

    if abs(ref_tss - pred_tss) > tss_window:
        return False
    if abs(ref_tes - pred_tes) > tes_window:
        return False

    if len(ref.introns) != len(pred.introns):
        return False

    for (ref_start, ref_end), (pred_start, pred_end) in zip(ref.introns, pred.introns):
        if abs(ref_start - pred_start) > intron_window:
            return False
        if abs(ref_end - pred_end) > intron_window:
            return False

    return True


def pred_pred_match(
    a: FalsePositivePred,
    b: FalsePositivePred,
    tss_window: int,
    tes_window: int,
    intron_window: int,
) -> bool:
    """Check whether two false-positive predictions belong to the same FP cluster."""
    if a.seqid != b.seqid:
        return False
    if a.strand != b.strand:
        return False

    a_tss, a_tes = tss_tes(a.start, a.end, a.strand)
    b_tss, b_tes = tss_tes(b.start, b.end, b.strand)

    if abs(a_tss - b_tss) > tss_window:
        return False
    if abs(a_tes - b_tes) > tes_window:
        return False

    if len(a.introns) != len(b.introns):
        return False

    for (a_start, a_end), (b_start, b_end) in zip(a.introns, b.introns):
        if abs(a_start - b_start) > intron_window:
            return False
        if abs(a_end - b_end) > intron_window:
            return False

    return True


# ----------------------------
# Reference TSV parsing
# ----------------------------

def parse_reference_tsv(path: str) -> List[TranscriptRef]:
    """Load reference transcripts from TSV."""
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []

        col_isoform = resolve_column(["isoform", "transcript_id", "id"], fieldnames)
        col_seqid = resolve_column(["chromosome", "seqid", "chrom", "chr"], fieldnames)
        col_start = resolve_column(["start"], fieldnames)
        col_end = resolve_column(["end"], fieldnames)
        col_strand = resolve_column(["strand"], fieldnames)
        col_exon = resolve_column(["exon_composition", "exon_comp"], fieldnames)
        col_reads = resolve_column(["reads", "read_count"], fieldnames)

        required = {
            "isoform": col_isoform,
            "chromosome": col_seqid,
            "start": col_start,
            "end": col_end,
            "strand": col_strand,
            "exon_composition": col_exon,
            "reads": col_reads,
        }

        missing = [label for label, column in required.items() if column is None]
        if missing:
            raise ValueError(
                f"Missing reference TSV column(s): {missing}. Available columns: {fieldnames}"
            )

        refs: List[TranscriptRef] = []
        auto_index = 0

        for row in reader:
            auto_index += 1

            isoform = str(row[col_isoform]).strip()
            seqid = str(row[col_seqid]).strip()
            strand = norm_strand(row[col_strand])

            start = parse_int(row[col_start])
            end = parse_int(row[col_end])
            if start is None or end is None:
                continue
            if start > end:
                start, end = end, start

            exon_comp = str(row[col_exon]).strip()
            reads = parse_reads(row[col_reads])

            exons = parse_exon_composition(exon_comp)
            if not exons:
                exons = [(start, end)]

            introns = exons_to_introns(exons)
            ref_id = isoform if isoform else f"REF_{auto_index:06d}"

            refs.append(
                TranscriptRef(
                    ref_id=ref_id,
                    isoform=isoform,
                    seqid=seqid,
                    strand=strand,
                    start=start,
                    end=end,
                    exon_comp=exon_comp,
                    reads=reads,
                    exons=tuple(exons),
                    introns=tuple(introns),
                )
            )

    return refs


# ----------------------------
# GFF3 parsing
# ----------------------------

def parse_gff3_transcripts(path: str) -> List[TranscriptPred]:
    """Load transcript predictions from a GFF3 file."""
    mrnas: Dict[str, Dict[str, object]] = {}
    exons_by_parent: Dict[str, List[Tuple[int, int]]] = defaultdict(list)

    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            cols = line.split("\t")
            if len(cols) != 9:
                continue

            seqid, _, feature_type, start_raw, end_raw, score, strand, _, attrs = cols
            start = parse_int(start_raw)
            end = parse_int(end_raw)
            if start is None or end is None:
                continue
            if start > end:
                start, end = end, start

            strand = norm_strand(strand)
            attr_dict = parse_gff3_attrs(attrs)

            if feature_type == "mRNA":
                transcript_id = attr_dict.get("ID")
                if not transcript_id:
                    continue

                mrnas[transcript_id] = {
                    "seqid": seqid,
                    "strand": strand,
                    "start": start,
                    "end": end,
                    "reads": get_pred_reads(score, attr_dict),
                }

            elif feature_type == "exon":
                parent = attr_dict.get("Parent")
                if not parent:
                    continue

                for transcript_id in parent.split(","):
                    transcript_id = transcript_id.strip()
                    if transcript_id:
                        exons_by_parent[transcript_id].append((start, end))

    preds: List[TranscriptPred] = []
    for transcript_id, mrna in mrnas.items():
        exons = sorted(exons_by_parent.get(transcript_id, []))
        if not exons:
            exons = [(int(mrna["start"]), int(mrna["end"]))]

        introns = exons_to_introns(exons)

        preds.append(
            TranscriptPred(
                pred_id=transcript_id,
                seqid=str(mrna["seqid"]),
                strand=str(mrna["strand"]),
                start=int(mrna["start"]),
                end=int(mrna["end"]),
                reads=float(mrna["reads"]),
                exons=tuple(exons),
                introns=tuple(introns),
            )
        )

    return preds


# ----------------------------
# FP clustering
# ----------------------------

def cluster_false_positives(
    fps: List[FalsePositivePred],
    tss_window: int,
    tes_window: int,
    intron_window: int,
) -> List[FalsePositiveCluster]:
    """Cluster unmatched predictions across files."""
    if not fps:
        return []

    grouped: Dict[Tuple[str, str], List[FalsePositivePred]] = defaultdict(list)
    for fp in fps:
        grouped[(fp.seqid, fp.strand)].append(fp)

    clusters: List[FalsePositiveCluster] = []

    for (seqid, strand), items in grouped.items():
        items_sorted = sorted(items, key=lambda x: tss_tes(x.start, x.end, x.strand)[0])
        used = [False] * len(items_sorted)

        for i, seed in enumerate(items_sorted):
            if used[i]:
                continue

            members: List[FalsePositivePred] = []
            for j, candidate in enumerate(items_sorted):
                if used[j]:
                    continue
                if pred_pred_match(seed, candidate, tss_window, tes_window, intron_window):
                    used[j] = True
                    members.append(candidate)

            clusters.append(
                FalsePositiveCluster(
                    seqid=seqid,
                    strand=strand,
                    start=seed.start,
                    end=seed.end,
                    exons=seed.exons,
                    introns=seed.introns,
                    members=members,
                )
            )

    clusters.sort(key=lambda c: (c.seqid, c.strand, c.start, c.end))
    return clusters


# ----------------------------
# Evaluation
# ----------------------------

def evaluate_one_file(
    refs: List[TranscriptRef],
    preds: List[TranscriptPred],
    tss_window: int,
    tes_window: int,
    intron_window: int,
) -> Tuple[Dict[str, List[str]], Dict[str, List[float]], int, int, int, List[TranscriptPred]]:
    """
    Evaluate one GFF3 file against the reference set.
    Returns matches, matched read counts, TP/FP/FN counts, and unmatched predictions.
    """
    ref_to_preds: Dict[str, List[str]] = {ref.ref_id: [] for ref in refs}
    ref_to_pred_reads: Dict[str, List[float]] = {ref.ref_id: [] for ref in refs}
    pred_matched: Dict[str, bool] = {pred.pred_id: False for pred in preds}

    refs_by_key: Dict[Tuple[str, str], List[TranscriptRef]] = defaultdict(list)
    preds_by_key: Dict[Tuple[str, str], List[TranscriptPred]] = defaultdict(list)

    for ref in refs:
        refs_by_key[(ref.seqid, ref.strand)].append(ref)
    for pred in preds:
        preds_by_key[(pred.seqid, pred.strand)].append(pred)

    for key, ref_list in refs_by_key.items():
        pred_list = preds_by_key.get(key, [])
        if not pred_list:
            continue

        for ref in ref_list:
            for pred in pred_list:
                if transcript_match(ref, pred, tss_window, tes_window, intron_window):
                    ref_to_preds[ref.ref_id].append(pred.pred_id)
                    ref_to_pred_reads[ref.ref_id].append(pred.reads)
                    pred_matched[pred.pred_id] = True

    tp = sum(1 for ref_id, matches in ref_to_preds.items() if matches)
    fn = sum(1 for ref_id, matches in ref_to_preds.items() if not matches)
    fp_preds = [pred for pred in preds if not pred_matched[pred.pred_id]]
    fp = len(fp_preds)

    return ref_to_preds, ref_to_pred_reads, tp, fp, fn, fp_preds


# ----------------------------
# Main
# ----------------------------

def main() -> int:
    parser = argparse.ArgumentParser(
        description="Benchmark transcript predictions from GFF3 files against a reference TSV."
    )
    parser.add_argument("--ref", required=True, help="Reference TSV file")
    parser.add_argument("--gff3-dir", required=True, help="Directory containing GFF3 files")
    parser.add_argument("--glob", default="*.gff3", help="GFF3 glob pattern (default: *.gff3)")

    parser.add_argument("--tss-window", type=int, default=10, help="TSS window size (default: 10)")
    parser.add_argument("--tes-window", type=int, default=10, help="TES window size (default: 10)")
    parser.add_argument(
        "--intron_window", "--intron-window",
        dest="intron_window",
        type=int,
        default=10,
        help="Intron window size (default: 10)",
    )

    parser.add_argument("--out", required=True, help="Output prefix")

    args = parser.parse_args()

    if args.tss_window < 0 or args.tes_window < 0 or args.intron_window < 0:
        print("[ERROR] Window sizes cannot be negative.", file=sys.stderr)
        return 2

    refs = parse_reference_tsv(args.ref)

    gff3_paths = sorted(glob.glob(os.path.join(args.gff3_dir, args.glob)))
    if not gff3_paths:
        print("[ERROR] No GFF3 files were found.", file=sys.stderr)
        return 2

    per_file_ref_matches: Dict[str, Dict[str, List[str]]] = {}
    per_file_ref_reads: Dict[str, Dict[str, List[float]]] = {}
    stats_rows: List[Tuple[str, int, int, int, float, float, float]] = []
    false_positive_rows: List[FalsePositivePred] = []

    ref_found_anywhere: Dict[str, bool] = {ref.ref_id: False for ref in refs}

    for path in gff3_paths:
        filename = os.path.basename(path)
        preds = parse_gff3_transcripts(path)

        ref_to_preds, ref_to_pred_reads, tp, fp, fn, fp_preds = evaluate_one_file(
            refs=refs,
            preds=preds,
            tss_window=args.tss_window,
            tes_window=args.tes_window,
            intron_window=args.intron_window,
        )

        for ref_id, pred_ids in ref_to_preds.items():
            if pred_ids:
                ref_found_anywhere[ref_id] = True

        precision = safe_float_div(tp, tp + fp)
        recall = safe_float_div(tp, tp + fn)
        f1 = f1_score(precision, recall)

        per_file_ref_matches[filename] = ref_to_preds
        per_file_ref_reads[filename] = ref_to_pred_reads
        stats_rows.append((filename, tp, fp, fn, precision, recall, f1))

        for pred in fp_preds:
            false_positive_rows.append(
                FalsePositivePred(
                    gff3_file=filename,
                    pred_id=pred.pred_id,
                    seqid=pred.seqid,
                    strand=pred.strand,
                    start=pred.start,
                    end=pred.end,
                    reads=pred.reads,
                    exons=pred.exons,
                    introns=pred.introns,
                )
            )

    matches_path = args.out + "_matches.tsv"
    stats_path = args.out + "_stats.tsv"
    false_positives_path = args.out + "_false_positives.tsv"
    ref_not_found_path = args.out + "_reference_not_found.tsv"

    file_names = [os.path.basename(path) for path in gff3_paths]
    refs_sorted = sorted(refs, key=lambda r: (r.seqid, r.strand, r.start, r.end, r.ref_id))

    # ---------------- matches.tsv ----------------
    with open(matches_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        header = [
            "ref_id",
            "isoform",
            "seqid",
            "strand",
            "start",
            "end",
            "exon_composition",
            "reads",
            "intron_count",
        ]
        header += [f"{fn}_match" for fn in file_names]
        header += [f"{fn}_pred_ids" for fn in file_names]
        header += [f"{fn}_pred_reads" for fn in file_names]
        writer.writerow(header)

        for ref in refs_sorted:
            row = [
                ref.ref_id,
                ref.isoform,
                ref.seqid,
                ref.strand,
                ref.start,
                ref.end,
                ref.exon_comp,
                ref.reads,
                len(ref.introns),
            ]

            for fn in file_names:
                matches = per_file_ref_matches.get(fn, {}).get(ref.ref_id, [])
                row.append(1 if matches else 0)

            for fn in file_names:
                matches = per_file_ref_matches.get(fn, {}).get(ref.ref_id, [])
                row.append(",".join(matches) if matches else ".")

            for fn in file_names:
                values = per_file_ref_reads.get(fn, {}).get(ref.ref_id, [])
                row.append(fmt_reads(sum(values)) if values else "0")

            writer.writerow(row)

    # ---------------- stats.tsv ----------------
    with open(stats_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["gff3_file", "TP", "FP", "FN", "precision", "recall", "F1"])

        for filename, tp, fp, fn, precision, recall, f1 in stats_rows:
            writer.writerow([
                filename,
                tp,
                fp,
                fn,
                f"{precision:.6f}",
                f"{recall:.6f}",
                f"{f1:.6f}",
            ])

    # ---------------- false_positives.tsv ----------------
    fp_clusters = cluster_false_positives(
        false_positive_rows,
        tss_window=args.tss_window,
        tes_window=args.tes_window,
        intron_window=args.intron_window,
    )

    with open(false_positives_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")

        header = [
            "ID",
            "seqid",
            "strand",
            "start",
            "end",
            "tss",
            "tes",
            "exon_composition",
            "intron_count",
            "members_total",
            "files_present",
        ]
        header += [f"{fn}_present" for fn in file_names]
        header += [f"{fn}_count" for fn in file_names]
        header += [f"{fn}_pred_ids" for fn in file_names]
        header += [f"{fn}_pred_reads" for fn in file_names]
        writer.writerow(header)

        for i, cluster in enumerate(fp_clusters, start=1):
            cluster_id = f"FP_{i:06d}"
            tss, tes = tss_tes(cluster.start, cluster.end, cluster.strand)
            exon_comp = exons_to_string(cluster.exons)

            file_to_ids: Dict[str, List[str]] = {fn: [] for fn in file_names}
            file_to_reads: Dict[str, List[float]] = {fn: [] for fn in file_names}

            for member in cluster.members:
                file_to_ids[member.gff3_file].append(member.pred_id)
                file_to_reads[member.gff3_file].append(member.reads)

            files_present = [fn for fn in file_names if file_to_ids[fn]]

            row = [
                cluster_id,
                cluster.seqid,
                cluster.strand,
                cluster.start,
                cluster.end,
                tss,
                tes,
                exon_comp,
                len(cluster.introns),
                len(cluster.members),
                ",".join(files_present) if files_present else ".",
            ]
            row += [1 if file_to_ids[fn] else 0 for fn in file_names]
            row += [len(file_to_ids[fn]) for fn in file_names]
            row += [",".join(file_to_ids[fn]) if file_to_ids[fn] else "." for fn in file_names]
            row += [fmt_reads(sum(file_to_reads[fn])) if file_to_reads[fn] else "0" for fn in file_names]

            writer.writerow(row)

    # ---------------- reference_not_found.tsv ----------------
    with open(ref_not_found_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "ref_id",
            "isoform",
            "seqid",
            "strand",
            "start",
            "end",
            "exon_composition",
            "reads",
            "intron_count",
        ])

        for ref in refs_sorted:
            if not ref_found_anywhere.get(ref.ref_id, False):
                writer.writerow([
                    ref.ref_id,
                    ref.isoform,
                    ref.seqid,
                    ref.strand,
                    ref.start,
                    ref.end,
                    ref.exon_comp,
                    ref.reads,
                    len(ref.introns),
                ])

    print(
        f"[OK] refs={len(refs)} gff3_files={len(gff3_paths)} "
        f"tss_window={args.tss_window} tes_window={args.tes_window} intron_window={args.intron_window}"
    )
    print(f"[OK] matches             : {matches_path}")
    print(f"[OK] stats               : {stats_path}")
    print(f"[OK] false positives     : {false_positives_path}")
    print(f"[OK] reference not found : {ref_not_found_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())