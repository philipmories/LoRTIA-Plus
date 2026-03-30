#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compare TSS, TES, and intron features from multiple GFF/GFF3 files
against a reference GFF using a positional tolerance window.

The script:
- loads canonical TSS/TES/intron features from a main GFF
- matches features from all other input files to the reference set
- sums matching scores per sample
- writes a TSV summary and a merged GFF3 output
"""

from __future__ import annotations

import argparse
import os
import sys
from collections import defaultdict
from typing import Dict, Iterable, List, Optional, Tuple

ALLOWED_TYPES = {"TSS", "TES", "intron"}


def parse_gff_attributes(attr_field: str) -> Dict[str, str]:
    """Parse GFF-style attributes into a dictionary."""
    attributes: Dict[str, str] = {}
    if not attr_field or attr_field == ".":
        return attributes

    for part in [p.strip() for p in attr_field.split(";") if p.strip()]:
        if "=" in part:
            key, value = part.split("=", 1)
            attributes[key.strip()] = value.strip()
        else:
            tokens = part.split()
            if len(tokens) >= 2:
                attributes[tokens[0].strip()] = " ".join(tokens[1:]).strip().strip('"')

    return attributes


def iter_gff_records(path: str) -> Iterable[Dict[str, object]]:
    """Yield parsed records from a GFF/GFF3 file, including .gz input."""
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


def load_main_features(main_gff_path: str) -> Tuple[
    List[Dict[str, object]],
    Dict[Tuple[str, str, str], List[Dict[str, object]]],
]:
    """
    Load canonical TSS/TES/intron features from the reference GFF.

    Returns:
    - canonical_features: ordered list of reference features
    - canonical_index: keyed by (seqid, strand, type)
    """
    raw_features: List[Dict[str, object]] = []

    for record in iter_gff_records(main_gff_path):
        if record["type"] in ALLOWED_TYPES:
            raw_features.append(record)

    raw_features.sort(
        key=lambda record: (
            record["seqid"],
            record["strand"],
            record["type"],
            record["start"],
            record["end"],
        )
    )

    canonical_features: List[Dict[str, object]] = []
    canonical_index: Dict[Tuple[str, str, str], List[Dict[str, object]]] = defaultdict(list)

    for idx, record in enumerate(raw_features, start=1):
        feature = {
            "feat_id": f"FT_{idx:06d}",
            "type": record["type"],
            "seqid": record["seqid"],
            "strand": record["strand"],
            "start": record["start"],
            "end": record["end"],
        }
        canonical_features.append(feature)
        canonical_index[(feature["seqid"], feature["strand"], feature["type"])] .append(feature)  # type: ignore

    return canonical_features, canonical_index


def feature_distance(
    record: Dict[str, object],
    canonical_feature: Dict[str, object],
    pos_window: int,
) -> Optional[int]:
    """
    Return the distance between an input feature and a reference feature.

    Returns None if the feature does not match within the positional window.

    Rules:
    - TSS/TES: compare only the start coordinate
    - intron: both start and end must be within the window;
      distance is the sum of the two endpoint differences
    """
    feature_type = record["type"]
    start = int(record["start"])
    end = int(record["end"])
    canonical_start = int(canonical_feature["start"])
    canonical_end = int(canonical_feature["end"])

    if feature_type in ("TSS", "TES"):
        distance = abs(start - canonical_start)
        return distance if distance <= pos_window else None

    start_diff = abs(start - canonical_start)
    end_diff = abs(end - canonical_end)
    if start_diff <= pos_window and end_diff <= pos_window:
        return start_diff + end_diff
    return None


def accumulate_counts_for_sample(
    sample_path: str,
    sample_name: str,
    canonical_index: Dict[Tuple[str, str, str], List[Dict[str, object]]],
    pos_window: int,
    counts: Dict[str, Dict[str, float]],
) -> None:
    """
    Match all valid features in one sample to the canonical reference set.

    If one feature matches multiple canonical features within the window,
    its value is added to all of them.
    """
    for record in iter_gff_records(sample_path):
        feature_type = record["type"]
        if feature_type not in ALLOWED_TYPES:
            continue

        key = (record["seqid"], record["strand"], feature_type)  # type: ignore
        candidates = canonical_index.get(key)
        if not candidates:
            continue

        matched_features: List[Dict[str, object]] = []
        for canonical_feature in candidates:
            if feature_distance(record, canonical_feature, pos_window) is not None:
                matched_features.append(canonical_feature)

        if not matched_features:
            continue

        score = record["score"]
        value = 1.0 if score is None else float(score)  # type: ignore[arg-type]

        for canonical_feature in matched_features:
            feature_id = str(canonical_feature["feat_id"])
            counts.setdefault(feature_id, {})
            counts[feature_id][sample_name] = counts[feature_id].get(sample_name, 0.0) + value


def int_if_integer(value: float):
    """Convert a float to int if it has no fractional part."""
    return int(value) if float(value).is_integer() else float(value)


def format_plain_number(value: float, decimals: int = 6) -> str:
    """Format numbers without unnecessary trailing zeros."""
    if float(value).is_integer():
        return str(int(value))
    text = f"{value:.{decimals}f}".rstrip("0").rstrip(".")
    return text if text else "0"


def write_tsv(
    tsv_path: str,
    canonical_features: List[Dict[str, object]],
    counts: Dict[str, Dict[str, float]],
    sample_names: List[str],
) -> None:
    """Write the reference feature table as TSV."""
    import csv

    headers = ["feature_id", "type", "chromosome", "start", "end", "strand"] + sample_names
    with open(tsv_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(headers)

        for feature in canonical_features:
            feature_id = str(feature["feat_id"])
            row = [
                feature_id,
                feature["type"],
                feature["seqid"],
                feature["start"],
                feature["end"],
                feature["strand"],
            ] + [
                int_if_integer(counts.get(feature_id, {}).get(sample_name, 0.0))
                for sample_name in sample_names
            ]
            writer.writerow(row)


def write_gff3(
    gff_path: str,
    canonical_features: List[Dict[str, object]],
    counts: Dict[str, Dict[str, float]],
    sample_names: List[str],
    write_sample_attrs: bool = False,
) -> None:
    """Write the reference feature set as GFF3 with merged sample scores."""
    with open(gff_path, "w", encoding="utf-8", newline="") as handle:
        handle.write("##gff-version 3\n")

        for feature in canonical_features:
            feature_id = str(feature["feat_id"])
            seqid = str(feature["seqid"])
            strand = str(feature["strand"])
            start = int(feature["start"])
            end = int(feature["end"])
            feature_type = str(feature["type"])

            total = sum(float(v) for v in counts.get(feature_id, {}).values())
            score_str = "." if total == 0 else format_plain_number(total)

            attrs = [
                f"ID={feature_id}",
                f"Name={feature_id}",
                f"feature_type={feature_type}",
            ]

            if write_sample_attrs:
                sample_attr = ";".join(
                    f"{sample_name}={format_plain_number(counts.get(feature_id, {}).get(sample_name, 0.0))}"
                    for sample_name in sample_names
                )
                if sample_attr:
                    attrs.append(sample_attr)

            handle.write(
                "\t".join(
                    [
                        seqid,
                        "compare",
                        feature_type,
                        str(start),
                        str(end),
                        score_str,
                        strand,
                        ".",
                        ";".join(attrs),
                    ]
                )
                + "\n"
            )


def collect_gff_files(indir: str, pattern: str = "auto") -> List[str]:
    """Recursively collect non-empty input GFF/GFF3 files."""
    patt = pattern.lower()

    if patt == "auto":
        extensions: Tuple[str, ...] = (".gff3", ".gff", ".gff3.gz", ".gff.gz")
    elif patt in (".gff3", "gff3"):
        extensions = (".gff3", ".gff3.gz")
    elif patt in (".gff", "gff"):
        extensions = (".gff", ".gff.gz")
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
        description="Compare TSS, TES, and intron features against a reference GFF."
    )
    parser.add_argument(
        "-i",
        "--indir",
        required=True,
        help="Input directory containing comparison GFF/GFF3 files (searched recursively).",
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
        help="File extension filter: auto | .gff3 | .gff (default: auto).",
    )
    parser.add_argument(
        "--main-gff",
        required=True,
        help="Reference GFF/GFF3 file containing TSS/TES/intron features.",
    )
    parser.add_argument(
        "--pos-window",
        type=int,
        default=0,
        help=(
            "Positional tolerance window in bp. For TSS/TES only the start position is checked; "
            "for introns both boundaries must fall within the window."
        ),
    )
    parser.add_argument(
        "--write-sample-attrs",
        action="store_true",
        help="Include per-sample counts in GFF3 attributes.",
    )
    args = parser.parse_args()

    indir = os.path.abspath(args.indir)
    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)

    base_name = os.path.basename(outdir.rstrip("/"))
    tsv_path = os.path.join(outdir, f"{base_name}.tsv")
    gff_path = os.path.join(outdir, f"{base_name}.gff3")

    main_gff_path = os.path.abspath(args.main_gff)
    if not os.path.isfile(main_gff_path):
        sys.exit(f"Reference GFF file not found: {main_gff_path}")

    main_sample_name = os.path.basename(main_gff_path)
    if main_sample_name.endswith(".gz"):
        main_sample_name = main_sample_name[:-3]

    print(f"Reference GFF: {main_gff_path}")
    canonical_features, canonical_index = load_main_features(main_gff_path)

    if not canonical_features:
        sys.exit("No TSS/TES/intron features were found in the reference GFF.")

    sample_names: List[str] = [main_sample_name]
    counts: Dict[str, Dict[str, float]] = {}

    accumulate_counts_for_sample(
        main_gff_path,
        main_sample_name,
        canonical_index,
        args.pos_window,
        counts,
    )

    files = collect_gff_files(indir, args.pattern)
    files = [path for path in files if os.path.abspath(path) != main_gff_path]

    for path in files:
        sample_name = os.path.basename(path)
        if sample_name.endswith(".gz"):
            sample_name = sample_name[:-3]
        sample_names.append(sample_name)
        print(f"Processing: {path}")
        accumulate_counts_for_sample(
            path,
            sample_name,
            canonical_index,
            args.pos_window,
            counts,
        )

    print(f"Total samples including reference: {len(sample_names)}")

    write_tsv(tsv_path, canonical_features, counts, sample_names)
    write_gff3(
        gff_path,
        canonical_features,
        counts,
        sample_names,
        write_sample_attrs=args.write_sample_attrs,
    )

    print(f"Done: {tsv_path}")
    print(f"Done: {gff_path}")
    print(f"Reference feature count: {len(canonical_features)}")
    print(f"Samples: {', '.join(sample_names)}")


if __name__ == "__main__":
    main()
