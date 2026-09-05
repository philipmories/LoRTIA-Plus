#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "Usage: $0 SOURCE_RECOVERY_MANIFEST [OUTPUT_DIR]" >&2
    exit 2
fi

SOURCE_MANIFEST="$1"
OUTDIR="${2:-output}"
PYTHON_BIN="${PYTHON_BIN:-python}"
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="$ROOT_DIR/scripts"

mkdir -p \
    "$OUTDIR/Figure_input" \
    "$OUTDIR/category_and_length_outputs" \
    "$OUTDIR/isoforms_per_gene_outputs" \
    "$OUTDIR/endpoint_support" \
    "$OUTDIR/manifests"

BASIC_MANIFEST="$OUTDIR/manifests/manifest_basic_stats_6annotators.tsv"

"$PYTHON_BIN" "$SCRIPTS/build_basic_stats_manifest.py" \
    -i "$SOURCE_MANIFEST" \
    -o "$BASIC_MANIFEST"

"$PYTHON_BIN" "$SCRIPTS/build_sqanti_basic_stat_table.py" \
    -i "$BASIC_MANIFEST" \
    -o "$OUTDIR/Figure_input/SQANTI3_basic_stat.tsv"

"$PYTHON_BIN" "$SCRIPTS/extract_sqanti_length_and_category_summaries.py" \
    --manifest "$BASIC_MANIFEST" \
    -o "$OUTDIR/category_and_length_outputs"

cp "$OUTDIR/category_and_length_outputs/sqanti_length_boxplot_stats.tsv" \
   "$OUTDIR/Figure_input/sqanti_length_boxplot_stats.tsv"

"$PYTHON_BIN" "$SCRIPTS/summarize_isoforms_per_gene.py" \
    -d "$BASIC_MANIFEST" \
    -o "$OUTDIR/isoforms_per_gene_outputs/isoformsgene.tsv" \
    --category-output "$OUTDIR/isoforms_per_gene_outputs/sqanti_category_stats.tsv"

cp "$OUTDIR/isoforms_per_gene_outputs/isoformsgene.tsv" \
   "$OUTDIR/Figure_input/isoformsgene.tsv"

"$PYTHON_BIN" "$SCRIPTS/summarize_sqanti_endpoint_support.py" \
    -d "$BASIC_MANIFEST" \
    -o "$OUTDIR/endpoint_support/sqanti_endpoint_support.tsv"

echo "[OK] SQANTI3 basic-statistics preprocessing completed."
echo "[OK] Canonical manifest: $BASIC_MANIFEST"
echo "[OK] Plot-ready inputs: $OUTDIR/Figure_input"
