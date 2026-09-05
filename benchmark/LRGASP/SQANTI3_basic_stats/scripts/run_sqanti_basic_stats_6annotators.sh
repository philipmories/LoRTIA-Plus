#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY_DIR="$ROOT_DIR/scripts/preprocessing"
R_DIR="$ROOT_DIR/scripts/plotting"
DEFAULT_MANIFEST="$ROOT_DIR/manifest/manifest_reference.tsv"
DEFAULT_INPUT_DIR="$ROOT_DIR/input"
DEFAULT_OUTPUT_DIR="$ROOT_DIR/output"
PYTHON_BIN="${PYTHON_BIN:-python}"
RSCRIPT_BIN="${RSCRIPT_BIN:-Rscript}"

usage() {
    cat <<'USAGE'
SQANTI3 basic-statistics workflow (six annotators)

FULL REGENERATION FROM SQANTI3 classification files:
  ./run_sqanti_basic_stats_6annotators.sh [MANIFEST] [OUTPUT_DIR]

  Defaults:
    MANIFEST   = manifest/manifest_reference.tsv
    OUTPUT_DIR = output

PLOT-ONLY REPRODUCTION FROM ARCHIVED TSV INPUTS:
  ./run_sqanti_basic_stats_6annotators.sh --plot-only [INPUT_DIR] [OUTPUT_DIR]

  Defaults:
    INPUT_DIR  = input
    OUTPUT_DIR = output

Examples:
  ./run_sqanti_basic_stats_6annotators.sh
  ./run_sqanti_basic_stats_6annotators.sh manifest/manifest_reference.tsv results/full_run
  ./run_sqanti_basic_stats_6annotators.sh --plot-only
  ./run_sqanti_basic_stats_6annotators.sh --plot-only input results/plot_only
USAGE
}

run_plots() {
    local outdir="$1"
    mkdir -p "$outdir/logs"
    "$RSCRIPT_BIN" "$R_DIR/run_sqanti_basic_stats_all.R" "$outdir" \
        2>&1 | tee "$outdir/logs/plotting.log"
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

if [[ "${1:-}" == "--plot-only" ]]; then
    INPUT_DIR="${2:-$DEFAULT_INPUT_DIR}"
    OUTDIR="${3:-$DEFAULT_OUTPUT_DIR}"

    mkdir -p "$OUTDIR/Figure_input" "$OUTDIR/logs"

    required_inputs=(
        SQANTI3_basic_stat.tsv
        sqanti_length_boxplot_stats.tsv
        isoformsgene.tsv
    )

    for f in "${required_inputs[@]}"; do
        if [[ ! -f "$INPUT_DIR/$f" ]]; then
            echo "[ERROR] Missing archived plot input: $INPUT_DIR/$f" >&2
            exit 1
        fi
        cp "$INPUT_DIR/$f" "$OUTDIR/Figure_input/$f"
    done

    echo "[1/1] Plotting from archived Figure_input TSV files..."
    run_plots "$OUTDIR"
    echo "[OK] Plot-only reproduction completed: $OUTDIR"
    exit 0
fi

SOURCE_MANIFEST="${1:-$DEFAULT_MANIFEST}"
OUTDIR="${2:-$DEFAULT_OUTPUT_DIR}"

if [[ ! -f "$SOURCE_MANIFEST" ]]; then
    echo "[ERROR] Manifest not found: $SOURCE_MANIFEST" >&2
    exit 1
fi

mkdir -p \
    "$OUTDIR/Figure_input" \
    "$OUTDIR/category_and_length_outputs" \
    "$OUTDIR/isoforms_per_gene_outputs" \
    "$OUTDIR/endpoint_support" \
    "$OUTDIR/manifests" \
    "$OUTDIR/logs"

BASIC_MANIFEST="$OUTDIR/manifests/manifest_basic_stats_6annotators.tsv"

echo "[1/6] Building and validating the canonical 6 x 5 x 3 manifest..."
"$PYTHON_BIN" "$PY_DIR/build_basic_stats_manifest.py" \
    -i "$SOURCE_MANIFEST" \
    -o "$BASIC_MANIFEST" \
    2>&1 | tee "$OUTDIR/logs/01_build_manifest.log"

echo "[2/6] Building the SQANTI3 basic-statistics table..."
"$PYTHON_BIN" "$PY_DIR/build_sqanti_basic_stat_table.py" \
    -i "$BASIC_MANIFEST" \
    -o "$OUTDIR/Figure_input/SQANTI3_basic_stat.tsv" \
    2>&1 | tee "$OUTDIR/logs/02_basic_stat.log"

echo "[3/6] Summarizing transcript lengths and SQANTI3 categories..."
"$PYTHON_BIN" "$PY_DIR/extract_sqanti_length_and_category_summaries.py" \
    --manifest "$BASIC_MANIFEST" \
    -o "$OUTDIR/category_and_length_outputs" \
    2>&1 | tee "$OUTDIR/logs/03_length_category.log"
cp "$OUTDIR/category_and_length_outputs/sqanti_length_boxplot_stats.tsv" \
   "$OUTDIR/Figure_input/sqanti_length_boxplot_stats.tsv"

echo "[4/6] Summarizing isoforms per gene..."
"$PYTHON_BIN" "$PY_DIR/summarize_isoforms_per_gene.py" \
    -d "$BASIC_MANIFEST" \
    -o "$OUTDIR/isoforms_per_gene_outputs/isoformsgene.tsv" \
    --category-output "$OUTDIR/isoforms_per_gene_outputs/sqanti_category_stats.tsv" \
    2>&1 | tee "$OUTDIR/logs/04_isoforms_per_gene.log"
cp "$OUTDIR/isoforms_per_gene_outputs/isoformsgene.tsv" \
   "$OUTDIR/Figure_input/isoformsgene.tsv"

echo "[5/6] Generating the auxiliary endpoint-support summary..."
"$PYTHON_BIN" "$PY_DIR/summarize_sqanti_endpoint_support.py" \
    -d "$BASIC_MANIFEST" \
    -o "$OUTDIR/endpoint_support/sqanti_endpoint_support.tsv" \
    2>&1 | tee "$OUTDIR/logs/05_endpoint_support.log"

echo "[6/6] Generating the four plot blocks and master figure..."
run_plots "$OUTDIR"

echo "[OK] Full SQANTI3 basic-statistics workflow completed."
echo "[OK] Canonical manifest: $BASIC_MANIFEST"
echo "[OK] Plot-ready inputs: $OUTDIR/Figure_input"
echo "[OK] Master figure: $OUTDIR/Figure2_MASTER/Figure2_ABCD_MASTER.pdf"
