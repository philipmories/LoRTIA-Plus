# KSHV benchmark companion dataset

## Overview

This archive is a **KSHV benchmark companion dataset** associated with the manuscript **“LoRTIA Plus: a chemistry-agnostic, feature-first software package for long-read transcriptome annotation”**. The KSHV component was designed to compare **LoRTIA Plus** with **bambu**, **FLAIR**, **IsoQuant**, and **NAGATA** in a compact, transcriptionally dense viral genome context using two ONT library types: **dcDNA** and **dRNA**. According to the manuscript, the viral benchmark evaluated **TSS**, **TES**, and **transcript-level** recovery separately, and the final KSHV reference elements were retained only if they were recovered by at least one annotator within the predefined positional tolerance. 

The Figshare version is intended to contain the **main KSHV benchmark and reference files**. The complete benchmarking workflow, downstream statistical summaries, and figure-generation / post-processing code are maintained separately in the project GitHub repository. The goal of this README is to describe both the KSHV data package structure and the role of the Python and R scripts that were used to generate and visualize the benchmark. fileciteturn53file19turn53file7

## Figshare directory structure

```text
KSHV/
├── README.md
├── annotator_gff3/
├── cell_active_reference/
├── references/
└── figure_input/
```

## Meaning of the directories

### `annotator_gff3/`

This directory contains the **annotator-specific GFF3 predictions** that served as the direct input of the viral benchmark. It should contain the KSHV output GFF3 files for **LoRTIA Plus**, **bambu**, **FLAIR**, **IsoQuant**, and **NAGATA**, separated by the two library types:

- **dcDNA**
- **dRNA**

A clear internal structure is for example:

```text
annotator_gff3/
├── dcDNA/
│   ├── LoRTIA.gff3
│   ├── bambu.gff3
│   ├── FLAIR.gff3
│   ├── IsoQuant.gff3
│   └── NAGATA.gff3
└── dRNA/
    ├── LoRTIA.gff3
    ├── bambu.gff3
    ├── FLAIR.gff3
    ├── IsoQuant.gff3
    └── NAGATA.gff3
```

The transcript-level benchmark was performed by comparing annotator GFF3 predictions against a reference TSV/GFF3 transcript set using strand-aware **TSS**, **TES**, and **intron-chain** windows. fileciteturn53file0turn53file2

### `cell_active_reference/`

This directory contains the **filtered / active benchmark reference files** that were actually retained for the final KSHV evaluation. The directory name follows the same general terminology used elsewhere in the project, but here it should be interpreted operationally: these are the reference subsets that entered the final benchmark for the corresponding KSHV dataset or library type.

Because the manuscript states that KSHV reference features were retained only if they were recovered by at least one annotator within the allowed positional tolerance, this directory is the appropriate place to store the final benchmark-restricted reference side for **dcDNA** and **dRNA**. fileciteturn53file19

A practical internal layout is for example:

```text
cell_active_reference/
├── dcDNA/
└── dRNA/
```

### `references/`

This directory contains the **core KSHV reference transcript and feature resources** that define the benchmark space. It should contain:

- the **KSHV reference transcript list** in TSV and/or GFF3 format,
- the reference **TSS**, **TES**, and **intron** feature files,
- and any merged reference files from which the benchmark was derived.

The uploaded scripts show that the benchmark could be built from merged transcript annotations together with separately written TSS/TES/intron GFF3 outputs. `merge_transcripts.py` can generate a merged transcript set together with separate TSS, TES, and intron files, while `compare_gff_tss_tes_introns.py` compares annotator-derived GFF/GFF3 feature files to a reference feature set.

### `figure_input/`

This directory contains the **plot-ready TSV files** used directly by `KSHV.R`. The script expects fixed input file names and specific required columns.

For the top metrics panels, the script expects:

- `dcDNA_stats_TSS.tsv`
- `dcDNA_stats_TES.tsv`
- `dcDNA_stats.tsv`
- `dRNA_stats_TSS.tsv`
- `dRNA_stats_TES.tsv`
- `dRNA_stats.tsv`

Each of these files must contain at least the columns:

- `TP`
- `FP`
- `FN`

The category panel requires:

- `category_dcdna.tsv`
- `category_dRNA.tsv`

with a `Category` column and at least one of:
- `isoform`
- `ref_id`

The UpSet section requires a combined false-positive TSV with:

- `Start`
- `End`
- `Strand`
- `Original_Transcript_ID` or `Cluster_ID`

and annotator-specific binary `*_present` columns for both **dRNA** and **dcDNA** false positives. fileciteturn53file9turn53file16

------------------------------------------

## KSHV benchmark workflow summary

```text
Multiple annotation GFF/GFF3/GTF files
    -> merge_transcripts.py
    -> merged transcript reference + TSS/TES/intron GFF3 files

Reference feature set + annotator-specific GFF3 files
    -> compare_gff_tss_tes_introns.py
    -> feature-level summary per reference element

Reference TSS/TES positions + comparison TSS/TES files
    -> compare_tsv.py
    -> position-level match counts

Reference transcript TSV + annotator GFF3 directory
    -> evaluate_transcripts_tsv_vs_gff3.py
    -> TP / FP / FN / precision / recall / F1 + matches + FP tables

Transcript TSV + category TSV
    -> Compare_category_with_windows.py
    -> matched transcript IDs and categories

Transcript TSV + numeric annotator-specific TSV
    -> compare_isoforms_with_windows.py
    -> numeric values transferred to matched transcript models

Multiple false-positive TSV files
    -> merge_tsv_transcripts_stats_by_windows.py
    -> combined FP cluster table with *_present columns

The resulting TSV files
    -> KSHV.R
    -> composite KSHV benchmark figure
```

------------------------------------------

## `merge_transcripts.py`

### What does it do?

This script reads transcript and exon annotations from multiple **GFF/GFF3/GTF** files and merges them into a unified transcript set based on **exon composition**. It can optionally apply a terminal window for fuzzy merging of nearly identical transcript models and preserves / aggregates sample-level count or score information. fileciteturn53file15turn53file11

### Example run command

```bash
python merge_transcripts.py \
  --indir input_gff/ \
  --outdir merged_reference/ \
  --pattern auto \
  --window 10 \
  --exon-score-mode max \
  --collapse-sites \
  --write-sample-attrs
```

### What input does it require?

- an input directory (`--indir`) searched recursively for **GFF/GFF3/GTF** files
- optional file-extension filtering (`--pattern`)
- transcript and exon annotations in the input files

### What is its operating logic?

- it loads transcript and exon data from each input file
- it merges identical transcript models by exon composition
- it can use transcript scores or exon-derived scores as count values
- it can generate separate **TSS**, **TES**, and **intron** outputs
- it can optionally collapse identical TSS/TES/intron sites across transcripts

### What is the output?

Using the output directory name as the base prefix, the script writes:

- `<outdir_name>.tsv`
- `<outdir_name>.gff3`
- `<outdir_name>.TSS.gff3`
- `<outdir_name>.TES.gff3`
- `<outdir_name>.introns.gff3` fileciteturn53file15

### Why is it useful in the KSHV benchmark?

It is suitable for building a unified **KSHV transcript reference** and the corresponding **TSS / TES / intron** reference feature sets. 

------------------------------------------

## `compare_gff_tss_tes_introns.py`

### What does it do?

This script compares **TSS**, **TES**, and **intron** features in one or more GFF/GFF3 files against a **reference GFF/GFF3 feature set** using a positional tolerance window. It sums feature scores per sample and writes both a **TSV** summary and a **merged GFF3** output. 

### Example run command

```bash
python compare_gff_tss_tes_introns.py \
  --main-gff references/kshv_reference_features.gff3 \
  --indir annotator_features/ \
  --out-tsv results/kshv_features.tsv \
  --out-gff3 results/kshv_features_merged.gff3 \
  --pos-window 10
```

### What input does it require?

- a reference GFF/GFF3 file (`--main-gff`)
- an input directory (`--indir`) containing annotator GFF/GFF3 files
- reference-side `TSS`, `TES`, and `intron` feature records

### What is its operating logic?

- it loads canonical reference features from the main GFF
- it scans all other input GFF/GFF3 files for the same feature types
- it matches input features to the reference set within the specified positional window
- it accumulates scores per reference feature and per sample

### What is the output?

- a feature-level **TSV**
- a merged **GFF3**

The merged GFF3 stores the reference features with aggregated scores. fileciteturn46file1

### Why is it useful in the KSHV benchmark?

It is useful for feature-level **TSS/TES/intron** benchmarking and for defining which reference-side features are recovered by the annotators.

------------------------------------------

## `compare_tsv.py`

### What does it do?

This script compares two **TSV or GFF/GFF3** files by **TSS or TES positions**, using strand-aware interpretation of `start` and `end`. For each row in the first input, it counts how many positions in the second input fall within the specified window and writes the result as an additional count column in the output TSV. fileciteturn46file3turn53file10

### Example run command

```bash
python compare_tsv.py \
  -i references/kshv_tss.tsv \
  -g annotator_tss.gff3 \
  -o results/kshv_tss_overlap.tsv \
  --mode TSS \
  --window 10
```

### What input does it require?

The first and second input can be:

- TSV
- GFF
- GFF3

For TSV input, the script requires at least:

```text
seqnames    strand    start    end
```

### What is its operating logic?

- it converts `start`, `end`, and `strand` into strand-aware **TSS** or **TES** positions
- it indexes the positions from the second file
- it counts how many positions fall within the specified window for each row of the first file

### What is the output?

- a **TSV** that preserves all columns from the first input
- plus one additional count column representing the number of positional matches. fileciteturn46file3

### Why is it useful in the KSHV benchmark?

It is useful for rapid reference-side TSS/TES matching and for determining how many predictions support a given reference site.

------------------------------------------

## `evaluate_transcripts_tsv_vs_gff3.py`

### What does it do?

This script compares a **reference TSV transcript set** against a directory of annotator **GFF3** predictions. It decides transcript-level matches using **strand + TSS window + TES window + intron-chain window**, and then calculates **TP**, **FP**, **FN**, **precision**, **recall**, and **F1** per GFF3 file. It also writes detailed match tables, false-positive tables, and a list of reference transcripts not recovered by any GFF3 prediction. fileciteturn53file2turn53file0

### Example run command

```bash
python evaluate_transcripts_tsv_vs_gff3.py \
  --ref references/KSHV_reference_transcripts.tsv \
  --gff3-dir annotator_gff3/dcDNA/ \
  --glob "*.gff3" \
  --tss-window 10 \
  --tes-window 10 \
  --intron-window 10 \
  --out results/dcDNA
```

### What input does it require?

#### Reference TSV
Required columns:

```text
isoform
chromosome
start
end
strand
exon_composition
reads
```

#### Prediction input
- a GFF3 directory (`--gff3-dir`)
- GFF3 files with `mRNA` and `exon` features

### What is its operating logic?

- it parses the reference TSV into transcript objects
- it parses each GFF3 file into transcript predictions
- it matches predictions to references using TSS, TES, and intron windows
- one reference may match multiple predictions, but still counts as one TP
- unmatched predictions are collected as false positives
- false positives can be clustered across files

### What is the output?

Using the `--out` prefix, the script writes:

- `<out>_matches.tsv`
- `<out>_stats.tsv`
- `<out>_false_positives.tsv`
- `<out>_reference_not_found.tsv` 

### Why is it useful in the KSHV benchmark?

This is the main **transcript-level benchmarking script**.

------------------------------------------

## `Compare_category_with_windows.py`

### What does it do?

This script compares two **TSV** files at the transcript level using **TSS**, **TES**, and **intron** windows. For each row in the base file, it identifies matching transcript models in the query file and appends their **transcript IDs** and **categories** as new output columns. 

### Example run command

```bash
python Compare_category_with_windows.py \
  --base transcript_table.tsv \
  --query category_table.tsv \
  --out category_annotated.tsv \
  --tss-window 10 \
  --tes-window 10 \
  --intron-window 5
```

### What input does it require?

#### Base TSV
Required columns:
- `chromosome` / `chrom` / `seqid`
- `start`
- `end`
- `strand`
- `exon_composition` / `exon_comp`

#### Query TSV
Required columns:
- `transcript_id` / `isoform` / `id`
- `chromosome` / `chrom` / `seqid`
- `start`
- `end`
- `strand`
- `exon_composition` / `exon_comp`
- `category`

### What is its operating logic?

- it searches for transcript-level matches between the base and query tables
- matching is decided using TSS/TES/intron boundary windows
- it records the matching query transcript IDs and categories

### What is the output?

- the original base TSV
- plus two new columns:
  - `Matched_Transcript_IDs`
  - `Matched_Categories`

(or the names specified with `--id-out-col` and `--cat-out-col`) 

### Why is it useful in the KSHV benchmark?

It is useful for generating the **category-annotated transcript table** used in the category-composition panel.

------------------------------------------

## `compare_isoforms_with_windows.py`

### What does it do?

This script compares two **TSV** files at the transcript level using **TSS**, **TES**, and **intron** windows. It preserves all columns from the base file and adds the **numeric value columns** from the query file by summing the values of all matching transcript models. 

### Example run command

```bash
python compare_isoforms_with_windows.py \
  --base transcript_table.tsv \
  --query quant_table.tsv \
  --out transcript_table_with_counts.tsv \
  --tss-window 10 \
  --tes-window 10 \
  --intron-window 5
```

### What input does it require?

Both TSV files must support transcript-level matching through the following fields:

```text
isoform
chromosome
start
end
strand
exon_composition
```

All query columns other than the fixed transcript-coordinate columns are treated as numeric values and are summed onto the base rows. fileciteturn46file2

### What is its operating logic?

- it finds transcript-level matches using TSS/TES/intron windows
- it sums all numeric query columns for matching transcript models
- it writes those values back onto the base table

### What is the output?

- a **TSV** containing all base columns
- plus the transferred numeric query columns. 

### Why is it useful in the KSHV benchmark?

It is useful for transferring matched quantitative / support values onto a common transcript table.

------------------------------------------

## `merge_tsv_transcripts_stats_by_windows.py`

### What does it do?

This script merges multiple **TSV transcript lists** using separate **TSS**, **TES**, and **intron** windows to decide matches. It builds common **Cluster_ID** values and preserves prefixed `*_present` columns from the input files, producing one row per merged cluster. 

### Example run command

```bash
python merge_tsv_transcripts_stats_by_windows.py \
  -i fp_loRTIA.tsv fp_bambu.tsv fp_flair.tsv fp_isoquant.tsv fp_nagata.tsv \
  -o FP.tsv \
  -t 10 \
  -e 10 \
  -n 5 \
  -x
```

### What input does it require?

Multiple input TSV files (`-i / --inputs`) containing at least:

- `Transcript_ID` / `transcript_id` / `ID`
- `chromosome` / `chrom` / `seqid`
- `Start`
- `End`
- `Strand`
- optionally `Exon_Composition`
- optionally pre-existing `*_present` columns

### What is its operating logic?

- it loads transcript rows from each input TSV
- it searches for pairwise matches using TSS/TES/intron windows
- it clusters matched transcripts using a union-find approach
- it assigns new shared `Cluster_ID` values
- it preserves the source-specific `*_present` columns with filename-based prefixes

### What is the output?

- a **cluster-level TSV**, for example `FP.tsv`, containing fields such as:
  - `Cluster_ID`
  - `Original_Transcript_ID`
  - `chrom`
  - `Start`
  - `End`
  - `Strand`
  - `Exon_Composition`
  - `Introns`
  - source-prefixed annotator-specific `*_present` columns. fileciteturn46file6

### Why is it useful in the KSHV benchmark?

It is used to create the combined **false-positive cluster table** required for the UpSet section of the KSHV figure. The `*_present` output pattern is consistent with the `dRNA_false_positives__..._present` and `dcDNA_false_positives__..._present` columns expected by `KSHV.R`. 

------------------------------------------

## `KSHV.R`

### What does it do?

This R script generates a **composite KSHV benchmark figure** with three major components:

- **Panel A:** TSS / TES / transcript precision–recall–F1 panels
- **Panel B:** transcript category stacked bar plot
- **Panel C:** UpSet plots based on one combined false-positive TSV. 

### Example run command

```bash
Rscript KSHV.R
```

### What input does it require?

#### Metrics panels
The script looks for the following fixed files:

- `dcDNA_stats_TSS.tsv`
- `dcDNA_stats_TES.tsv`
- `dcDNA_stats.tsv`
- `dRNA_stats_TSS.tsv`
- `dRNA_stats_TES.tsv`
- `dRNA_stats.tsv`

Required columns:
- `TP`
- `FP`
- `FN` fileciteturn53file9

#### Category panel
The script expects:

- `category_dcdna.tsv`
- `category_dRNA.tsv`

Required columns:
- `Category`
- and at least one of `isoform` or `ref_id`

If `.gff3_match` or `_match` annotator columns are present, the script converts them into annotator-specific category-composition summaries. 

#### UpSet panel
The script expects one combined false-positive TSV.

Required columns:
- `Start`
- `End`
- `Strand`
- `Original_Transcript_ID` or `Cluster_ID`

It also requires binary annotator-specific columns matching:

- `dRNA_false_positives__<annotator>_present`
- `dcDNA_false_positives__<annotator>_present` 

### What is its operating logic?

- it reads the metrics tables and computes precision / recall / F1 panels for the five annotators
- it reads the category TSV files and creates dcDNA and dRNA stacked category-composition panels
- it reads the combined false-positive TSV and builds annotator-overlap UpSet plots

### What is the output?

The script assembles the metrics, category, and UpSet components into a combined KSHV benchmark figure and writes PDF outputs for the figure panels and the combined figure. The script structure shows that this is the main KSHV summary figure used for the benchmark visualization. 

### What does the figure show?

- how accurately each annotator recovers **TSS**, **TES**, and full **transcript models**
- how the recovered matched transcripts are distributed across broad transcript categories
- how shared or annotator-specific the **false-positive transcript models** are

------------------------------------------

## Short summary

In this KSHV benchmark workflow:

- `merge_transcripts.py` generates merged transcript and feature reference files
- `compare_gff_tss_tes_introns.py` performs feature-level TSS/TES/intron comparison against a reference GFF
- `compare_tsv.py` performs positional TSS or TES matching between TSV/GFF inputs
- `evaluate_transcripts_tsv_vs_gff3.py` is the main transcript-level benchmarking script
- `Compare_category_with_windows.py` adds matched transcript categories
- `compare_isoforms_with_windows.py` transfers numeric values onto matched transcript models
- `merge_tsv_transcripts_stats_by_windows.py` builds the combined false-positive cluster table
- `KSHV.R` assembles the final composite KSHV benchmark figure.
