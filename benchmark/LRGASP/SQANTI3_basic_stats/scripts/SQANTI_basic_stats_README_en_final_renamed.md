# SQANTI3 basic-statistics workflow for the six-annotator LRGASP benchmark

## 1. Purpose

This package contains the complete SQANTI3 basic-statistics and plotting workflow used to compare six transcript annotators across human LRGASP long-read RNA-seq datasets.

The complete experimental design contains:

- **6 annotators:** LoRTIA-Plus, bambu, FLAIR, IsoQuant, NAGATA, StringTie2
- **5 chemistries:** ONT-CapTrap, ONT-cDNA, ONT-dRNA, PacBio, PacBio-CapTrap
- **3 cell lines:** H1, H1-endo, WTC11
- **90 SQANTI3 `classification.txt` files** in total

The workflow is intentionally separated into modules so that users can either reproduce the complete analysis from SQANTI3 classification files or regenerate the figure directly from the archived plot-ready TSV files.

---

## 2. Directory structure

```text
SQANTI_basic_stats_Figshare_FINAL/
├── README.md
├── requirements.txt
├── R_PACKAGES.txt
├── SCRIPT_LIST.tsv
├── run_sqanti_basic_stats_6annotators.sh
│
├── manifest/
│   ├── README.md
│   └── manifest_reference.tsv
│
├── input/
│   ├── README.md
│   ├── SQANTI3_basic_stat.tsv
│   ├── sqanti_length_boxplot_stats.tsv
│   └── isoformsgene.tsv
│
├── scripts/
│   ├── preprocessing/
│   │   ├── build_basic_stats_manifest.py
│   │   ├── build_sqanti_basic_stat_table.py
│   │   ├── extract_sqanti_length_and_category_summaries.py
│   │   ├── summarize_isoforms_per_gene.py
│   │   └── summarize_sqanti_endpoint_support.py
│   │
│   └── plotting/
│       ├── run_sqanti_basic_stats_all.R
│       ├── plot_sqanti_yield_and_category_composition.R
│       ├── plot_sqanti_isoform_length_distributions.R
│       ├── plot_sqanti_coding_noncoding_composition.R
│       ├── plot_sqanti_isoforms_per_gene.R
│       └── assemble_sqanti_basic_stat_master_figure.R
│
└── output/
    └── README.md
```

The raw SQANTI3 `classification.txt` files are not duplicated in this package. They may be archived separately in the main SQANTI3 data tree. The source manifest records which classification file belongs to each annotator × chemistry × cell-line combination.

---

## 3. Workflow overview

```text
                    FULL REGENERATION

90 SQANTI3 classification.txt files
                │
                ▼
manifest/manifest_reference.tsv
                │
                ▼
MODULE 1 — Manifest preparation
build_basic_stats_manifest.py
                │
                ▼
canonical 90-row manifest
                │
     ┌──────────┼───────────────┐
     ▼          ▼               ▼
MODULE 2A   MODULE 2B       MODULE 2C
basic stat  length/category isoforms/gene
     │          │               │
     ▼          ▼               ▼
SQANTI3_     sqanti_length_   isoformsgene.tsv
basic_stat   boxplot_stats
.tsv         .tsv
     └──────────┼───────────────┘
                ▼
output/Figure_input/
                │
                ▼
MODULE 3 — R plotting driver
run_sqanti_basic_stats_all.R
                │
     ┌──────────┼──────────┬──────────┐
     ▼          ▼          ▼          ▼
   Block A    Block B    Block C    Block D
     └──────────┼──────────┴──────────┘
                ▼
MODULE 4 — Master assembly
Figure2_ABCD_MASTER.pdf
```

An additional auxiliary module, `summarize_sqanti_endpoint_support.py`, generates detailed endpoint-support statistics. It is not an input to the A-D master figure.

---

## 4. Dependencies

### Python

Python 3 with:

```text
pandas
numpy
```

Install with:

```bash
python -m pip install -r requirements.txt
```

### R

Required R packages:

```text
readr
dplyr
tidyr
ggplot2
stringr
scales
patchwork
ggh4x
```

For example:

```r
install.packages(c(
  "readr", "dplyr", "tidyr", "ggplot2",
  "stringr", "scales", "patchwork", "ggh4x"
))
```

`grid` is also used by the plotting scripts and is part of the standard R distribution.

---

# PART I — Reproduce the archived figure only

## 5. Fastest reproducible route: plot-only mode

The three files already stored under `input/` are the exact plot-ready tables needed by the R plotting workflow:

```text
input/SQANTI3_basic_stat.tsv
input/sqanti_length_boxplot_stats.tsv
input/isoformsgene.tsv
```

Therefore, users who only want to regenerate the figure do **not** need access to the 90 original SQANTI3 classification files.

From the package root:

```bash
chmod +x run_sqanti_basic_stats_6annotators.sh
./run_sqanti_basic_stats_6annotators.sh --plot-only
```

This performs the following operations:

1. creates `output/Figure_input/`;
2. copies the three archived TSV files into that directory;
3. calls `scripts/plotting/run_sqanti_basic_stats_all.R`;
4. generates the four plotting blocks;
5. combines them into the master figure.

The principal final output is:

```text
output/Figure2_MASTER/Figure2_ABCD_MASTER.pdf
```

The `Figure2` label is retained as a legacy code/output filename and does not necessarily correspond to the final manuscript figure numbering.

To use another output directory:

```bash
./run_sqanti_basic_stats_6annotators.sh --plot-only input results/plot_only
```

---

# PART II — Full regeneration from SQANTI3 classification files

## 6. Source manifest

The full workflow starts from:

```text
manifest/manifest_reference.tsv
```

It contains the complete **90-row** six-annotator design with columns:

```text
Chemistry    Cell-line    GTF
```

### Legacy naming of the `GTF` column

In this reused manifest, `GTF` does **not** contain transcript GTF files. It contains the paths to the corresponding SQANTI3 `classification.txt` files. The column name is retained because the manifest originated from an earlier transcript-recovery workflow.

### Absolute paths

The archived manifest preserves the original machine-specific absolute paths used for the analysis. For a full rerun on another computer, make a working copy and replace those paths with the locations of the corresponding archived SQANTI3 classification files.

The archived file itself should be retained unchanged for provenance.

---

## 7. One-command full workflow

If all `GTF` entries in the manifest resolve to the corresponding SQANTI3 classification files, run:

```bash
chmod +x run_sqanti_basic_stats_6annotators.sh
./run_sqanti_basic_stats_6annotators.sh \
  manifest/manifest_reference.tsv \
  output
```

The first argument is the source manifest. The second is the output directory.

The defaults are already `manifest/manifest_reference.tsv` and `output`, so the same complete workflow can also be started with:

```bash
./run_sqanti_basic_stats_6annotators.sh
```

The top-level shell script is the workflow orchestrator. It executes Modules 1-4 in sequence and writes one log file for each major step.

---

# MODULE 1 — Manifest preparation

## 8. `build_basic_stats_manifest.py`

Purpose: converts the source manifest into the canonical manifest used by all Python statistics modules and validates the complete factorial design.

Input:

```text
manifest/manifest_reference.tsv
```

Required input columns:

```text
Chemistry
Cell-line
GTF
```

Output:

```text
output/manifests/manifest_basic_stats_6annotators.tsv
```

Canonical output columns:

```text
Annotator    Chemical    Cell    Classification
```

For the complete analysis the output must contain:

```text
6 annotators × 5 chemistries × 3 cell lines = 90 rows
```

The manifest-based design is important because it guarantees that the same defined set of SQANTI3 files is propagated through all downstream summary modules rather than discovering files recursively from arbitrary directories.

---

# MODULE 2 — SQANTI3 summary statistics

## 9. Module 2A: basic SQANTI3 statistics

Script:

```text
scripts/preprocessing/build_sqanti_basic_stat_table.py
```

Input:

```text
output/manifests/manifest_basic_stats_6annotators.tsv
```

Principal output:

```text
output/Figure_input/SQANTI3_basic_stat.tsv
```

This table contains summaries by chemistry × cell line × SQANTI3 structural category × annotator, including isoform counts and support-related rates. It is used by two plotting branches:

- Block A — isoform yield and SQANTI3 class composition
- Block C — coding/non-coding composition

---

## 10. Module 2B: transcript-length and category summaries

Script:

```text
scripts/preprocessing/extract_sqanti_length_and_category_summaries.py
```

The workflow calls this script in **manifest mode** so that exactly the same set of classification files is used as in the other modules.

Outputs:

```text
output/category_and_length_outputs/sqanti_category_summary.tsv
output/category_and_length_outputs/sqanti_length_boxplot_stats.tsv
```

The plot-ready length table is also copied to:

```text
output/Figure_input/sqanti_length_boxplot_stats.tsv
```

This is the input to Block B — transcript-length distributions.

---

## 11. Module 2C: isoforms per gene

Script:

```text
scripts/preprocessing/summarize_isoforms_per_gene.py
```

Outputs:

```text
output/isoforms_per_gene_outputs/isoformsgene.tsv
output/isoforms_per_gene_outputs/sqanti_category_stats.tsv
```

The main table is also copied to:

```text
output/Figure_input/isoformsgene.tsv
```

The gene-level isoform bins are:

```text
1
2-3
4-5
>=6
```

This table is the input to Block D — number of isoforms per gene.

---

## 12. Auxiliary endpoint-support module

Script:

```text
scripts/preprocessing/summarize_sqanti_endpoint_support.py
```

Output:

```text
output/endpoint_support/sqanti_endpoint_support.tsv
```

This module summarizes detailed endpoint-related support statistics by annotator × chemistry × cell line and structural class. It includes CAGE/poly(A)-related and TSS/TES support summaries.

It is archived as part of the basic-statistics analysis but is **not required to construct the A-D master figure**.

---

# MODULE 3 — Plotting

## 13. R plotting driver

The R plotting sub-workflow is controlled by:

```text
scripts/plotting/run_sqanti_basic_stats_all.R
```

It expects exactly these three files under the selected output directory:

```text
Figure_input/SQANTI3_basic_stat.tsv
Figure_input/sqanti_length_boxplot_stats.tsv
Figure_input/isoformsgene.tsv
```

It then executes the four plot scripts and the master-assembly script in the same R session.

This same-session execution is important because the master assembly does not read the TSV files directly. Instead, it combines four R plot objects generated by the preceding modules.

---

## 14. Block A — isoform yield and SQANTI3 class composition

Script:

```text
plot_sqanti_yield_and_category_composition.R
```

Input:

```text
SQANTI3_basic_stat.tsv
```

Main R object:

```text
Ablock_MAIN
```

The plot summarizes reconstructed isoform yield and the relative composition of SQANTI3 structural classes across annotators and long-read chemistries.

---

## 15. Block B — transcript-length distributions

Script:

```text
plot_sqanti_isoform_length_distributions.R
```

Input:

```text
sqanti_length_boxplot_stats.tsv
```

Main R object:

```text
IsoLen_MAIN
```

The plot compares transcript-length distributions across annotators and chemistries.

---

## 16. Block C — coding/non-coding composition

Script:

```text
plot_sqanti_coding_noncoding_composition.R
```

Input:

```text
SQANTI3_basic_stat.tsv
```

Main R object:

```text
Coding_MAIN
```

The plot compares coding and non-coding proportions for the main SQANTI3 transcript classes.

---

## 17. Block D — isoforms per gene

Script:

```text
plot_sqanti_isoforms_per_gene.R
```

Input:

```text
isoformsgene.tsv
```

Main R object:

```text
IsoformsPerGene_MAIN
```

The plot summarizes gene-level reconstruction complexity using the bins 1, 2-3, 4-5, and >=6 isoforms per gene.

---

# MODULE 4 — Master figure assembly

## 18. `assemble_sqanti_basic_stat_master_figure.R`

This script combines the four objects created by the plotting scripts:

```text
Ablock_MAIN
IsoLen_MAIN
Coding_MAIN
IsoformsPerGene_MAIN
```

They are stacked vertically and labeled A-D.

Principal output:

```text
output/Figure2_MASTER/Figure2_ABCD_MASTER.pdf
```

Because the master script depends on objects held in the current R session, it should normally be invoked through `run_sqanti_basic_stats_all.R` rather than run independently.

---

## 19. Output structure after a complete run

```text
output/
├── Figure_input/
│   ├── SQANTI3_basic_stat.tsv
│   ├── sqanti_length_boxplot_stats.tsv
│   └── isoformsgene.tsv
│
├── category_and_length_outputs/
│   ├── sqanti_category_summary.tsv
│   └── sqanti_length_boxplot_stats.tsv
│
├── isoforms_per_gene_outputs/
│   ├── isoformsgene.tsv
│   └── sqanti_category_stats.tsv
│
├── endpoint_support/
│   └── sqanti_endpoint_support.tsv
│
├── manifests/
│   └── manifest_basic_stats_6annotators.tsv
│
├── logs/
│   ├── 01_build_manifest.log
│   ├── 02_basic_stat.log
│   ├── 03_length_category.log
│   ├── 04_isoforms_per_gene.log
│   ├── 05_endpoint_support.log
│   └── plotting.log
│
├── Ablock_IsoformYield_SQANTI3/
├── Bblock_IsoLen_Boxplots/
├── coding_composition_plots_Astyle_header/
├── IsoformsPerGene_Astyle_header/
│
└── Figure2_MASTER/
    └── Figure2_ABCD_MASTER.pdf
```

The individual R scripts may also generate additional detailed or cell-line-resolved PDFs according to the settings stored in those scripts.

---

## 20. Which route should a user choose?

### To reproduce the published/archived figure

Use:

```bash
./run_sqanti_basic_stats_6annotators.sh --plot-only
```

This route is self-contained within this package because the three required Figure-input TSV files are archived under `input/`.

### To regenerate the statistics from SQANTI3 outputs

Use:

```bash
./run_sqanti_basic_stats_6annotators.sh manifest/manifest_reference.tsv output
```

This route additionally requires access to all 90 SQANTI3 `classification.txt` files referenced by the manifest. If the archive is moved to another machine, update a working copy of the manifest paths first.

---

## 21. Provenance note about the top-level runner

The analysis documentation used during the six-annotator revision referred to a one-command shell runner named `run_sqanti_basic_stats_6annotators.sh`. The individual Python analysis modules, the R plotting driver, and the R plotting scripts are archived in this package. The top-level shell runner supplied here reconstructs that documented orchestration using those exact modules and provides an additional `--plot-only` route for transparent reproduction from the archived derived tables.

No statistical calculations are implemented inside the shell runner itself; it only connects the archived analysis modules in the documented order.
