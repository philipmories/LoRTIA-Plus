# SQANTI-based transcript recovery benchmark

## Overview

This document describes the SQANTI-based transcript recovery benchmark workflow. The benchmark evaluates recovery of known, catalog-supported transcript structures across annotators, sequencing chemistries, and cell lines using SQANTI3 `classification.txt` outputs and chemistry × cell line-specific active reference GTF files derived from GENCODE v48.

The workflow consists of three main steps:

1. **Build chemistry × cell line-specific active reference transcript sets**
2. **Run the main transcript recovery benchmark**
3. **Summarize novelty burden, junction support, and plausibility metrics**

The benchmark produces not only summary tables, but also **plot-ready TSV files** that can be read directly by the accompanying R visualization scripts.

## Main file names

### Python scripts

- `build_sqanti_cell_active_references.py`
- `benchmark_sqanti_transcript_recovery.py`
- `summarize_sqanti_novelty_support.py`

### R scripts

- `FSM-ISM-visualization.R`
- `NNC-NIC-visualization.R`

## Figshare structure

The SQANTI benchmark files are archived on Figshare.

**Figshare DOI:** `INSERT_FIGSHARE_DOI_HERE`

The benchmark files are organized under the `SQANTI3/` directory.

```text
SQANTI3/
├── BAMBU/
│   └── <Chemistry>/<Cell-line>/
│       ├── *_classification.txt
│       ├── *_junctions.txt
│       └── *.gtf
├── FLAIR/
│   └── <Chemistry>/<Cell-line>/
│       ├── *_classification.txt
│       ├── *_junctions.txt
│       └── *.gtf
├── IsoQuant/
│   └── <Chemistry>/<Cell-line>/
│       ├── *_classification.txt
│       ├── *_junctions.txt
│       └── *.gtf
├── LoRTIA/
│   └── <Chemistry>/<Cell-line>/
│       ├── *_classification.txt
│       ├── *_junctions.txt
│       └── *.gtf
├── NAGATA/
│   └── <Chemistry>/<Cell-line>/
│       ├── *_classification.txt
│       ├── *_junctions.txt
│       └── *.gtf
├── cell_active_transcript_reference/
│   └── <Chemistry>-<Cell-line>-cell_active.gtf
└── Figure_input/
    ├── sqanti_fsmism_lortia_vs_others.tsv
    ├── SQANTI_Recall_FSM_ISM.tsv
    ├── SQANTI_FL_ratio.tsv
    ├── nic_nnc_burden_by_condition.tsv
    ├── novel_junction_support_by_annotator_chemistry.tsv
    ├── novel_junction_plausibility_by_annotator_chemistry.tsv
    ├── novel_isoform_canonicality_by_annotator_chemistry_category.tsv
    └── novel_locus_consistency_by_annotator_chemistry_category.tsv
```

### Meaning of the subdirectories

- **`<Annotator>/<Chemistry>/<Cell-line>/`**  
  Annotator-specific SQANTI3 files organized by chemistry and cell line:
  - `classification.txt`
  - `junctions.txt`
  - the annotator GTF file

- **`cell_active_transcript_reference/`**  
  Chemistry × cell line-specific active reference GTF files. These are derived from GENCODE v48 using FSM/ISM-supported `associated_transcript` identifiers from the corresponding chemistry × cell line group.

- **`Figure_input/`**  
  Plot-ready input tables used by the R visualization scripts.

## Workflow summary

```text
SQANTI3 classification.txt files
    -> build_sqanti_cell_active_references.py
GENCODE v48 GTF + chemistry × cell line grouping
    -> cell_active.gtf reference files

SQANTI3 classification.txt + cell_active.gtf
    -> benchmark_sqanti_transcript_recovery.py
    -> recovery metrics + ANOVA + LoRTIA-vs-others tables + plot-ready TSVs

SQANTI3 classification.txt + SQANTI3 junctions.txt
    -> summarize_sqanti_novelty_support.py
    -> NIC/NNC burden + junction support + plausibility + locus-level consistency
```

------------------------------------------

## `build_sqanti_cell_active_references.py`

### What does it do?

This script builds chemistry × cell line-specific `cell_active.gtf` reference files.  
For each chemistry × cell line group, it reads all SQANTI3 `classification.txt` files listed in the manifest, keeps transcripts classified as **FSM** or **ISM**, collects their `associated_transcript` identifiers, and subsets a GENCODE v48 GTF accordingly. The result is a sample-aware reference transcript set used by the downstream benchmark.

### Example run command

```bash
python build_sqanti_cell_active_references.py   -m manifests/manifest_reference.tsv   -g reference/gencode.v48.annotation.gtf.gz   -o SQANTI3/cell_active_transcript_reference/   --keep-structure   --strict   -v
```

### What input does it require?

1. A **3-column manifest**
2. A **GENCODE v48 GTF** file

### What does the manifest look like?

Required columns:

```text
Chemistry    Cell-line    GTF
```

The `GTF` column contains the path to a SQANTI3 `classification.txt` file.

### Example manifest

```text
Chemistry    Cell-line    GTF
ONT-CapTrap  H1           /SQANTI3/BAMBU/ONT-CapTrap/H1/H1-ONT-CapTrap-BAMBU_classification.txt
ONT-CapTrap  H1           /SQANTI3/FLAIR/ONT-CapTrap/H1/H1-ONT-CapTrap-FLAIR_classification.txt
ONT-CapTrap  H1           /SQANTI3/IsoQuant/ONT-CapTrap/H1/H1-ONT-CapTrap-IsoQuant_classification.txt
ONT-CapTrap  H1           /SQANTI3/LoRTIA/ONT-CapTrap/H1/H1-ONT-CapTrap-LoRTIA_classification.txt
ONT-CapTrap  H1           /SQANTI3/NAGATA/ONT-CapTrap/H1/H1-ONT-CapTrap-NAGATA_classification.txt
```

### Which columns are required in the classification input?

The script uses:

- `structural_category`
- `associated_transcript`

Only FSM/ISM-classified rows are used to define the active reference transcript set.

### What is the output?

The script writes one chemistry × cell line-specific GTF file per group:

```text
<Chemistry>-<Cell-line>-cell_active.gtf
```

Examples:

```text
ONT-CapTrap-H1-cell_active.gtf
ONT-CapTrap-H1-endo-cell_active.gtf
PacBio-WTC11-cell_active.gtf
```

These files are stored in:

```text
SQANTI3/cell_active_transcript_reference/
```

The output GTF files contain the subset of GENCODE v48 transcripts supported by FSM/ISM evidence in the corresponding chemistry × cell line context.

------------------------------------------

## `benchmark_sqanti_transcript_recovery.py`

### What does it do?

This script runs the main SQANTI-based transcript recovery benchmark.  
It compares each annotator’s SQANTI3 `classification.txt` output against the matching chemistry × cell line `cell_active.gtf` reference file and computes transcript recovery metrics based on FSM, ISM, and FSM∪ISM transcript sets. It also performs chemistry-stratified blocked ANOVA and planned LoRTIA-vs-other-annotator paired comparisons, and writes plot-ready TSV files for downstream visualization.

### Example run command

```bash
python benchmark_sqanti_transcript_recovery.py   -m manifests/manifest_reference_transcript.tsv   -o output/sqanti_benchmark/   --lortia-name LoRTIA   --strict   --write-composition   -v
```

### What input does it require?

A **4-column manifest** containing:

```text
Chemistry    Cell-line    GTF    Reference
```

- `GTF`: path to SQANTI3 `classification.txt`
- `Reference`: path to the corresponding `cell_active.gtf` file

### Example manifest

```text
Chemistry    Cell-line    GTF                                                                 Reference
ONT-CapTrap  H1           /SQANTI3/BAMBU/ONT-CapTrap/H1/H1-ONT-CapTrap-BAMBU_classification.txt      /SQANTI3/cell_active_transcript_reference/ONT-CapTrap-H1-cell_active.gtf
ONT-CapTrap  H1           /SQANTI3/FLAIR/ONT-CapTrap/H1/H1-ONT-CapTrap-FLAIR_classification.txt      /SQANTI3/cell_active_transcript_reference/ONT-CapTrap-H1-cell_active.gtf
ONT-CapTrap  H1           /SQANTI3/IsoQuant/ONT-CapTrap/H1/H1-ONT-CapTrap-IsoQuant_classification.txt /SQANTI3/cell_active_transcript_reference/ONT-CapTrap-H1-cell_active.gtf
ONT-CapTrap  H1           /SQANTI3/LoRTIA/ONT-CapTrap/H1/H1-ONT-CapTrap-LoRTIA_classification.txt    /SQANTI3/cell_active_transcript_reference/ONT-CapTrap-H1-cell_active.gtf
ONT-CapTrap  H1           /SQANTI3/NAGATA/ONT-CapTrap/H1/H1-ONT-CapTrap-NAGATA_classification.txt    /SQANTI3/cell_active_transcript_reference/ONT-CapTrap-H1-cell_active.gtf
```

### Which columns are required in the classification input?

The script uses:

- `structural_category`
- `associated_transcript`

From these it defines:
- the FSM set
- the ISM set
- the FSM∪ISM set

and intersects them with the reference transcript universe defined by `cell_active.gtf`.

### What metrics does it calculate?

The main benchmark metrics are:

- `Recall_FSM`
- `Recall_ISM`
- `Recall_FSM_ISM`
- `FL_ratio_FSM_over_FSM_ISM`
- `ISM_share`

In brief:
- `Recall_FSM`: fraction of active reference transcripts recovered as FSM
- `Recall_ISM`: fraction of active reference transcripts recovered as ISM
- `Recall_FSM_ISM`: fraction of active reference transcripts recovered as FSM or ISM
- `FL_ratio_FSM_over_FSM_ISM`: FSM proportion among recovered catalog transcripts
- `ISM_share`: ISM proportion among recovered catalog transcripts

### What statistical analyses does it perform?

For each chemistry, the script performs blocked ANOVA:

```text
value ~ C(annotator) + C(cell)
```

It also runs paired LoRTIA-vs-others comparisons across cell lines with BH correction and Cohen’s `dz` effect size.

### What is the output?

The script writes:

**Main benchmark outputs**
- `sqanti_fsmism_metrics.tsv`
- `sqanti_fsmism_anova.tsv`
- `sqanti_fsmism_lortia_vs_others.tsv`

**Plot-ready outputs**
- `SQANTI_Recall_FSM_ISM.tsv`
- `SQANTI_FL_ratio.tsv`
- `SQANTI_Recall_FSM_ISM__pairwise.tsv`
- `SQANTI_FL_ratio__pairwise.tsv`

**Optional output**
- `sqanti_structural_category_counts.tsv`  
  written only when `--write-composition` is used.

### What information do these files contain?

- `sqanti_fsmism_metrics.tsv`: per-sample benchmark metrics by annotator
- `sqanti_fsmism_anova.tsv`: chemistry-specific ANOVA results
- `sqanti_fsmism_lortia_vs_others.tsv`: planned pairwise LoRTIA-vs-other comparisons
- `SQANTI_Recall_FSM_ISM.tsv`: boxplot/scatter-compatible input for combined FSM+ISM recovery
- `SQANTI_FL_ratio.tsv`: boxplot/scatter-compatible input for full-lengthness
- `SQANTI_Recall_FSM_ISM__pairwise.tsv`: pairwise heatmap input for recovery comparisons
- `SQANTI_FL_ratio__pairwise.tsv`: pairwise heatmap input for full-lengthness comparisons

------------------------------------------

## `summarize_sqanti_novelty_support.py`

### What does it do?

This script summarizes the burden, support, and plausibility of novel transcript structures.  
It is separate from the main FSM/ISM recovery benchmark and performs a dedicated novelty/plausibility analysis based on SQANTI3 classification and junction outputs. The script produces NIC/NNC burden, novel junction support, splice-site plausibility, novel isoform canonicality, and locus-level consistency tables.

### Example run command

```bash
python summarize_sqanti_novelty_support.py   --classification-manifest manifests/sqanti_classification_manifest.tsv   --junction-manifest manifests/sqanti_junction_manifest.tsv   --outdir output/sqanti_novelty/
```

### What input does it require?

Two separate manifests:

1. **classification manifest**
2. **junction manifest**

### What does the classification manifest look like?

Required columns:

```text
Chemistry    Cell-line    GTF
```

Optional column:
- `Annotator`

If `Annotator` is not present, the script attempts to infer it from the file path.

### Example classification manifest

```text
Chemistry    Cell-line    GTF                                                      Annotator
ONT-dRNA     H1           /SQANTI3/LoRTIA/ONT-dRNA/H1/H1-ONT-dRNA-LoRTIA_classification.txt     LoRTIA
ONT-dRNA     H1           /SQANTI3/IsoQuant/ONT-dRNA/H1/H1-ONT-dRNA-IsoQuant_classification.txt   IsoQuant
PacBio       WTC11        /SQANTI3/FLAIR/PacBio/WTC11/WTC11-PacBio-FLAIR_classification.txt       FLAIR
```

### What does the junction manifest look like?

Required columns:

```text
Annotator    Chemical    Cell    Junctions
```

Important: this script expects the exact column names `Chemical` and `Cell` for the junction manifest.

### Example junction manifest

```text
Annotator    Chemical     Cell     Junctions
LoRTIA       ONT-dRNA     H1       /SQANTI3/LoRTIA/ONT-dRNA/H1/H1-ONT-dRNA-LoRTIA_junctions.txt
IsoQuant     ONT-dRNA     H1       /SQANTI3/IsoQuant/ONT-dRNA/H1/H1-ONT-dRNA-IsoQuant_junctions.txt
FLAIR        PacBio       WTC11    /SQANTI3/FLAIR/PacBio/WTC11/WTC11-PacBio-FLAIR_junctions.txt
```

### Which columns are required in the classification input?

The script uses:

- `isoform`
- `structural_category`
- `associated_gene`
- `associated_transcript`
- `length`
- `exons`
- `all_canonical`
- `predicted_NMD`
- `bite`

### Which columns are required in the junction input?

The script uses:

- `isoform`
- `chrom`
- `strand`
- `genomic_start_coord`
- `genomic_end_coord`
- `junction_category`
- `splice_site`
- `RTS_junction`
- `indel_near_junct`

### What is the output?

The script writes:

- `nic_nnc_burden_by_condition.tsv`
- `nic_nnc_burden_by_annotator.tsv`
- `novel_junction_support_by_annotator.tsv`
- `novel_junction_support_by_annotator_chemistry.tsv`
- `novel_junction_plausibility_by_annotator_chemistry.tsv`
- `novel_isoform_canonicality_by_annotator_chemistry_category.tsv`
- `novel_locus_consistency_by_annotator_chemistry_category.tsv`

### What information do these files contain?

- `nic_nnc_burden_*`: NIC/NNC and total novel isoform burden
- `novel_junction_support_*`: support structure of novel junctions across annotators, chemistries, and cell lines
- `novel_junction_plausibility_*`: splice-site canonicality, RTS, and indel-near-junction summaries
- `novel_isoform_canonicality_*`: isoform-level canonicality, predicted NMD, BITE, length, and exon-count summaries
- `novel_locus_consistency_*`: gene-level accumulation and support of novel transcript structures

------------------------------------------

## Plot-ready benchmark outputs

The workflow produces not only summary tables but also plot-ready TSV files. The transcript recovery benchmark writes recovery and pairwise comparison tables for FSM/ISM-based performance, while the novelty-support workflow writes burden, support, plausibility, and locus-consistency tables for NIC/NNC-focused visualization.

------------------------------------------

## `FSM-ISM-visualization.R`

### What does it do?

This R script generates a **combined SQANTI transcript recovery figure** composed of three panels:

- **Panel A:** forest plot of LoRTIA-vs-other-annotator effect sizes for **FSM+ISM recovery**
- **Panel B:** forest plot of LoRTIA-vs-other-annotator effect sizes for **full-lengthness**, defined as `FSM / (FSM + ISM)`
- **Panel C:** scatter plot showing the **trade-off between recovery and full-lengthness** across annotators and chemistries

The script combines these panels into a single PDF figure.

### Example run command

```bash
Rscript FSM-ISM-visualization.R
```

### Which files does it require?

This script uses files written by `benchmark_sqanti_transcript_recovery.py`:

- `sqanti_fsmism_lortia_vs_others.tsv`
- `SQANTI_Recall_FSM_ISM.tsv`
- `SQANTI_FL_ratio.tsv`

### What do the input files need to contain?

#### `sqanti_fsmism_lortia_vs_others.tsv`
Used for the two forest plots. The script expects at least:

```text
Metric    Chemistry    Annotator_2    mean_delta    sd_delta    n_cells    p_adj_BH
```

These values are used to calculate:
- effect size (`mean_delta`)
- standard error and 95% confidence interval
- significance marker (`p_adj_BH < 0.05`)

#### `SQANTI_Recall_FSM_ISM.tsv`
Used for the recovery axis of the scatter plot. It must contain at least:

```text
Chemistry    Cell line    Annotator    F1
```

Here, `F1` is used as the plotted value for **FSM+ISM recovery**.

#### `SQANTI_FL_ratio.tsv`
Used for the full-lengthness axis of the scatter plot. It must contain at least:

```text
Chemistry    Cell line    Annotator    F1
```

Here, `F1` is used as the plotted value for **FSM / (FSM + ISM)**.

### What is the output?

The script writes:

- `SQANTI_Figure_ForestA_ForestB_TradeoffC.pdf`

### What does the figure show?

**Panel A** shows whether LoRTIA has a positive or negative effect size relative to the other annotators for **active-catalog recovery (FSM+ISM recall)** within each chemistry.  
**Panel B** shows the same comparison for **full-lengthness among recovered transcripts**.  
**Panel C** places each annotator × chemistry × cell line observation in a **recovery vs full-lengthness** space, with larger centroid points summarizing mean positions by annotator and chemistry.

### What is the message of the figure?

This figure is designed to show two things at once:

1. whether LoRTIA outperforms the other annotators in **recovery** and **full-lengthness**, and  
2. how the methods are positioned in the **trade-off space** between recovering more catalog-supported transcripts and reconstructing them more completely.

------------------------------------------

## `NNC-NIC-visualization.R`

### What does it do?

This R script generates a **multi-panel NIC/NNC novelty figure** focused on catalog-external transcript discovery and support structure. It contains four panels:

- **Panel A:** NIC/NNC burden
- **Panel B:** novel junction support structure
- **Panel C:** canonicality of NNC novelty
- **Panel D:** locus-level consistency of NNC novelty

The complete figure is exported in both PDF and PNG format.

### Example run command

```bash
Rscript NNC-NIC-visualization.R
```

### Which files does it require?

This script uses files written by `summarize_sqanti_novelty_support.py`:

- `nic_nnc_burden_by_condition.tsv`
- `novel_junction_support_by_annotator_chemistry.tsv`
- `novel_junction_plausibility_by_annotator_chemistry.tsv`
- `novel_isoform_canonicality_by_annotator_chemistry_category.tsv`
- `novel_locus_consistency_by_annotator_chemistry_category.tsv`

### What do the input files need to contain?

#### `nic_nnc_burden_by_condition.tsv`
Used for **Panel A**. It must contain at least:

```text
annotator    chemistry    n_total_isoforms    n_nic    n_nnc
```

The script aggregates these counts and converts them to fractions of reconstructed isoforms.

#### `novel_junction_support_by_annotator_chemistry.tsv`
Used for **Panel B**. It must contain at least:

```text
annotator    chemistry    frac_cross_tool    frac_multi_chemistry_single_tool    frac_multi_cell_line_single_tool    frac_single_condition
```

These values are displayed as stacked support classes for novel junctions.

#### `novel_junction_plausibility_by_annotator_chemistry.tsv`
Used for **Panel C**. It must contain at least:

```text
annotator    chemistry    scope    canonical_fraction
```

The script specifically uses `scope == "nnc_only"` to summarize NNC junction canonicality.

#### `novel_isoform_canonicality_by_annotator_chemistry_category.tsv`
Also used for **Panel C**. It must contain at least:

```text
annotator    chemistry    category    canonical_isoform_fraction
```

The script specifically uses `category == "NNC"` to compare junction-level and isoform-level canonicality.

#### `novel_locus_consistency_by_annotator_chemistry_category.tsv`
Used for **Panel D**. It must contain at least:

```text
annotator    chemistry    category    genes_with_independent_supported_novel_junction    genes_with_novelty    frac_genes_with_ge_2_novel_isoforms
```

The script derives the fraction of genes with independently supported novel junctions from these columns.

### What is the output?

The script writes:

- `FigureX_NIC_NNC_discovery.pdf`
- `FigureX_NIC_NNC_discovery.png`

### What does the figure show?

**Panel A** shows the fraction of reconstructed isoforms classified as **NIC** or **NNC** for each annotator and chemistry.  
**Panel B** shows how novel junctions are distributed across support classes, including **cross-tool**, **multi-chemistry**, **multi-cell-line**, and **single-condition** support.  
**Panel C** compares **junction-level canonicality** and **isoform-level canonicality** for NNC novelty.  
**Panel D** summarizes whether NNC novelty accumulates consistently at the locus level, including the fraction of genes with at least two NNC isoforms and the fraction of genes with independently supported NNC junctions.

### What is the message of the figure?

This figure emphasizes that novelty should not be interpreted only by **quantity**, but also by **support structure**, **splice plausibility**, and **locus-level consistency**. In other words, it helps distinguish annotators that produce many novel calls from annotators whose novel calls are better supported and more structurally plausible.

------------------------------------------

## Short summary

In this SQANTI-based benchmark:

- annotator-specific SQANTI3 files are stored under `SQANTI3/<Annotator>/<Chemistry>/<Cell-line>/`
- chemistry × cell line-specific active reference GTF files are stored under `SQANTI3/cell_active_transcript_reference/`
- `build_sqanti_cell_active_references.py` constructs active reference GTFs from SQANTI3 FSM/ISM-supported `associated_transcript` identifiers and GENCODE v48
- `benchmark_sqanti_transcript_recovery.py` computes transcript recovery metrics, statistical comparisons, and plot-ready tables
- `summarize_sqanti_novelty_support.py` summarizes NIC/NNC burden and the support/plausibility of novel transcript structures
- `FSM-ISM-visualization.R` generates the main recovery/full-lengthness comparison figure
- `NNC-NIC-visualization.R` generates the novelty burden, support structure, and plausibility summary figure.
