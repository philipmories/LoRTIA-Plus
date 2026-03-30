# LoRTIA Plus benchmark companion dataset

This Figshare record accompanies the manuscript **“LoRTIA Plus: a chemistry-agnostic, feature-first software package for long-read transcriptome annotation”** and provides the benchmark-associated annotation resources used to compare **LoRTIA Plus** with **bambu**, **FLAIR**, **IsoQuant**, and **NAGATA** across two complementary long-read transcriptomics systems:

- the **KSHV transcriptome**
- **human LRGASP long-read RNA-seq datasets**

Because of platform file-count limitations, the dataset is distributed here in a simplified archive-based form.

## Files included in this Figshare record

- **`KSHV.rar`**  
  Contains the KSHV-specific benchmark annotation resources used in the viral comparison workflow. This package is centered on the compact, transcriptionally dense KSHV transcriptome benchmark performed on **ONT dcDNA** and **ONT dRNA** datasets.

- **`LRGASP_SQANTI3.rar`**  
  Contains the human benchmark annotation resources used in the LRGASP-based analyses. This package includes annotator-specific **GTF** files together with the corresponding **SQANTI3 `classification.txt`** and **`junctions.txt`** outputs for the human transcriptome benchmark across **three cell lines** and **five long-read sequencing chemistries**.

## Benchmark scope

The study benchmarked **LoRTIA Plus** against **bambu**, **FLAIR**, **IsoQuant**, and **NAGATA** in two major settings:

### 1. KSHV benchmark
The viral benchmark evaluated transcript annotation performance in a compact and transcriptionally dense herpesvirus genome. Analyses were performed separately for:

- **ONT direct-cDNA (dcDNA)**
- **ONT direct-RNA (dRNA)**

This part of the study focused on:

- **TSS recovery**
- **TES recovery**
- **transcript-level recovery**
- false-positive transcript-model overlap

### 2. Human LRGASP benchmark
The human benchmark was based on long-read RNA-seq data from the **LRGASP** framework and covered:

- **three cell lines**
  - **H1-hES**
  - **H1-DE**
  - **WTC11**
- **five sequencing chemistries**
  - **ONT-dRNA**
  - **ONT-cDNA**
  - **ONT-CapTrap**
  - **PacBio cDNA**
  - **PacBio CapTrap**

This part of the study included:

- recovery of known transcript structures
- SQANTI3-based structural classification
- FSM/ISM-based transcript recovery
- NIC/NNC novelty analyses
- splice-junction support and plausibility analyses
- isoform-level and gene-level summary analyses

## Package contents

### `KSHV.rar`
This archive contains the KSHV-related benchmark files. Depending on the final packaged version, it may include:

- annotator-specific **GFF3** predictions
- benchmark-restricted / active reference files
- KSHV reference transcript and feature resources
- optional figure-input tables
- a detailed KSHV-specific README describing the associated scripts and workflow logic

### `LRGASP_SQANTI3.rar`
This archive contains the human benchmark files. Depending on the final packaged version, it may include:

- annotator-specific **GTF** files
- **SQANTI3 `classification.txt`** files
- **SQANTI3 `junctions.txt`** files
- chemistry × cell line-specific reference files
- figure-input tables
- workflow-specific README files for the SQANTI-based benchmark and summary-statistics modules

## Code and workflow availability

This Figshare record is a **benchmark companion dataset** and is not intended to serve as the primary archive for the complete computational workflows.

All benchmarking scripts, workflow documentation, statistical summarization code, and figure-generation scripts are available separately in the associated **GitHub repository**.

The detailed README files describing the individual workflows and scripts are also included within the corresponding archive packages and/or the GitHub repository.

## Raw sequencing data

Raw sequencing datasets are not duplicated in this Figshare record. They remain available from their original public repositories, as described in the manuscript Data Availability section.

## Recommended use

Use this Figshare record for:

- accessing the benchmarked transcript annotation files
- inspecting the KSHV and human SQANTI3 structural outputs
- reusing the benchmark input resources associated with the manuscript

Use the associated GitHub repository for:

- running the workflows
- reproducing the benchmark summaries and figures
- inspecting script-level input/output definitions and execution details

## Citation

Please cite the associated manuscript when using these data:

**Torma G, Balázs Z, Fülöp Á, Tombácz D, Boldogkői Z. LoRTIA Plus: a chemistry-agnostic, feature-first software package for long-read transcriptome annotation.**
