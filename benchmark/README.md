# Benchmark workflows

This directory contains the benchmark workflows used in the **LoRTIA Plus** study to evaluate transcript annotation performance across **viral** and **human long-read RNA-seq** benchmark systems. The benchmark code is organized into two top-level modules:

- **`LRGASP/`** — human long-read transcriptome benchmarks based on LRGASP datasets
- **`KSHV/`** — viral benchmark workflows for the KSHV transcriptome

The **`LRGASP/`** module contains three complementary workflows:

- **`boundary_benchmark/`** — benchmarking of **TSS** and **TES** detection across annotators, sequencing chemistries, and cell lines
- **`SQANTI3_benchmark/`** — SQANTI-based evaluation of recovery of known catalog-supported transcript structures, together with novelty-support and plausibility analyses
- **`SQANTI3_basic_stats/`** — summary statistics and plotting workflows for isoform yield, structural class composition, isoform length distributions, coding/non-coding composition, and gene-level isoform complexity

The **`KSHV/`** module contains the benchmark workflow used to compare **LoRTIA Plus** with **bambu**, **FLAIR**, **IsoQuant**, and **NAGATA** in the compact, transcriptionally dense **KSHV** transcriptome using **ONT dcDNA** and **ONT dRNA** datasets. This benchmark evaluates **TSS**, **TES**, and **transcript-level** recovery and includes figure-generation code for the composite KSHV benchmark figure.

This GitHub repository is the primary location for the **benchmarking scripts, workflow documentation, manifests, and figure-generation code**. Large benchmark input/output files and archived companion datasets are stored separately on **Figshare**.
