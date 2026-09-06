# LoRTIA Plus

<p align="center">
  <img src="LoRTIA_Plus/LoRTIA.png" alt="LoRTIA" width="500">
</p>

**Chemistry-agnostic, feature-first long-read transcriptome annotation**

LoRTIA Plus is a long-read RNA-seq transcript annotation workflow for detecting
transcription start sites (TSSs), transcription end sites (TESs), splice
junctions, transcript isoforms, and candidate fusion events from genome-aligned
long reads.

The software is designed around a **feature-first** strategy: transcript
boundaries and introns are evaluated before transcript reconstruction. LoRTIA
Plus does not require a reference transcript annotation to define, correct, or
constrain the TSS, TES, splice-junction, or transcript-model universe.

LoRTIA Plus supports Oxford Nanopore Technologies (ONT) and Pacific Biosciences
(PacBio) long-read RNA-seq data, including cDNA, CapTrap, direct RNA, and
full-length PacBio cDNA workflows when the relevant terminal sequence evidence
has been retained.

- **Software repository:** https://github.com/philipmories/LoRTIA-Plus
- **Archived release:** https://doi.org/10.5281/zenodo.19348065
- **LoRTIA Plus preprint:** https://doi.org/10.64898/2026.04.03.716279
- **Benchmark workflows:** [`benchmark/`](benchmark/)
**Benchmark scripts, documentation, input files, and output files:** https://doi.org/10.6084/m9.figshare.31891948

> **Important**
>
> LoRTIA Plus extends the original LoRTIA implementation substantially.
> The historical 2019 LoRTIA documentation describes the original software and
> should not be used as the primary documentation for LoRTIA Plus. The present
> README, the current `LoRTIA -h` interface, and the supplementary software
> documentation accompanying the LoRTIA Plus manuscript describe the current
> implementation.

---

## Contents

- [Overview](#overview)
- [What is new in LoRTIA Plus?](#what-is-new-in-lortia-plus)
- [Workflow](#workflow)
- [Installation](#installation)
- [Installation verification](#installation-verification)
- [Input requirements](#input-requirements)
- [Quick start](#quick-start)
- [Terminal adapter-search geometry](#terminal-adapter-search-geometry)
- [Adapter score and parameter interactions](#adapter-score-and-parameter-interactions)
- [Library-specific starting configurations](#library-specific-starting-configurations)
- [Sample-level feature detection](#sample-level-feature-detection)
- [Multi-sample feature merging](#multi-sample-feature-merging)
- [Transcript reconstruction](#transcript-reconstruction)
- [Fusion detection](#fusion-detection)
- [Main modules](#main-modules)
- [Outputs](#outputs)
- [Dependency specification](#dependency-specification)
- [Benchmark reproducibility](#benchmark-reproducibility)
- [Documentation](#documentation)
- [Citation](#citation)
- [Legacy LoRTIA](#legacy-lortia)

---

## Overview

LoRTIA Plus separates long-read transcript annotation into explicit feature
detection and reconstruction stages.

The principal workflow is:

```text
genome-aligned SAM/BAM
        |
        v
read-level terminal QC and orientation
        |
        v
sample-level TSS / TES / intron detection
        |
        +----> processed and strand-resolved alignments
        |
        v
cross-sample feature merging with Sum_gffs.py
        |
        v
optional study-specific feature selection
        |
        v
accepted TSS / TES / intron GFF3 feature set
        |
        v
Transcript_Annotator.py
        |
        v
transcript models
```

Fusion-event detection is handled by the integrated
`Fusion_Annotator.py` path, which evaluates multiple primary/supplementary
alignment segments originating from the same read in query-coordinate space.

### Key design principles

- **Feature-first reconstruction.** TSSs, TESs, and introns are evaluated before
  transcript models are assembled.
- **No transcript annotation required for discovery.** LoRTIA Plus reconstructs
  features directly from aligned reads and the reference genome.
- **Terminal-sequence-aware QC.** cDNA workflows can use 5′ adapter and
  poly(A)-associated evidence to evaluate transcript ends and orientation.
- **Chemistry-aware configuration.** Terminal sequences and score thresholds are
  configurable for the library chemistry and preprocessing strategy.
- **Cross-sample consolidation.** Sample-level TSS, TES, and intron calls can be
  merged across replicate libraries before transcript reconstruction.
- **Optional downstream feature selection.** Cross-sample support tables can be
  filtered according to study-specific recurrence, cumulative-support, or
  feature-quality criteria when required.
- **Direct RNA support.** `--dRNA` uses mapper-derived strand information and
  evaluates the true 3′ terminus without requiring 5′ adapter detection.
- **Fusion-event support.** Query-space split alignments can be evaluated as
  candidate structural fusion events.

---
## What is new in LoRTIA Plus?

LoRTIA Plus extends the original LoRTIA software with several major
implementation and functionality updates.

### Polars-based transcript reconstruction

Transcript-processing and reconstruction steps have been modernized using
**Polars-based dataframe processing**. This reduces dataframe-processing
overhead and improves the performance and scalability of transcript
reconstruction, particularly for large long-read RNA-seq datasets.

### Modernized adapter-alignment backend

The previous `Bio.pairwise2`-based adapter-alignment implementation has been
replaced with `Bio.Align.PairwiseAligner`.

This removes dependence on the deprecated `Bio.pairwise2` interface and uses a
maintained alignment backend for terminal adapter and poly(A)-associated
sequence recognition.

### Chemistry-specific processing modes

LoRTIA Plus introduces library- and chemistry-specific processing modes so that
terminal-sequence recognition and strand assignment can be adapted to the
sequence information retained by different long-read RNA-seq protocols.

The current workflow supports:

- ONT cDNA / direct-cDNA;
- ONT CapTrap;
- PacBio cDNA / Iso-Seq with retained 5′ primer;
- PacBio CapTrap;
- ONT direct RNA through `--dRNA`;
- full-length PacBio cDNA after 5′ primer removal through `--pacbio`, provided
  that detectable terminal poly(A)-associated sequence is retained.

In `--dRNA` mode, LoRTIA Plus does not require 5′ adapter detection and uses
mapper-derived strand information. In `--pacbio` mode, both read ends are
searched for the configured 3′ signal; a single accepted 3′ call determines
orientation and the opposite end is treated as the 5′ end.

### Integrated fusion-event analysis

LoRTIA Plus includes an integrated fusion-analysis workflow through
`Fusion_Annotator.py`.

Primary and supplementary alignment segments originating from the same read
are evaluated together in query-coordinate space. Candidate split-read chains
are assessed using configurable thresholds for query overlap, unaligned query
gap, breakpoint clustering tolerance, and minimum distinct-read support.

The current default fusion settings are:

```text
maximum split-query overlap: 15 nt
maximum query gap:           20 nt
breakpoint wobble:           50 nt
minimum supporting reads:     1
```

The fusion workflow can report structural fusion candidates and, when
compatible transcript-boundary, intron, and orientation information is
available, reconstruct complete fusion-transcript models.

---

## Installation

LoRTIA Plus is intended for UNIX-like systems, including Linux and Windows
Subsystem for Linux (WSL).

### Recommended reproducible installation with Conda

Clone the repository and check out the corresponding release:

```bash
git clone https://github.com/philipmories/LoRTIA-Plus.git
cd LoRTIA-Plus
git checkout v1.0.0
```

Create and activate the pinned environment:

```bash
conda env create -f environment.yml
conda activate lortia-plus
```

Install LoRTIA Plus into the environment:

```bash
python -m pip install . --no-deps --no-build-isolation
```

The `--no-deps` and `--no-build-isolation` options keep the installation within
the versions already provided by the pinned Conda environment.

### Standard pip installation

After cloning or downloading the source repository:

```bash
python -m pip install .
```

A direct installation from the tagged GitHub release is also supported by the
package metadata:

```bash
python -m pip install "git+https://github.com/philipmories/LoRTIA-Plus.git@v1.0.0"
```

When Conda is not used, the external command-line dependencies `samtools` and
GNU `awk` must be installed separately and available in `PATH`.

---

## Installation verification

After installation:

```bash
python verify_install.py
LoRTIA -h
```

The supplied `verify_install.py` checks:

- Python version
- pinned Python dependency versions
- `samtools`
- GNU `awk`
- Python user-site isolation
- availability of the installed `LoRTIA` command
- successful execution of the LoRTIA Plus help interface

A fully matching reference environment reports:

```text
LoRTIA Plus installation check passed.
```

For complete validation, run LoRTIA Plus end-to-end on a representative SAM or
BAM dataset using the library-appropriate terminal configuration.

---

## Input requirements

LoRTIA Plus takes:

1. a genome-aligned **SAM or BAM** file,
2. an output directory, and
3. the **reference FASTA** used for mapping.

General invocation:

```bash
LoRTIA [options] input_file output_path reference_fasta
```

The reference FASTA must correspond to the reference used to generate the
alignment.

### Preserve terminal query sequence

Terminal adapter/poly(A)-associated evidence is evaluated from query sequence
retained at alignment ends. Therefore, alignments used for terminal recognition
should preserve terminal query sequence.

For minimap2-based BAM workflows, a suitable starting configuration is:

```bash
minimap2 -ax splice -Y --secondary=no -t <threads> \
  reference.fasta reads.fastq |
  samtools sort -@ <threads> -o sample.sorted.bam -

samtools index sample.sorted.bam
```

For ONT direct RNA:

```bash
minimap2 -ax splice -uf -k14 -Y --secondary=no -t <threads> \
  reference.fasta reads.fastq |
  samtools sort -@ <threads> -o sample.sorted.bam -

samtools index sample.sorted.bam
```

`-Y` preserves supplementary-alignment query sequence as soft clipping rather
than hard clipping, while `--secondary=no` suppresses secondary mappings but
retains supplementary split alignments.

---

## Quick start

The following commands are **recommended starting configurations**, not
universal presets. Adapter sequences must match the actual library preparation
and the sequence retained after preprocessing.

### ONT cDNA / direct-cDNA

```bash
LoRTIA \
  -5 TGCCATTAGGCCGGG \
  --five_score 16 \
  --check_in_soft 15 \
  --check_in_match 10 \
  --check_from_alignment 3 \
  -3 AAAAAAAAAAAAAAA \
  --three_score 16 \
  --shs_for_ts 3 \
  -s poisson \
  sample.sorted.bam output_folder reference.fasta
```

### ONT CapTrap

```bash
LoRTIA \
  -5 AACGCAGAGTAC \
  --five_score 16 \
  --check_in_soft 15 \
  --check_in_match 10 \
  --check_from_alignment 3 \
  -3 AAAAAAAAAAAAAAA \
  --three_score 16 \
  --shs_for_ts 3 \
  -s poisson \
  sample.sorted.bam output_folder reference.fasta
```

### PacBio cDNA / Iso-Seq with 5′ primer retained

```bash
LoRTIA \
  -5 AGAGTACATGGG \
  --five_score 16 \
  --check_in_soft 15 \
  --check_in_match 10 \
  --check_from_alignment 3 \
  -3 AAAAAAAAAAAAAAA \
  --three_score 18 \
  --shs_for_ts 3 \
  -s poisson \
  sample.sorted.bam output_folder reference.fasta
```

### PacBio full-length cDNA after 5′ primer removal

Use this mode only when the 5′ primer has been removed but detectable terminal
poly(A)-associated sequence remains:

```bash
LoRTIA \
  --pacbio \
  --check_in_soft 15 \
  --check_in_match 10 \
  --check_from_alignment 3 \
  -3 AAAAAAAAAAAAAAA \
  --three_score 18 \
  --shs_for_ts 3 \
  -s poisson \
  sample.sorted.bam output_folder reference.fasta
```

In `--pacbio` mode, both ends are searched for the configured 3′ signal.
Exactly one correct 3′ call is required to infer orientation; the opposite end
is then accepted as the 5′ end without requiring a 5′ adapter hit.

Full-length status is an input assumption in this mode. If preprocessing has
also removed the terminal poly(A) signal, `--pacbio` cannot restore the missing
terminal evidence.

### PacBio CapTrap

```bash
LoRTIA \
  -5 AACGCAGAGTAC \
  --five_score 16 \
  --check_in_soft 15 \
  --check_in_match 10 \
  --check_from_alignment 3 \
  -3 AAAAAAAAAAAAAAA \
  --three_score 16 \
  --shs_for_ts 3 \
  -s poisson \
  sample.sorted.bam output_folder reference.fasta
```

### ONT direct RNA

```bash
LoRTIA \
  --dRNA \
  --check_in_soft 6 \
  --check_in_match 10 \
  --check_from_alignment 3 \
  -3 AAAA \
  --three_score 6 \
  --shs_for_ts 12 \
  -s poisson \
  sample.sorted.bam output_folder reference.fasta
```

In `--dRNA` mode, LoRTIA Plus does not perform 5′ adapter detection. Strand is
taken from mapper-derived information, the true 5′ end is accepted, and
terminal quality control is applied to the true 3′ end.

---

## Terminal adapter-search geometry

The three terminal-search parameters below control **different parts of the
same read-space search geometry**.

```text
S S S S S S S S S S S S S S S S | M M M M M M M M M M M M M M M M
                                  ^
                                  mapped boundary

      <------ --check_in_soft ---->|<-- --check_in_match -->
      <----------- adapter-search window ------------------>
```

- `S` = terminal soft-clipped query sequence
- `M` = query sequence already included in the mapped alignment

### `--check_in_soft`

Defines the maximum number of terminal **soft-clipped query nucleotides**
included in the local adapter-search window immediately adjacent to the mapped
boundary.

A smaller search depth reduces exposure to internal adapter-like sequence but
can remove usable terminal sequence if set too stringently.

**Generic default:** `30 nt`

### `--check_in_match`

Extends the same search window into the immediately adjacent **aligned read
sequence inside the mapped boundary**.

This is important when adapter-derived bases have been partially aligned to the
reference rather than represented entirely as soft clipping.

These bases do **not** need to be perfect reference matches, and
`--check_in_match` is **not** a required adapter-match count.

**Default:** `10 nt`

### `--check_from_alignment`

After local alignment of the expected terminal sequence,
`--check_from_alignment` constrains the physical placement of the adapter hit.

It is the maximum number of **read bases** by which the **boundary-facing end**
of an accepted local adapter hit may stop **before the mapped boundary on the
clipped side**.

It is:

- a read-space positional constraint,
- **not** a genomic distance,
- **not** a distance from the read tip.

Adapter hits extending into the mapped side are limited by the aligned-side
sequence made available through `--check_in_match`.

**Default:** `3 nt`

### Worked example

With:

```text
--check_in_soft 15
--check_in_match 10
--check_from_alignment 3
```

a 20-nt terminal soft clip contributes the 15 bases nearest the mapped
boundary, and up to 10 adjacent aligned read bases are appended to the search
window.

- adapter hit ends **3 read bases before** the mapped boundary → positional test
  passes;
- adapter hit ends **4 read bases before** the mapped boundary → positional test
  fails regardless of alignment score;
- adapter hit extends **2 bases into the aligned side** → geometrically allowed
  because those bases are inside the `--check_in_match` portion of the search
  window, although the terminal artefact checks still apply;
- hard-clipped query sequence cannot be recovered from the alignment record.

The conventional adapter-search path requires terminal soft clipping:
`--check_in_match` alone does not turn an entirely unclipped read end into an
adapter-searchable end.

---

## Adapter score and parameter interactions

Terminal recognition uses local Smith-Waterman alignment implemented with
`Bio.Align.PairwiseAligner`.

The generic scoring configuration is:

| Component | Score |
|---|---:|
| Match | +2 |
| Mismatch | -3 |
| Gap opening | -3 |
| Gap extension | -3 |

The best local-alignment score must be **strictly greater** than
`--five_score` or `--three_score`. Equality fails.

Example with a score threshold of `16`:

- 8 perfectly paired bases → score `16` → **fails**
- 9 perfectly paired bases → score `18` → **passes the score test**

The positional and terminal QC criteria must also pass.

Parameter interactions matter:

- increasing `--check_in_soft` exposes more terminal soft-clipped sequence to
  possible adapter hits;
- increasing `--check_in_match` can recover partly aligned adapters but also
  changes the mapped-side sequence available to terminal short-homology checks;
- increasing `--check_from_alignment` permits more distant terminal hits;
- adapter length and score thresholds should be selected together;
- score, search-window size, and positional tolerance are not interchangeable.

### Terminal short homology: `--shs_for_ts`

`--shs_for_ts` controls a read-end short-homology test used to flag potential
reverse-transcription template-switching configurations.

The test evaluates adapter similarity extending into adjacent mapped query
sequence. Potential-template-switching observations are retained as a separate
read-end class for downstream QC.

The threshold should generally be smaller than `--check_in_match`; values above
`--check_in_match` make this read-end short-homology test inactive.

---

## Library-specific starting configurations

| Library | 5′ terminal search | 3′ terminal search | `check_in_soft` | Notes |
|---|---|---|---:|---|
| ONT cDNA / direct-cDNA | `TGCCATTAGGCCGGG`, score 16 | `A`×15, score 16 | 15 | `shs_for_ts=3` |
| ONT CapTrap | `AACGCAGAGTAC`, score 16 | `A`×15, score 16 | 15 | `shs_for_ts=3` |
| PacBio cDNA, 5′ primer retained | `AGAGTACATGGG`, score 16 | `A`×15, score 18 | 15 | `shs_for_ts=3` |
| PacBio CapTrap | `AACGCAGAGTAC`, score 16 | `A`×15, score 16 | 15 | `shs_for_ts=3` |
| PacBio cDNA, 5′ primer removed | no 5′ adapter search; use `--pacbio` | `A`×15, score 18 | 15 | requires retained terminal poly(A)-associated evidence |
| ONT direct RNA | no 5′ adapter search | `AAAA`, score 6 | 6 | `--dRNA`, `shs_for_ts=12` |

These values are starting configurations. Terminal sequences should always be
matched to the actual library preparation and preprocessing history.

---

## Sample-level feature detection

The top-level LoRTIA workflow performs read-level terminal QC and then evaluates
TSSs, TESs, and introns at sample level.

### TSS and TES

Candidate transcript boundaries are evaluated at nucleotide resolution.

Core starting settings include:

| Parameter | Default / recommended starting value | Role |
|---|---:|---|
| `wobble` | ±10 nt | local neighbourhood for representative boundary selection |
| `minimum` | 2 reads | minimum absolute feature support |
| `ratio` | 0.001 | feature support relative to local transcript coverage |
| `distance` | 15 nt | inward offset used for local coverage sampling |
| `cov_sample` | 5 positions | number of positions averaged for local coverage |
| Poisson TSS test | optional; recommended with `-s poisson` | local-background significance test for TSSs |

When `-s poisson` is enabled, the Poisson significance test is applied to TSS
generation. TES generation does not use the TSS Poisson filter.

### Introns

Splice junctions are extracted from `N` CIGAR operations after read-level
artefact checks.

Core starting settings include:

| Parameter | Default / starting value | Role |
|---|---:|---|
| `intron_wobble` | ±15 nt at donor and acceptor | defines nearby alternative junctions |
| `rare_intron` | 0.05 | minimum support relative to the strongest nearby junction |
| `force_consensus` | False | optionally restricts retained introns to consensus splice-site pairs |

GT-AG, GC-AG, AT-AC, and corresponding reverse-strand motifs are annotated as
consensus splice-site pairs. Consensus annotation and consensus-only exclusion
are separate operations. When `--force_consensus` is disabled, splice-site
motif status can still be annotated but is not used as a hard exclusion
criterion.

Intron-boundary short-homology information is also recorded as junction-level
QC information and can be used in study-specific downstream feature selection
when required.

---

## Multi-sample feature merging

When multiple libraries or replicates represent the same biological dataset,
sample-level TSS, TES, and intron features can be consolidated with
`Sum_gffs.py`.

`listfile.txt` contains one sample-level LoRTIA output prefix per line.

### TSS

```bash
python3 LoRTIA_Plus/Sum_gffs.py listfile.txt merged_output tss
```

### TES

```bash
python3 LoRTIA_Plus/Sum_gffs.py listfile.txt merged_output tes
```

### Introns

```bash
python3 LoRTIA_Plus/Sum_gffs.py listfile.txt merged_output intron
```

For TSSs and TESs, nearby coordinates are consolidated within the configured
endpoint wobble. Introns are merged by donor-acceptor coordinate identity.

The resulting feature-specific support tables record the merged coordinate and
sample-specific support.

### Optional study-specific feature selection

The merged support tables allow the investigator to define the feature universe
used for reconstruction.

Possible criteria include, for example:

- minimum number of supporting samples,
- cumulative read support,
- feature-specific quality annotations.

These criteria are **study-specific**, not mandatory LoRTIA Plus defaults.

For the LRGASP and KSHV benchmark analyses reported with LoRTIA Plus, no
additional post-merge minimum sample-recurrence, cumulative-support, SHS,
consensus-motif, or other feature-quality threshold was applied: all merged
features were retained for transcript reconstruction.

---

## Transcript reconstruction

`Transcript_Annotator.py` reconstructs transcript models from a processed
LoRTIA Plus BAM and a common accepted-feature prefix.

General form:

```bash
python3 LoRTIA_Plus/Transcript_Annotator.py \
  sample_processed.bam accepted_feature_prefix sample_output_prefix
```

The accepted feature prefix refers to:

```text
accepted_feature_prefix_tss.gff3
accepted_feature_prefix_tes.gff3
accepted_feature_prefix_intron.gff3
```

Read boundaries are matched to accepted TSS/TES coordinates within the
configured endpoint tolerance. Mapper-defined introns must match accepted
donor-acceptor coordinates exactly.

Compatible reads are collapsed according to:

- reference sequence,
- strand,
- accepted TSS,
- accepted TES,
- ordered intron chain.

Standalone reconstruction defaults include:

| Parameter | Default |
|---|---:|
| endpoint `wobble` | 10 nt |
| unsupported `gap` | 1 nt |
| `mintr_count` | 1 read |

---

## Fusion detection

LoRTIA Plus includes structural fusion-event detection through
`Fusion_Annotator.py`.

Multiple primary/supplementary alignment segments sharing the same read name
are evaluated in query-coordinate space.

Current top-level defaults:

| Parameter | Default | Meaning |
|---|---:|---|
| `--max_split_query_overlap` | 15 nt | maximum allowed overlap between retained split segments |
| `--fusion_max_query_gap` | 20 nt | maximum unaligned query gap between adjacent retained segments |
| `--fusion_wobble` | 50 nt | breakpoint-clustering tolerance |
| `--fusion_min_reads` | 1 | minimum number of distinct supporting reads |

Structural event reporting is adapter-independent. Complete reconstructed
fusion-transcript models additionally require compatible external TSS/TES and
intron features and resolved transcript orientation.

The fusion detector was technically validated using the JAFFAL ONT_85
simulation; benchmark-specific commands and outputs are documented in the
supplementary material and benchmark resources.

---

## Main modules

| File / command | Principal role |
|---|---|
| `LoRTIA` | top-level command-line workflow and orchestration |
| `LoRTIA_Plus/Samprocessor.py` | read-level alignment handling, terminal-sequence QC, orientation, and primary feature tables |
| `LoRTIA_Plus/Stats.py` | sample-level TSS, TES, and intron evaluation |
| `LoRTIA_Plus/Gff_creator.py` | feature GFF3 generation |
| `LoRTIA_Plus/Sum_gffs.py` | cross-sample feature consolidation and support matrices |
| `LoRTIA_Plus/Transcript_Annotator.py` | transcript reconstruction from accepted features |
| `LoRTIA_Plus/Fusion_Annotator.py` | split-alignment fusion-event detection and fusion-transcript reconstruction |

Use:

```bash
LoRTIA -h
```

for the complete current top-level CLI parameter reference, defaults, parameter
interactions, geometry explanation, and worked examples.

---

## Outputs

Depending on the selected workflow and library type, LoRTIA Plus produces
sample-level and reconstruction-level outputs including:

- processed alignments,
- strand-resolved / oriented alignments,
- terminal read-end and QC tables,
- local coverage tables,
- TSS statistics and TSS GFF3,
- TES statistics and TES GFF3,
- intron statistics and intron GFF3,
- merged cross-sample feature-support tables,
- accepted-feature GFF3 files,
- reconstructed transcript GFF3,
- transcript TSV output,
- annotated BAM output,
- structural fusion-candidate tables,
- reconstructed fusion-transcript outputs.

The exact files generated depend on the selected library mode and analysis path.

---

## Dependency specification

The complete pinned reference environment is provided in `environment.yml`.

### Python runtime

| Dependency | Pinned version |
|---|---:|
| Python | 3.10.14 |
| NumPy | 1.26.4 |
| pandas | 2.2.3 |
| Polars | 1.8.2 |
| pysam | 0.22.1 |
| SciPy | 1.14.1 |
| Biopython | 1.88 |

### External command-line dependencies

| Dependency | Pinned version | Purpose |
|---|---:|---|
| samtools | 1.20 | SAM/BAM processing |
| GNU awk | 5.3.1 | text processing used by the workflow |

### Build / installation tools

| Dependency | Pinned version |
|---|---:|
| pip | 24.2 |
| setuptools | 75.1.0 |
| wheel | 0.44.0 |

`bedtools` is **not required** by the current LoRTIA Plus implementation.

The previous `Bio.pairwise2` adapter-alignment implementation has been replaced
with `Bio.Align.PairwiseAligner`; the current software does not depend on the
deprecated `Bio.pairwise2` interface.

Files supporting reproducible installation:

- `environment.yml`
- `requirements.txt`
- `pyproject.toml`
- `verify_install.py`

---

## Benchmark reproducibility

The repository also contains benchmark workflows under:

```text
benchmark/
```

The LoRTIA Plus manuscript evaluates the software across:

- KSHV long-read transcriptome datasets,
- human LRGASP datasets,
- multiple ONT and PacBio library chemistries,
- independent simulation-based analyses,
- a targeted JAFFAL fusion-event validation.

Benchmark-specific effective parameters and commands should be taken from the
benchmark documentation and manuscript supplementary material rather than
assuming that every benchmark setting is a universal software default.

In particular, the published benchmark workflow used the recommended
`-s poisson` TSS configuration but did **not** enable
`--force_consensus`, and no additional post-merge feature recurrence/count/QC
filter was applied to the merged LRGASP or KSHV feature sets.

---

## Documentation

Current LoRTIA Plus documentation is provided in three complementary forms:

1. **This README** — overview, installation, workflow, parameter geometry,
   chemistry-specific starting configurations, and main modules.
2. **`LoRTIA -h`** — complete command-line reference with current defaults,
   parameter interactions, geometry explanation, and worked examples.
3. **LoRTIA Plus supplementary software documentation** — detailed
   implementation-level methods, feature-selection logic, benchmark-specific
   settings, and reproducibility information.

The historical original-LoRTIA wiki is not the authoritative documentation for
LoRTIA Plus.

---

## Citation

If you use LoRTIA Plus, please cite the LoRTIA Plus work and the archived
software release.

### LoRTIA Plus manuscript / preprint

Torma G, Balázs Z, Fülöp Á, Tombácz D, Boldogkői Z.  
**LoRTIA Plus: a chemistry-agnostic, feature-first software package for
long-read transcriptome annotation.**  
bioRxiv (2026).  
https://doi.org/10.64898/2026.04.03.716279

### Software archive

**LoRTIA-Plus**  
Zenodo: https://doi.org/10.5281/zenodo.19348065

---

## Legacy LoRTIA

LoRTIA Plus is derived from the original LoRTIA toolkit developed by
Zsolt Balázs and collaborators.

Original repository:

https://github.com/zsolt-balazs/LoRTIA

The original repository and its 2019 wiki document the historical LoRTIA
implementation. They remain useful for provenance and software history, but
users of LoRTIA Plus should follow the current LoRTIA Plus documentation and
CLI help described above.
