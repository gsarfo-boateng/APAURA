# APAURA

**Alternative Polyadenylation Architecture, Usage, and Regulatory Analysis**

[![Version](https://img.shields.io/badge/version-0.10.1-blue.svg)](#citation)
[![Java](https://img.shields.io/badge/Java-8%2B-orange.svg)](#portable-java-jar)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

APAURA is a command-line toolkit for detecting, quantifying, classifying, and visualizing **alternative polyadenylation (APA)** from multiple RNA-sequencing data types.

APAURA supports:

- protocol-aware 3'-end sequencing;
- direct RNA / true 3'-end sequencing;
- IsoQuant long-read corrected reads;
- conventional short-read RNA-seq;
- single-cell APA analysis;
- RNA-binding-protein motif enrichment;
- publication-oriented APA locus visualization.

```text
                                               ___    ____   ___   _   _  ____    ___
                                              / _ \  |  _ \ / _ \ | | | ||  _ \  / _ \
                                             / /_\ \ | |_) / /_\ \| | | || |_) |/ /_\ \
                                            /  _  _ \|  __/  _  _ | |_| ||  _ < |  _  |
                                            \_/ \_/ /|_|   \_/ \_/ \___/ |_| \_\|_| |_|

                                            APAURA v0.10.1
                                            Alternative Polyadenylation Architecture,
                                            Usage, and Regulatory Analysis

                                            modes:
                                            run-3seq | run-direct | run-isoquant
                                            short-read | sc
                                            plot-region | rbp-enrich | ends
```

---

## Installation

APAURA v0.10.1 is distributed as **precompiled native executables for macOS and Linux** and as a **portable Java JAR**.

For most users, the native executable is recommended because **Java, Maven, and GraalVM are not required**.

Download the appropriate file from:

**[APAURA v0.10.1 GitHub Release](https://github.com/gsarfo-boateng/APAURA/releases/tag/v0.10.1)**

### macOS Apple Silicon

For Apple Silicon Macs including M1, M2, M3, and M4:

Download:

```text
apaura-v0.10.1-macos-arm64.tar.gz
```

Extract:

```bash
tar -xzf apaura-v0.10.1-macos-arm64.tar.gz
```

Make APAURA executable:

```bash
chmod +x apaura
```

Check the installation:

```bash
./apaura version
```

Show help:

```bash
./apaura --help
```

If macOS blocks the downloaded executable because it is not notarized, remove the quarantine attribute:

```bash
xattr -d com.apple.quarantine apaura
```

Then run:

```bash
./apaura --help
```

### macOS Intel

Download:

```text
apaura-v0.10.1-macos-x64.tar.gz
```

Extract:

```bash
tar -xzf apaura-v0.10.1-macos-x64.tar.gz
```

Make APAURA executable:

```bash
chmod +x apaura
```

Check the installation:

```bash
./apaura version
```

### Linux x86-64

For standard Intel/AMD Linux systems:

Download:

```text
apaura-v0.10.1-linux-x64.tar.gz
```

Extract:

```bash
tar -xzf apaura-v0.10.1-linux-x64.tar.gz
```

Make APAURA executable:

```bash
chmod +x apaura
```

Check the installation:

```bash
./apaura version
```

### Linux ARM64

Download:

```text
apaura-v0.10.1-linux-arm64.tar.gz
```

Extract:

```bash
tar -xzf apaura-v0.10.1-linux-arm64.tar.gz
```

Make APAURA executable:

```bash
chmod +x apaura
```

Check the installation:

```bash
./apaura version
```

### Add APAURA to your PATH

To run APAURA from any directory without typing `./apaura`, place the executable in a directory on your `PATH`.

For a user-level installation:

```bash
mkdir -p "$HOME/bin"
```

```bash
mv apaura "$HOME/bin/"
```

Add the directory to your shell configuration if necessary:

```bash
export PATH="$HOME/bin:$PATH"
```

For `zsh`, this can be added permanently with:

```bash
echo 'export PATH="$HOME/bin:$PATH"' >> ~/.zshrc
```

Reload the shell:

```bash
source ~/.zshrc
```

APAURA can then be called directly:

```bash
apaura version
```

```bash
apaura --help
```

For a system-wide Linux installation:

```bash
sudo mv apaura /usr/local/bin/apaura
```

### Verify SHA-256 checksums

Each release archive and JAR is accompanied by a `.sha256` file.

On macOS:

```bash
shasum -a 256 -c apaura-v0.10.1-macos-arm64.tar.gz.sha256
```

On Linux:

```bash
sha256sum -c apaura-v0.10.1-linux-x64.tar.gz.sha256
```

---

## Portable Java JAR

APAURA is also distributed as:

```text
APAURA-v0.10.1.jar
```

The portable JAR requires **Java 8 or newer**.

Check Java:

```bash
java -version
```

Check APAURA:

```bash
java -jar APAURA-v0.10.1.jar version
```

Show help:

```bash
java -jar APAURA-v0.10.1.jar --help
```

For example:

```bash
java -jar APAURA-v0.10.1.jar run-3seq --samples samples.tsv --gtf annotation.gtf --protocol quantseq-rev --control control --case treated --out results
```

The native executable and Java JAR expose the same APAURA command interface.

---

## Build from source

Most users do **not** need to build APAURA from source.

Developers can clone the repository:

```bash
git clone https://github.com/gsarfo-boateng/APAURA.git
```

```bash
cd APAURA
```

Build and test:

```bash
mvn clean package
```

APAURA also provides:

```bash
./build.sh
```

The standard Maven build produces:

```text
target/APAURA.jar
```

To build a native executable, use a GraalVM distribution providing `native-image`:

```bash
mvn clean package -Pnative
```

The native build produces:

```text
target/apaura
```

---

## Quick start

After installation:

```bash
apaura --help
```

Available top-level commands include:

```text
ends
run-direct
run-3seq
run-isoquant
short-read
run-short
sc
rbp-enrich
plot-region
citation
version
```

Command-specific help is available with:

```bash
apaura run-3seq --help
```

```bash
apaura run-direct --help
```

```bash
apaura run-isoquant --help
```

```bash
apaura short-read --help
```

```bash
apaura sc --help
```

```bash
apaura rbp-enrich --help
```

```bash
apaura plot-region --help
```

---

## Main capabilities

| Task | Input | APAURA command |
|---|---|---|
| Protocol-aware 3'-end sequencing APA | Aligned BAM/SAM + GTF | `run-3seq` |
| Direct RNA / true 3'-end APA | BAM/SAM + GTF | `run-direct` |
| IsoQuant long-read APA | IsoQuant corrected reads + GTF | `run-isoquant` |
| Conventional short-read RNA-seq APA | BAM/SAM + GTF | `short-read` |
| Single-cell APA | Single-cell alignments + PAS reference + metadata | `sc` |
| Publication-oriented locus visualization | APAURA results or BAM/SAM + GTF | `plot-region` |
| RBP motif enrichment | APAURA results + FASTA + MEME motifs | `rbp-enrich` |
| Transcript-end extraction | BAM/SAM | `ends` |

---

## APA model

APAURA follows a layered model:

```text
poly(A) sites
      ↓
poly(A)-site usage
      ↓
terminal APA architecture
      ↓
mechanism-specific differential APA
      ↓
candidate regulatory mechanisms
```

APAURA separates **terminal APA architecture** from the **APA mechanism that changes between biological conditions**.

This allows tandem 3' UTR APA, alternative last exon usage, and intronic polyadenylation to remain biologically distinct rather than being collapsed into a single proximal-versus-distal score.

### Terminal APA architecture

| Architecture | Interpretation |
|---|---|
| `SINGLE_PAS` | One resolved terminal exon group with one retained terminal PAS |
| `TUTR` | One resolved terminal exon group with multiple retained PAS |
| `ALE` | Multiple resolved terminal exon groups without tandem PAS within a group |
| `mixed` / `MIXED` | Multiple terminal exon groups with tandem PAS in at least one group |
| `UNRESOLVED` | Retained PAS cannot be assigned completely to reconciled terminal structure |

**IPA is not an `apa_architecture` value.**

Intronic polyadenylation is quantified as a separate termination mechanism. A gene can therefore show IPA gain or loss independently of its terminal architecture.

### Direction of change

```text
TUTR → shortening / lengthening / no_directional_shift
ALE  → proximal_switch / distal_switch / switch
IPA  → gain / loss
multiple changing mechanisms → complex
```

Mechanism-specific fields such as `tutr_direction`, `ale_direction`, and `ipa_direction` preserve simultaneous effects occurring within the same gene.

---

# Analysis modes

## 1. Protocol-aware 3'-end sequencing

Use `run-3seq` for aligned libraries in which read-end geometry represents or closely tracks transcript cleavage/polyadenylation sites.

Examples include:

```text
3SEQ
3P-seq
PAS-seq
PolyA-seq
3'READS
3'READS+
PAT-seq
A-seq
A-seq2
SAPAS
3T-fill
PAPERCLIP
2P-seq
PA-seq
WTTS-seq
MAPS
3PC
EXPRESS
QuantSeq
```

Supported protocol names are:

```text
3seq
3p-seq
pas-seq
polyaseq
3reads
3reads+
pat-seq
a-seq
a-seq2
sapas
3t-fill
paperclip
2p-seq
pa-seq
wtts-seq
maps
3pc
exprss
quantseq-fwd
quantseq-rev
quantseq-ffpe
generic
```

Example:

```bash
apaura run-3seq --samples samples.tsv --gtf annotation.gtf --protocol quantseq-rev --control control --case treated --out results_3seq
```

Generic endpoint mode:

```bash
apaura run-3seq --samples samples.tsv --gtf annotation.gtf --protocol generic --control control --case treated --out results_3seq
```

### Sample sheet

```text
sample	condition	bam
CTL_1	control	/absolute/path/CTL_1.bam
CTL_2	control	/absolute/path/CTL_2.bam
CTL_3	control	/absolute/path/CTL_3.bam
TRT_1	treated	/absolute/path/TRT_1.bam
TRT_2	treated	/absolute/path/TRT_2.bam
TRT_3	treated	/absolute/path/TRT_3.bam
```

Accepted alignment path column names include:

```text
bam
file
path
```

Relative paths are resolved relative to the sample sheet.

`run-3seq` performs:

- protocol-aware endpoint extraction;
- PAS clustering;
- PAS annotation;
- PAS counting;
- poly(A)-site usage estimation;
- TUTR/ALE/IPA interpretation;
- PAS-level differential testing;
- gene-level compositional testing.

### Gene-level inference

For replicated genes with at least two retained PAS and at least two evaluable samples per condition, the primary gene-level test is a **moderated Dirichlet-multinomial likelihood-ratio test**.

For a gene with `K` retained PAS, biological replicate `s` contributes:

```text
C_gs = (C_g1s, C_g2s, ..., C_gKs)
```

APAURA compares:

```text
H0: control and case share the same PAS composition
H1: control and case have different PAS compositions
```

using:

```text
LR = 2 × (logLik_full - logLik_null)
```

Gene-level p-values are adjusted using Benjamini-Hochberg FDR.

In canonical `apa.tsv`:

```text
effect_size = APD
pvalue = moderated Dirichlet-multinomial gene-level p-value
fdr = BH-adjusted gene-level DM FDR
```

`apa_sites.tsv` remains PAS-level and contains PAU, ΔPAU, and PAS-level statistics.

### Current `run-3seq` boundaries

Protocol presets currently primarily encode **protocol-specific endpoint geometry**.

Protocol-specific UMI parsing and FASTA-based internal-priming filtering are not yet implemented.

Not every named protocol preset has been independently benchmarked. Assay-specific validation is therefore recommended when applying APAURA to a new protocol.

See:

[`docs/TUTORIAL_3SEQ.md`](docs/TUTORIAL_3SEQ.md)

---

## 2. Direct RNA / true 3'-end analysis

Use `run-direct` when aligned read ends represent true or direct transcript 3' ends.

Example:

```bash
apaura run-direct --samples samples.tsv --gtf annotation.gtf --control control --case treated --out direct_results
```

To preserve exact single-nucleotide cleavage coordinates:

```bash
apaura run-direct --samples samples.tsv --gtf annotation.gtf --control control --case treated --cluster-distance 0 --out direct_results
```

---

## 3. IsoQuant long-read analysis

APAURA can analyze IsoQuant corrected-read output.

Example sample sheet:

```text
sample	condition	isoquant_dir
CTL_1	control	/absolute/path/CTL_1_isoquant
CTL_2	control	/absolute/path/CTL_2_isoquant
TRT_1	treated	/absolute/path/TRT_1_isoquant
TRT_2	treated	/absolute/path/TRT_2_isoquant
```

Run:

```bash
apaura run-isoquant --samples isoquant_samples.tsv --gtf annotation.gtf --control control --case treated --out isoquant_results
```

---

## 4. Conventional short-read RNA-seq

Use `short-read` for conventional fragmented RNA-seq.

Example:

```bash
apaura short-read --samples samples.tsv --gtf annotation.gtf --control control --case treated --library reverse --out short_read_results
```

Valid library modes are:

```text
unstranded
forward
reverse
```

Short-read fragment ends are **not interpreted directly as cleavage sites**.

APAURA instead uses:

- transcript-aware annotation;
- CIGAR-aware coverage;
- annotated terminal PAS structure;
- terminal exon group usage;
- intronic termination evidence;
- TUTR, ALE, and IPA-specific models.

---

# Single-cell APA

APAURA provides a platform-agnostic single-cell APA workflow through:

```bash
apaura sc --help
```

The workflow is organized as:

```text
annotation or external PAS calls
              ↓
      gtf-sites / import-sites
              ↓
          merge-sites
              ↓
           quantify
              ↓
          reference
              ↓
           compare
              ↓
cell-context-aware APA results
```

| Command | Purpose |
|---|---|
| `sc gtf-sites` | Build a strand-aware candidate PAS catalog from GTF/GFF |
| `sc import-sites` | Import PAS calls generated externally |
| `sc merge-sites` | Merge PAS catalogs into a consensus reference |
| `sc quantify` | Assign single-cell molecules to PAS and perform UMI-aware quantification |
| `sc reference` | Build grouped PAS-usage reference tables |
| `sc compare` | Compare PAS usage between biological conditions |

## Biological replication principle

APAURA keeps **cellular resolution** separate from **biological replication**.

The hierarchy is:

```text
molecule
   ↓
cell
   ↓
cell group / cell type
   ↓
biological sample
   ↓
condition
```

Cells from the same biological sample are not treated as independent biological replicates.

For grouped analyses, PAS counts are summarized within biological sample before condition-level inference.

When adequate biological replication is unavailable, APAURA reports descriptive `OBSERVED_*` remodeling classes rather than treating cells as independent biological replicates.

---

## Build a PAS catalog

```bash
apaura sc gtf-sites --annotation annotation.gtf --out sites.tsv
```

The default clustering distance is 24 nt.

A custom value can be specified:

```bash
apaura sc gtf-sites --annotation annotation.gtf --cluster-bp 24 --out sites.tsv
```

Typical outputs include:

```text
sites.tsv
sites.site_architecture.tsv
sites.gene_architecture.tsv
sites.architecture_summary.tsv
```

Single-cell architecture classes include:

```text
SINGLE_PAS
TUTR
ALE
MIXED
UNRESOLVED
```

---

## Import external PAS calls

```bash
apaura sc import-sites --help
```

This permits experimentally derived or externally detected PAS catalogs to be incorporated into APAURA.

---

## Merge PAS catalogs

```bash
apaura sc merge-sites --help
```

Possible PAS sources include:

```text
annotation-derived transcript ends
+
3'-end sequencing
+
direct RNA sequencing
+
external APA callers
```

---

## Quantify PAS usage

Example:

```bash
apaura sc quantify --sites sites.tsv --alignment sample.bam --sample-id SAMPLE_1 --gene-tag GX --three-prime-by gene --site-window 50 --out SAMPLE_1.counts.tsv
```

The core count representation is:

```text
sample_id	barcode	site_id	count
```

For example:

```text
sample_id	barcode	site_id	count
SAMPLE_1	CELL001	GENE1_PAS1	1
SAMPLE_1	CELL001	GENE1_PAS2	3
SAMPLE_1	CELL002	GENE1_PAS1	2
```

Molecules are UMI-deduplicated before PAS counts are written.

---

## Build grouped reference tables

Example cell metadata:

```text
barcode	sample_id	cell_type
CELL001	SAMPLE_1	B_cell
CELL002	SAMPLE_1	B_cell
CELL003	SAMPLE_1	T_cell
```

Run:

```bash
apaura sc reference --sites sites.tsv --counts counts.tsv --cells cells.tsv --group-by cell_type --out reference_results
```

Multiple grouping variables can be used:

```bash
apaura sc reference --sites sites.tsv --counts counts.tsv --cells cells.tsv --group-by cell_type,cell_subtype --out reference_results
```

The conceptual aggregation is:

```text
sample × cell group × gene × PAS
```

---

## Compare biological conditions

Example biological sample table:

```text
sample_id	condition
CTRL_1	control
CTRL_2	control
TRT_1	treated
TRT_2	treated
```

Run:

```bash
apaura sc compare --sites sites.tsv --samples samples.tsv --counts counts.tsv --control control --treated treated --out comparison_results
```

Cell-context-aware comparison:

```bash
apaura sc compare --sites sites.tsv --samples samples.tsv --counts counts.tsv --cells cells.tsv --group-by cell_type --control control --treated treated --out comparison_results
```

Important options include:

```text
--min-sample-count
--min-condition-count
--min-delta
--min-jsd
--fdr
--distribution-permutations
--min-interaction
--no-interactions
```

### Single-cell remodeling metrics

APAURA summarizes:

```text
UTR index
    directional proximal ↔ distal movement

entropy
    concentration or dispersion of PAS usage

Jensen-Shannon divergence
    redistribution of the complete PAS-usage distribution

dominant PAS switch
    change in the most frequently used PAS
```

### Unreplicated classifications

Possible descriptive classes include:

```text
OBSERVED_SHORTENING
OBSERVED_LENGTHENING
OBSERVED_COMPLEX_SWITCH
NO_MAJOR_OBSERVED_SHIFT
LOW_COVERAGE
```

`OBSERVED_*` indicates that the classification describes observed samples and is not presented as replicate-supported population inference.

---

# RBP motif enrichment

APAURA can test sequence motifs around regulated PAS and ALE splice-acceptor regions.

Example:

```bash
apaura rbp-enrich --sites results/apa_sites.tsv --genes results/apa.tsv --fasta genome.fa --motifs rbp_motifs.meme --out rbp_results
```

RBP enrichment identifies **candidate regulatory associations**.

It does not by itself establish direct RBP binding or causality.

Typical outputs include:

```text
rbp.tsv
rbp_hits.tsv
rbp_profile.tsv
rbp_windows.tsv
rbp_run_info.txt
```

---

# Locus visualization

`plot-region` is APAURA's publication-oriented locus renderer.

It can operate using completed APAURA results or directly from BAM/SAM alignments and a GTF.

## Plot completed APAURA results

```bash
apaura plot-region --results results --gene HNRNPA2B1 --samples samples.tsv --gtf annotation.gtf --control U1_WT --case U1_MUT --library reverse --show-replicates --show-overview --format pdf --out plots/HNRNPA2B1
```

Depending on the analysis, figures may include:

- whole-gene coverage;
- terminal-region coverage;
- biological replicate tracks;
- tested PAS;
- ΔPAU;
- TUTR direction;
- IPA usage;
- ALE changes;
- terminal APA architecture;
- annotated terminal-exon structure.

## Standalone BAM/SAM + GTF visualization

A completed APAURA analysis is not required.

```bash
apaura plot-region --samples samples.tsv --gtf annotation.gtf --gene HNRNPA2B1 --control U1_WT --case U1_MUT --library reverse --show-replicates --show-overview --format pdf --out plots/HNRNPA2B1
```

Standalone mode does **not** perform differential APA testing.

Without completed APAURA results, formal ΔPAU, TUTR, ALE, IPA, and APA architecture calls are not inferred from the plot alone.

## Multiple genes

```bash
apaura plot-region --samples samples.tsv --gtf annotation.gtf --gene HNRNPA2B1,PLA2G12A,TMEM254 --control U1_WT --case U1_MUT --library reverse --show-overview --format pdf --out plots
```

Using a gene list:

```bash
apaura plot-region --samples samples.tsv --gtf annotation.gtf --genes-file genes.txt --control U1_WT --case U1_MUT --library reverse --out plots
```

## Genomic region

```bash
apaura plot-region --samples samples.tsv --gtf annotation.gtf --region chr7:26170000-26197000 --control U1_WT --case U1_MUT --library reverse --format pdf --out plots/chr7_region
```

### Optional plotting dependencies

`plot-region` additionally requires:

- `Rscript`;
- `samtools`;
- `ggplot2`;
- `dplyr`;
- `readr`;
- `stringr`;
- `purrr`;
- `tibble`;
- `patchwork`;
- `rtracklayer`;
- `GenomicRanges`.

Install R dependencies with:

```r
install.packages(c("ggplot2", "dplyr", "readr", "stringr", "purrr", "tibble", "patchwork"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("rtracklayer", "GenomicRanges"))
```

---

# Primary outputs

## Canonical APA output

APAURA analyses use a canonical output structure centered on:

```text
apa.tsv
apa_events.tsv
apa_sites.tsv
```

### `apa.tsv`

Gene-level APA summary containing overall architecture, direction, effect sizes, statistical evidence, and mechanism-level summaries.

### `apa_events.tsv`

Event-level TUTR, ALE, and IPA interpretation.

### `apa_sites.tsv`

PAS-level usage and differential evidence.

Additional detailed and diagnostic outputs may be written under:

```text
details/
```

---

## Direct RNA / 3'-end / IsoQuant supporting outputs

Depending on mode, supporting files can include:

```text
pas_all.bed
pas.bed
pas_counts.tsv
pau.tsv
apa_architecture.tsv
run_info.txt
```

See:

[`docs/OUTPUT_SCHEMA.md`](docs/OUTPUT_SCHEMA.md)

---

# Core quantities

For PAS `k` in gene `g`, sample `s`:

```text
PAU_gks = C_gks / Σ_j C_gjs
```

PAS-level differential usage is:

```text
delta_PAU = mean(PAU_case) - mean(PAU_control)
```

Within a TUTR-capable terminal branch, PAS are ordered proximal-to-distal in transcript orientation and assigned normalized positions `x_k`.

The weighted positional index is:

```text
TUTR_WPI = Σ_k PAU_k × x_k
```

Therefore:

```text
delta_tutr_wpi > 0 → distal movement / lengthening
delta_tutr_wpi < 0 → proximal movement / shortening
```

Usage redistribution can also be summarized using total-variation distance:

```text
TVD = 0.5 × Σ_k |PAU_case,k - PAU_control,k|
```

ALE and IPA are quantified independently so branch switching and intronic termination are not collapsed into a TUTR positional statistic.

---

# Statistical methods

APAURA uses mechanism-appropriate statistical models and effect-size summaries.

Methods include:

- moderated Dirichlet-multinomial likelihood-ratio testing for replicated `run-3seq` gene-level PAS composition;
- beta-binomial likelihood-ratio testing for replicated PAS-usage comparisons;
- Fisher's exact testing only for genuine experimental 1-vs-1 designs;
- Benjamini-Hochberg multiple-testing correction;
- weighted positional indices for TUTR and ALE;
- total-variation distance;
- Jensen-Shannon divergence;
- weighted isotonic regression for short-read terminal-coverage deconvolution;
- tie-aware distributional statistics where appropriate.

Replicated designs lacking sufficient evaluable replicate support remain `not_tested` rather than being converted to a 1-vs-1 fallback test.

Statistical significance should be interpreted together with effect size and biological direction.

---

# `run-3seq` output contract

For replicated `run-3seq` analyses, the primary hierarchy is:

```text
apa.tsv
    gene-level APA summary

apa_events.tsv
    TUTR/ALE/IPA event interpretation

apa_sites.tsv
    PAS-level PAU, delta_PAU, and site-level statistics

details/three_seq_dm.tsv
    moderated Dirichlet-multinomial diagnostics
```

For evaluable multi-PAS genes, canonical `apa.tsv` reports:

```text
effect_size = APD
pvalue = moderated Dirichlet-multinomial gene-level p-value
fdr = BH-adjusted gene-level DM FDR
```

Use `apa_sites.tsv` to identify the individual PAS contributing to a gene-level compositional change.

---

# Testing and validation

Run Java unit tests with:

```bash
mvn test
```

Run the synthetic regression suite with:

```bash
bash tests/run_tests.sh
```

APAURA v0.10.1 includes synthetic regression coverage for:

- SAM and BAM reading;
- true 3'-end/direct-RNA analysis;
- protocol-aware `run-3seq`;
- IsoQuant input;
- short-read TUTR;
- ALE;
- IPA;
- mixed TUTR + ALE architecture;
- multi-mechanism TUTR + IPA changes;
- no-directional-shift cases;
- RBP motif enrichment;
- single-cell analysis;
- statistical utilities.

For the validated release:

```text
All APAURA v0.10.1 synthetic tests passed.
```

Synthetic tests verify implementation against constructed truth. They do not independently establish sensitivity, specificity, or calibration across every biological dataset or sequencing protocol.

Real-data benchmarking remains important.

---

# Current boundaries

APAURA v0.10.1 has several intentional boundaries:

- conventional fragmented RNA-seq endpoints are not treated as direct cleavage-site observations;
- short-read analysis is annotation-guided rather than a general de novo PAS breakpoint-discovery engine;
- `run-3seq` does not currently perform protocol-specific UMI parsing;
- `run-3seq` does not currently provide a complete FASTA-based internal-priming filtering framework;
- not every named 3'-end protocol preset has been independently benchmarked;
- standalone `plot-region` is a visualization mode and does not perform differential APA testing;
- standalone locus plotting requires GTF annotation to resolve transcript structure;
- RBP motif enrichment reports candidate sequence associations rather than direct binding or causality;
- the built-in sequential BAM reader implements the BAM fields required by APAURA rather than providing a complete random-access BAM API.

---

# Documentation

Detailed documentation is available in:

- [`docs/INSTALLATION.md`](docs/INSTALLATION.md)
- [`docs/QUICKSTART.md`](docs/QUICKSTART.md)
- [`docs/ALGORITHM.md`](docs/ALGORITHM.md)
- [`docs/ARCHITECTURE_MODEL.md`](docs/ARCHITECTURE_MODEL.md)
- [`docs/OUTPUT_SCHEMA.md`](docs/OUTPUT_SCHEMA.md)
- [`docs/STRANDEDNESS.md`](docs/STRANDEDNESS.md)
- [`docs/VALIDATION.md`](docs/VALIDATION.md)
- [`docs/TROUBLESHOOTING.md`](docs/TROUBLESHOOTING.md)
- [`docs/DEVELOPMENT.md`](docs/DEVELOPMENT.md)
- [`docs/TUTORIAL_3SEQ.md`](docs/TUTORIAL_3SEQ.md)
- [`docs/TUTORIAL_SHORT_READ.md`](docs/TUTORIAL_SHORT_READ.md)
- [`docs/TUTORIAL_DIRECT_RNA.md`](docs/TUTORIAL_DIRECT_RNA.md)
- [`docs/TUTORIAL_ISOQUANT.md`](docs/TUTORIAL_ISOQUANT.md)
- [`docs/TUTORIAL_PLOT_REGION.md`](docs/TUTORIAL_PLOT_REGION.md)
- [`docs/TUTORIAL_RBP_ENRICHMENT.md`](docs/TUTORIAL_RBP_ENRICHMENT.md)
- [`docs/TUTORIAL_BATCH_COMPARISONS.md`](docs/TUTORIAL_BATCH_COMPARISONS.md)
- [`docs/RELEASE_CHECKLIST.md`](docs/RELEASE_CHECKLIST.md)

---

# Citation

Print the citation with:

```bash
apaura citation
```

Please cite:

> **Boateng-Sarfo G, E Niewold, Sigor S, Wang L. APAURA: Alternative Polyadenylation Architecture, Usage, and Regulatory Analysis. Version 0.10.1, 2026.**

A machine-readable citation is provided in:

[`CITATION.cff`](CITATION.cff)

---

# License

APAURA is distributed under the MIT License.

See:

[`LICENSE`](LICENSE)
