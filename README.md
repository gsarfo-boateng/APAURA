# APAURA

**Alternative Polyadenylation Architecture, Usage, and Regulatory Analysis**

[![Version](https://img.shields.io/badge/version-0.10.1-blue.svg)](#citation)
[![Java](https://img.shields.io/badge/Java-8%2B-orange.svg)](#installation)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

APAURA is a Java toolkit for detecting, quantifying, classifying, and visualizing alternative polyadenylation (APA) from RNA-sequencing data. It supports true 3′-end/direct-RNA evidence, IsoQuant corrected reads, conventional short-read RNA-seq, RNA-binding-protein motif enrichment, and publication-style locus visualization.

**Version:** `0.10.1`

> If you use APAURA, please cite:
>
> **Boateng-Sarfo G. APAURA: Alternative Polyadenylation Architecture, Usage, and Regulatory Analysis. Version 0.10.1, 2026.**

---

## Overview

APAURA follows a layered biological model:

```text
poly(A) sites (PAS)
        ↓
poly(A)-site usage (PAU)
        ↓
terminal APA architecture
        ↓
mechanism-specific differential APA
        ↓
candidate RBP regulators
```

APAURA separates **terminal APA architecture** from **the mechanism that changes between conditions**. This allows tandem 3′ UTR APA (TUTR), alternative last exon usage (ALE), and intronic polyadenylation (IPA) to be represented without collapsing distinct biological processes into one proximal/distal score.

### Terminal APA architecture

| Architecture | Interpretation |
|---|---|
| `TUTR` | At least one terminal-exon branch contains multiple PAS |
| `ALE` | Multiple alternative terminal-exon branches are present |
| `mixed` | Both TUTR and ALE architecture are present |
| `none` | No retained TUTR/ALE terminal architecture is identified |

**IPA is not an `apa_architecture` value.** Intronic polyadenylation is quantified as a separate termination mechanism and can co-occur with TUTR or ALE.

### Direction of change

```text
TUTR → shortening / lengthening / no_directional_shift
ALE  → proximal_switch / distal_switch / switch
IPA  → gain / loss
multiple changing mechanisms → complex
```

Mechanism-specific columns such as `tutr_direction`, `ale_direction`, and `ipa_direction` preserve simultaneous effects within the same gene.

---

## Main capabilities

| Task | Input | APAURA command |
|---|---|---|
| True 3′-end / direct RNA APA | BAM/SAM + GTF | `run-direct` |
| IsoQuant long-read APA | IsoQuant corrected reads + GTF | `run-isoquant` |
| Conventional RNA-seq APA | BAM/SAM + GTF | `short-read` |
| Publication locus plots | APAURA results or BAM/SAM + GTF | `plot-region` |
| RBP motif enrichment | APAURA APA results + FASTA + MEME motifs | `rbp-enrich` |
| Transcript-end extraction | BAM/SAM | `ends` |

---

## Installation

### Core requirements

- Java 8 or newer
- Maven for the standard source build

For `plot-region` you also need:

- `Rscript`
- `samtools`
- R packages: `ggplot2`, `dplyr`, `readr`, `stringr`, `purrr`, `tibble`, `patchwork`
- Bioconductor packages: `rtracklayer`, `GenomicRanges`

Install plotting dependencies in R:

```r
install.packages(c(
  "ggplot2", "dplyr", "readr", "stringr",
  "purrr", "tibble", "patchwork"
))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("rtracklayer", "GenomicRanges"))
```

### Build from source

```bash
git clone https://github.com/gsarfo-boateng/APAURA.git
cd APAURA

mvn clean package
./build.sh
```

The portable build creates:

```text
APAURA.jar
```

Check the version:

```bash
java -jar APAURA.jar --version
```

Expected:

```text
APAURA v0.10.1
Alternative Polyadenylation Architecture, Usage, and Regulatory Analysis
```

---

## Quick start

### Conventional short-read RNA-seq

```bash
java -jar APAURA.jar short-read \
  --samples samples.tsv \
  --gtf annotation.gtf \
  --control control \
  --case treated \
  --library reverse \
  --out results
```

Short-read fragment ends are **not** interpreted as cleavage sites. APAURA instead uses transcript-aware annotation, CIGAR-aware coverage, annotated terminal PAS structure, terminal-exon branch usage, and intronic termination evidence to quantify TUTR, ALE, and IPA.

Valid library modes are:

```text
unstranded
forward
reverse
```

### Direct RNA / true 3′-end evidence

```bash
java -jar APAURA.jar run-direct \
  --samples samples.tsv \
  --gtf annotation.gtf \
  --control control \
  --case treated \
  --out results
```

To preserve exact single-nucleotide cleavage coordinates:

```bash
--cluster-distance 0
```

### IsoQuant corrected reads

```bash
java -jar APAURA.jar run-isoquant \
  --samples isoquant_samples.tsv \
  --gtf annotation.gtf \
  --control control \
  --case treated \
  --out results
```

### RBP motif enrichment

```bash
java -jar APAURA.jar rbp-enrich \
  --sites results/apa_sites.tsv \
  --genes results/apa.tsv \
  --fasta genome.fa \
  --motifs rbp_motifs.meme \
  --out rbp_results
```

The RBP module tests sequence motifs around regulated PAS and, for ALE events, terminal-exon splice-acceptor windows. Enrichment identifies **candidate regulatory associations**; it is not evidence of direct binding or causality by itself.

---

## Sample sheets

### BAM/SAM input

```text
sample	condition	bam
CTL_1	control	/path/CTL_1.bam
CTL_2	control	/path/CTL_2.bam
TRT_1	treated	/path/TRT_1.bam
TRT_2	treated	/path/TRT_2.bam
```

Accepted alignment path columns include `bam`, `file`, and `path`. Relative paths are resolved relative to the sample sheet.

### IsoQuant input

```text
sample	condition	isoquant_dir
CTL_1	control	/path/CTL_1_isoquant
CTL_2	control	/path/CTL_2_isoquant
TRT_1	treated	/path/TRT_1_isoquant
TRT_2	treated	/path/TRT_2_isoquant
```

---

# `plot-region`

`plot-region` is APAURA's publication-oriented locus renderer. It can be used in two distinct modes.

## A. Plot a completed APAURA analysis

When `--results` is supplied, the plot can combine normalized coverage with tested APA results.

```bash
java -jar APAURA.jar plot-region \
  --results results \
  --gene HNRNPA2B1 \
  --samples samples.tsv \
  --gtf annotation.gtf \
  --control U1_WT \
  --case U1_MUT \
  --library reverse \
  --show-replicates \
  --show-overview \
  --format pdf \
  --out plots/HNRNPA2B1
```

Depending on the locus and analysis mode, the figure can include:

- whole-gene coverage overview;
- zoomed 3′ terminal-region coverage;
- replicate coverage;
- tested terminal PAS and ΔPAU;
- TUTR direction;
- IPA usage;
- ALE branch-usage changes;
- terminal APA architecture;
- annotated last-exon structure.

## B. Standalone BAM/SAM + GTF visualization

A completed APAURA result directory is **not required**.

```bash
java -jar APAURA.jar plot-region \
  --samples samples.tsv \
  --gtf annotation.gtf \
  --gene HNRNPA2B1 \
  --control U1_WT \
  --case U1_MUT \
  --library reverse \
  --show-replicates \
  --show-overview \
  --format pdf \
  --out plots/HNRNPA2B1
```

Standalone mode shows:

- normalized BAM/SAM coverage;
- control-versus-case coverage;
- optional replicate tracks;
- whole-gene overview;
- annotation-derived terminal-region zoom;
- GTF-derived alternative last-exon structure.

Standalone mode **does not perform differential APA testing**. Therefore it does not report formal ΔPAU, TUTR, ALE, IPA, or APA-architecture calls when no APAURA results are supplied.

### Multiple genes

```bash
java -jar APAURA.jar plot-region \
  --samples samples.tsv \
  --gtf annotation.gtf \
  --gene HNRNPA2B1,PLA2G12A,TMEM254 \
  --control U1_WT \
  --case U1_MUT \
  --library reverse \
  --show-overview \
  --format pdf \
  --out plots
```

Or:

```bash
java -jar APAURA.jar plot-region \
  --samples samples.tsv \
  --gtf annotation.gtf \
  --genes-file genes.txt \
  --control U1_WT \
  --case U1_MUT \
  --library reverse \
  --out plots
```

### Genomic region

```bash
java -jar APAURA.jar plot-region \
  --samples samples.tsv \
  --gtf annotation.gtf \
  --region chr7:26170000-26197000 \
  --control U1_WT \
  --case U1_MUT \
  --library reverse \
  --format pdf \
  --out plots/chr7_region
```

See [`docs/TUTORIAL_PLOT_REGION.md`](docs/TUTORIAL_PLOT_REGION.md) for details.

---

## Primary outputs

### Conventional short-read analysis

```text
apa.tsv
apa_events.tsv
apa_sites.tsv
```

- `apa.tsv` — gene-level APA summary and mechanism-specific direction
- `apa_events.tsv` — event-level TUTR/ALE/IPA results
- `apa_sites.tsv` — PAS/site-level differential results

Additional diagnostic/intermediate tables may be produced for detailed inspection and validation.

### Direct RNA / true 3′-end and IsoQuant analyses

Common outputs include:

```text
pas_all.bed
pas.bed
pas_counts.tsv
pau.tsv
apa_sites.tsv
apa.tsv
apa_architecture.tsv
run_info.txt
```

### RBP enrichment

```text
rbp.tsv
rbp_hits.tsv
rbp_profile.tsv
rbp_windows.tsv
rbp_run_info.txt
```

See [`docs/OUTPUT_SCHEMA.md`](docs/OUTPUT_SCHEMA.md) for column definitions.

---

## Core quantities

For PAS `k` in gene `g`, sample `s`:

```text
PAU_gks = C_gks / sum_j C_gjs
```

Site-level change:

```text
delta_PAU = mean(PAU_case) - mean(PAU_control)
```

Within a TUTR-capable terminal branch, sites are ordered proximal-to-distal in transcript orientation and assigned normalized positions `x_k`:

```text
TUTR_WPI = sum_k PAU_k * x_k
```

Therefore:

```text
delta_tutr_wpi > 0  → distal movement / lengthening
delta_tutr_wpi < 0  → proximal movement / shortening
```

Usage redistribution can also be summarized by total-variation distance:

```text
TVD = 0.5 * sum_k |PAU_case,k - PAU_control,k|
```

ALE and IPA are quantified independently so branch switching and intronic termination do not get collapsed into a TUTR positional statistic.

---

## Statistical methods

APAURA uses mechanism-appropriate statistical and effect-size summaries, including:

- beta-binomial likelihood-ratio testing for replicated PAS-usage comparisons;
- Fisher's exact test for the 1-vs-1 fallback case;
- Benjamini-Hochberg FDR correction;
- weighted positional indices for TUTR and ALE direction;
- total-variation distance and Jensen-Shannon divergence for distributional shifts;
- weighted isotonic regression for short-read terminal-coverage deconvolution;
- tie-aware distributional statistics where appropriate.

Interpret significance together with effect size and mechanism-specific direction.

---

## Testing and validation

Run Java unit tests with:

```bash
mvn test
```

Run the complete synthetic regression suite with:

```bash
bash tests/run_tests.sh
```

APAURA v0.10.1 includes synthetic regression coverage for:

- SAM and BAM reading;
- true 3′-end/direct-RNA analysis;
- IsoQuant input;
- short-read TUTR;
- ALE;
- IPA;
- mixed TUTR + ALE architecture;
- multi-mechanism TUTR + IPA changes;
- no-directional-shift edge cases;
- RBP motif enrichment.

For the current validated release:

```text
All APAURA v0.10.1 synthetic tests passed.
```

Synthetic testing verifies implementation against known constructed truth. It does not by itself establish sensitivity, specificity, or calibration across all biological datasets. Real-data benchmarking remains important before v1.0.

---

## Current boundaries

APAURA v0.10.1 has several intentional boundaries:

- Conventional fragmented RNA-seq endpoints are not treated as direct cleavage-site observations.
- Short-read analysis is annotation-guided rather than a general de novo PAS breakpoint-discovery engine.
- Standalone `plot-region` is a visualization mode; without completed APAURA results it does not perform differential APA testing.
- Standalone locus plotting currently uses a GTF annotation to resolve gene/terminal structure.
- RBP motif enrichment reports candidate sequence associations, not direct RBP binding or causality.
- The built-in sequential BAM reader implements the fields APAURA requires rather than the complete BAM optional-tag/random-access API.
- Sequence-based internal-priming filtering is not yet a complete general-purpose module.

---

## Help

```bash
java -jar APAURA.jar --help

java -jar APAURA.jar run-direct --help
java -jar APAURA.jar run-isoquant --help
java -jar APAURA.jar short-read --help
java -jar APAURA.jar plot-region --help
java -jar APAURA.jar rbp-enrich --help
java -jar APAURA.jar ends --help
```

---

## Documentation

- [`docs/INSTALLATION.md`](docs/INSTALLATION.md)
- [`docs/QUICKSTART.md`](docs/QUICKSTART.md)
- [`docs/ALGORITHM.md`](docs/ALGORITHM.md)
- [`docs/ARCHITECTURE_MODEL.md`](docs/ARCHITECTURE_MODEL.md)
- [`docs/OUTPUT_SCHEMA.md`](docs/OUTPUT_SCHEMA.md)
- [`docs/STRANDEDNESS.md`](docs/STRANDEDNESS.md)
- [`docs/VALIDATION.md`](docs/VALIDATION.md)
- [`docs/TROUBLESHOOTING.md`](docs/TROUBLESHOOTING.md)
- [`docs/DEVELOPMENT.md`](docs/DEVELOPMENT.md)
- [`docs/TUTORIAL_SHORT_READ.md`](docs/TUTORIAL_SHORT_READ.md)
- [`docs/TUTORIAL_DIRECT_RNA.md`](docs/TUTORIAL_DIRECT_RNA.md)
- [`docs/TUTORIAL_ISOQUANT.md`](docs/TUTORIAL_ISOQUANT.md)
- [`docs/TUTORIAL_PLOT_REGION.md`](docs/TUTORIAL_PLOT_REGION.md)
- [`docs/TUTORIAL_RBP_ENRICHMENT.md`](docs/TUTORIAL_RBP_ENRICHMENT.md)
- [`docs/TUTORIAL_BATCH_COMPARISONS.md`](docs/TUTORIAL_BATCH_COMPARISONS.md)
- [`docs/RELEASE_CHECKLIST.md`](docs/RELEASE_CHECKLIST.md)

---

## Citation

```bash
java -jar APAURA.jar --citation
```

Please cite:

> **Boateng-Sarfo G. APAURA: Alternative Polyadenylation Architecture, Usage, and Regulatory Analysis. Version 0.10.1, 2026.**

A machine-readable citation is provided in [`CITATION.cff`](CITATION.cff).

---

## License

APAURA is distributed under the MIT License. See [`LICENSE`](LICENSE).
