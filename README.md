# APAURA

**Alternative Polyadenylation Architecture, Usage, and Regulatory Analysis**

[![Version](https://img.shields.io/badge/version-0.10.1-blue.svg)](#citation)
[![Java](https://img.shields.io/badge/Java-8%2B-orange.svg)](#installation)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

APAURA is a Java toolkit for detecting, quantifying, classifying, and visualizing alternative polyadenylation (APA) from RNA-sequencing data. It supports true 3'-end/direct-RNA evidence, IsoQuant corrected reads, conventional short-read RNA-seq, single-cell APA analysis, RNA-binding-protein motif enrichment, and publication-style locus visualization.

                                              apaura short-read ...
                                              apaura run-direct ...
                                              apaura run-isoquant ...
                                              apaura sc ...
                                              apaura rbp-enrich ...
                                              apaura plot-region ...
                                              
                                      short-read --------\\
                                      run-direct --------//
                                                             >---- APAURA ----> TUTR / ALE / MIXED/ IPA
                                      run-isoquant --------\\
                                       sc          --------//
                                      


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

APAURA separates **terminal APA architecture** from **the mechanism that changes between conditions**. This allows tandem 3' UTR APA (TUTR), alternative last exon usage (ALE), and intronic polyadenylation (IPA) to be represented without collapsing distinct biological processes into one proximal/distal score.

### Terminal APA architecture

| Architecture | Interpretation |
|---|---|
| `SINGLE_PAS` | One resolved terminal exon group with one retained terminal PAS |
| `TUTR` | One resolved terminal exon group with multiple retained PAS |
| `ALE` | Multiple resolved terminal exon groups, without tandem PAS within a group |
| `mixed` | Multiple resolved terminal exon groups with tandem PAS in at least one group |
| `UNRESOLVED` | Retained terminal PAS cannot be assigned completely to reconciled terminal structure |

**IPA is not an `apa_architecture` value.** Intronic polyadenylation is quantified as a separate termination mechanism. A gene can therefore have local IPA evidence independently of its terminal architecture, including when its overall terminal architecture is `UNRESOLVED`.

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
| True 3'-end / direct RNA APA | BAM/SAM + GTF | `run-direct` |
| IsoQuant long-read APA | IsoQuant corrected reads + GTF | `run-isoquant` |
| Conventional RNA-seq APA | BAM/SAM + GTF | `short-read` |
| Single-cell APA | Single-cell alignments + PAS reference + cell metadata | `sc` |
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
install.packages(c("ggplot2", "dplyr", "readr", "stringr", "purrr", "tibble", "patchwork"))

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
apaura --version
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
apaura short-read --samples samples.tsv  --gtf annotation.gtf --control control --case treated --library reverse --out results
```

Short-read fragment ends are **not** interpreted as cleavage sites. APAURA instead uses transcript-aware annotation, CIGAR-aware coverage, annotated terminal PAS structure, terminal exon group usage, and intronic termination evidence to quantify TUTR, ALE, and IPA.

Valid library modes are:

```text
unstranded
forward
reverse
```

### Direct RNA / true 3'-end evidence

```bash
apaura run-direct --samples samples.tsv --gtf annotation.gtf --control control --case treated --out results
```

To preserve exact single-nucleotide cleavage coordinates:

```bash
--cluster-distance 0
```

### IsoQuant corrected reads

```bash
apaura run-isoquant --samples isoquant_samples.tsv --gtf annotation.gtf --control control --case treated --out results
```
### Single-cell alternative polyadenylation

APAURA provides a dedicated `sc` mode for single-cell alternative polyadenylation analysis.

```bash
apaura sc --help
```

The single-cell workflow is organized into six commands:

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
cell-context-aware, architecture-aware APA results
```

| Command | Purpose |
|---|---|
| `sc gtf-sites` | Build a strand-aware candidate poly(A) site catalog from GTF/GFF annotation |
| `sc import-sites` | Import PAS calls generated by an external APA caller |
| `sc merge-sites` | Merge multiple PAS catalogs into a consensus reference |
| `sc quantify` | Assign single-cell molecules to PAS and perform UMI-aware quantification |
| `sc reference` | Build grouped PAS-usage reference tables |
| `sc compare` | Compare PAS usage between biological conditions |

### Biological replicate principle

APAURA keeps **cellular resolution** separate from **biological replication**.

Cells contribute molecular depth and allow APA usage to be evaluated within cell populations, but cells from the same biological sample are not treated as independent biological replicates.

The analysis hierarchy is therefore:

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

For grouped analyses, PAS counts are summarized within biological sample before condition-level inference.

With sufficient biological replication, APAURA can perform replicate-supported condition comparisons and multiple-testing correction. When replication is not available, APAURA reports descriptive `OBSERVED_*` remodeling classes rather than treating individual cells as biological replicates.

### 1. Build a PAS catalog from GTF/GFF annotation

```bash
apaura sc gtf-sites --annotation annotation.gtf --out sites.tsv
```

Nearby annotated transcript ends are clustered into candidate PAS. The default clustering distance is 24 nt and can be changed with:

```bash
--cluster-bp 24
```

`gtf-sites` also generates transcript-end architecture annotations associated with the PAS catalog.

Typical outputs include:

```text
sites.tsv
sites.site_architecture.tsv
sites.gene_architecture.tsv
sites.architecture_summary.tsv
```

The single-cell architecture layer classifies genes as:

```text
SINGLE_PAS
TUTR
ALE
MIXED
UNRESOLVED
```

`TUTR` represents multiple PAS within a terminal exon group, `ALE` represents alternative terminal-exon usage, and `MIXED` indicates that both mechanisms are present.

### 2. Import PAS calls from another APA tool

Existing PAS calls can be converted into the APAURA single-cell PAS format:

```bash
apaura sc import-sites --help
```

This allows experimentally derived or externally detected PAS catalogs to be used instead of, or together with, annotation-derived transcript ends.

### 3. Merge PAS catalogs

Multiple PAS catalogs can be combined into one consensus catalog:

```bash
apaura sc merge-sites --help
```

This is useful when PAS evidence comes from several sources, for example:

```text
annotation-derived transcript ends
+
3' end sequencing
+
direct RNA sequencing
+
external APA callers
```

The resulting consensus PAS catalog can then be used for molecule quantification.

### 4. Quantify PAS usage

Single-cell alignments are assigned to compatible PAS while retaining sample, cell-barcode, UMI, and gene information.

Example:

```bash
apaura sc quantify --sites sites.tsv --alignment sample.bam  --sample-id SAMPLE_1 --gene-tag GX --three-prime-by gene --site-window 50 --out SAMPLE_1.counts.tsv
```

The basic count representation is:

```text
sample_id    barcode    site_id    count
```

For example:

```text
sample_id   barcode   site_id       count
SAMPLE_1    CELL001   GENE1_PAS1    1
SAMPLE_1    CELL001   GENE1_PAS2    3
SAMPLE_1    CELL002   GENE1_PAS1    2
```

Molecules are UMI-deduplicated before PAS counts are written. Barcode, UMI, gene, strand, and transcript-oriented 3'-end information are retained during assignment so that repeated reads from the same molecule do not artificially increase PAS usage.

### 5. Build grouped reference tables

Cell metadata can be used to aggregate PAS usage within biological sample and cell population.

Example cell metadata:

```text
barcode    sample_id    cell_type
CELL001    SAMPLE_1     B_cell
CELL002    SAMPLE_1     B_cell
CELL003    SAMPLE_1     T_cell
```

Run:

```bash
apaura sc reference --sites sites.tsv --counts counts.tsv --cells cells.tsv --group-by cell_type --out reference_results
```

Multiple grouping variables can be combined:

```bash
--group-by cell_type,cell_subtype
```

The conceptual aggregation is:

```text
sample × cell group × gene × PAS
```

rather than treating individual cells as biological replicates.

When a matching `*.gene_architecture.tsv` file exists beside the PAS catalog, APAURA automatically propagates transcript-end architecture into compatible reference outputs.

### 6. Compare biological conditions

The `compare` command evaluates condition-dependent PAS redistribution.

A sample table contains biological sample identity and condition:

```text
sample_id    condition
CTRL_1       control
CTRL_2       control
TRT_1        treated
TRT_2        treated
```

Run:

```bash
apaura sc compare --sites sites.tsv --samples samples.tsv --counts counts.tsv --control control --treated treated --out comparison_results
```

For cell-context-aware comparisons, cell metadata and grouping variables can also be supplied:

```bash
apaura sc compare --sites sites.tsv --samples samples.tsv --counts counts.tsv --cells cells.tsv --group-by cell_type --control control --treated treated --out comparison_results
```

Important comparison options include:

```text
--min-sample-count INT
--min-condition-count INT
--min-delta FLOAT
--min-jsd FLOAT
--fdr FLOAT
--distribution-permutations INT
--min-interaction FLOAT
--no-interactions
```

### APA remodeling metrics

For each evaluable gene, APAURA reconstructs the PAS-usage distribution and summarizes several complementary properties of 3'-end remodeling.

```text
UTR index
    directional proximal ↔ distal movement

entropy
    concentration or dispersion of PAS usage

Jensen-Shannon divergence
    redistribution of the complete PAS-usage profile

dominant PAS switch
    change in the most frequently used PAS
```

These metrics are intentionally complementary. A gene can redistribute PAS usage substantially without producing a simple proximal-versus-distal shift.

### Remodeling classifications

For unreplicated comparisons, gene-level results can include:

```text
OBSERVED_SHORTENING
OBSERVED_LENGTHENING
OBSERVED_COMPLEX_SWITCH
NO_MAJOR_OBSERVED_SHIFT
LOW_COVERAGE
```

`OBSERVED_*` explicitly indicates that the classification describes the observed biological samples and is not being presented as replicate-supported population inference.

With adequate biological replication, statistical evidence and FDR can be evaluated at the biological-sample level.

### Cell-group interactions

When multiple cell groups are available, APAURA can test whether the treatment effect differs between groups.

Conceptually:

```text
interaction =
    treated(group1 - group2)
    -
    control(group1 - group2)
```

This asks whether condition-dependent APA remodeling is different between cell contexts rather than simply asking whether the two cell populations have different baseline PAS usage.

Interactions can be disabled with:

```bash
--no-interactions
```

### Single-cell output

Depending on the analysis and replication structure, `sc compare` can produce:

```text
all_results.tsv
unified_results.tsv
shortening.tsv
lengthening.tsv
complex_switches.tsv
observed_unreplicated_events.tsv
celltype_summary.tsv
celltype_enriched_events.tsv
classification_matrix.tsv
delta_utr_index_matrix.tsv
js_divergence_matrix.tsv
specificity_score_matrix.tsv
specificity_class_matrix.tsv
```

When transcript-end architecture is available, APAURA also produces architecture-aware copies such as:

```text
all_results.with_architecture.tsv
unified_results.with_architecture.tsv
observed_unreplicated_events.with_architecture.tsv
```

The architecture layer provides structural context for interpreting a PAS shift. It does **not** represent an additional treatment-response category.

For example, a proximal/distal redistribution may have different biological meaning depending on whether the underlying transcript structure is tandem 3' UTR (`TUTR`), alternative last exon (`ALE`), or `MIXED`.

For command-specific options:

```bash
apaura sc gtf-sites --help
apaura sc import-sites --help
apaura sc merge-sites --help
apaura sc quantify --help
apaura sc reference --help
apaura sc compare --help
```


### RBP motif enrichment

```bash
apaura rbp-enrich --sites results/apa_sites.tsv --genes results/apa.tsv --fasta genome.fa --motifs rbp_motifs.meme --out rbp_results
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
apaura plot-region --results results --gene HNRNPA2B1 --samples samples.tsv --gtf annotation.gtf --control U1_WT --case U1_MUT --library reverse --show-replicates --show-overview --format pdf --out plots/HNRNPA2B1
```

Depending on the locus and analysis mode, the figure can include:

- whole-gene coverage overview;
- zoomed 3' terminal-region coverage;
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
apaura plot-region --samples samples.tsv --gtf annotation.gtf --gene HNRNPA2B1 --control U1_WT --case U1_MUT --library reverse --show-replicates --show-overview --format pdf --out plots/HNRNPA2B1
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
apaura plot-region --samples samples.tsv --gtf annotation.gtf --gene HNRNPA2B1,PLA2G12A,TMEM254 --control U1_WT --case U1_MUT --library reverse --show-overview --format pdf --out plots
```

Or:

```bash
apaura plot-region --samples samples.tsv --gtf annotation.gtf --genes-file genes.txt --control U1_WT --case U1_MUT --library reverse --out plots
```

### Genomic region

```bash
apaura plot-region --samples samples.tsv --gtf annotation.gtf --region chr7:26170000-26197000 --control U1_WT --case U1_MUT --library reverse --format pdf --out plots/chr7_region
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

- `apa.tsv` - gene-level APA summary and mechanism-specific direction
- `apa_events.tsv` - event-level TUTR/ALE/IPA results
- `apa_sites.tsv` - PAS/site-level differential results

Additional diagnostic/intermediate tables may be produced for detailed inspection and validation.

### Direct RNA / true 3'-end and IsoQuant analyses

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
- Fisher's exact test only for genuine experimental 1-vs-1 designs; replicated designs with inadequate evaluable support remain `not_tested` with `insufficient_replicate_support`;
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
- true 3'-end/direct-RNA analysis;
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
apaura --help

apaura run-direct --help
apaura run-isoquant --help
apaura short-read --help
apaura plot-region --help
apaura rbp-enrich --help
apaura ends --help
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
apaura --citation
```

Please cite:

> **Boateng-Sarfo G. APAURA: Alternative Polyadenylation Architecture, Usage, and Regulatory Analysis. Version 0.10.1, 2026.**

A machine-readable citation is provided in [`CITATION.cff`](CITATION.cff).

---

## License

APAURA is distributed under the MIT License. See [`LICENSE`](LICENSE).
