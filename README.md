# APAURA

**APAURA - Alternative Polyadenylation Architecture, Usage, and Regulatory Analysis** is a Java toolkit for detecting and characterizing alternative polyadenylation (APA) changes from RNA-sequencing data.

Version: **0.10.1**

APAURA follows the biological workflow:

```text
PAS -> PAU -> APA architecture/change -> candidate RBP regulators
```

It uses different evidence models for long-read/true 3'-end data and conventional short-read RNA-seq, then reports results through a common APA framework.

## APA model

APAURA separates three concepts:

1. **PAS region** - where an individual poly(A) site occurs.
2. **APA architecture** - how terminal exons and terminal PAS are organized within a gene.
3. **APA direction** - how usage changes between conditions.

### PAS regions

Site-level `region` values are:

```text
3UTR
terminal_exon_noncoding
coding_exon
noncoding_exon
intron
downstream
intergenic
```

### Terminal APA architecture

Gene-level `apa_architecture` values are:

```text
TUTR multiple PAS within the same terminal-exon branch
ALE alternative terminal-exon branches
mixed both TUTR and ALE architecture
none  neither TUTR nor ALE is identified
```

Intronic polyadenylation (IPA) is tracked separately rather than treated as a terminal-exon architecture.

### Direction of change

Examples include:

```text
TUTR -> shortening / lengthening
ALE -> proximal_switch / distal_switch
IPA -> gain / loss
mixed or multi-mechanism change -> complex
```

Mechanism-specific columns (`tutr_direction`, `ale_direction`, `ipa_direction`) preserve simultaneous changes within the same gene.

## Supported analysis modes

### Direct RNA / true 3'-end evidence

Use `run-direct` for BAM/SAM data in which the transcript-oriented aligned 3' end is meaningful, including Nanopore direct RNA and suitable 3'-end-enriched assays.

```bash
java -jar APAURA.jar run-direct --samples samples.tsv --gtf annotation.gtf --control control --case treated --out results
```

To preserve exact single-nucleotide cleavage coordinates:

```bash
--cluster-distance 0
```

### IsoQuant corrected reads

```bash
java -jar APAURA.jar run-isoquant --samples isoquant_samples.tsv --gtf annotation.gtf --control control --case treated --out results
```

### Conventional short-read RNA-seq

```bash
java -jar APAURA.jar short-read --samples samples.tsv --gtf annotation.gtf --control control --case treated --library unstranded --out results
```

Short-read endpoints are not treated as cleavage sites. APAURA uses CIGAR-aware terminal-exon coverage and annotated PAS structure to estimate tandem-PAS usage.

### RBP motif enrichment

```bash
java -jar APAURA.jar rbp-enrich --sites results/apa_sites.tsv --genes results/apa.tsv --fasta genome.fa --motifs rbp_motifs.meme --out rbp_results
```

The RBP module analyzes motif enrichment around regulated PAS and, for ALE events, terminal-exon splice acceptors. Results are candidate regulatory associations, not proof of direct binding or causality.

## Sample sheets

Direct/short-read example:

```text
sample	condition	bam
CTL_1	control	/path/CTL_1.bam
CTL_2	control	/path/CTL_2.bam
TRT_1	treated	/path/TRT_1.bam
TRT_2	treated	/path/TRT_2.bam
```

Accepted path columns are `bam`, `file`, and `path`. Relative paths are resolved relative to the sample sheet.

IsoQuant example:

```text
sample	condition	isoquant_dir
CTL_1	control	/path/CTL_1_isoquant
CTL_2	control	/path/CTL_2_isoquant
TRT_1	treated	/path/TRT_1_isoquant
TRT_2	treated	/path/TRT_2_isoquant
```

## Core outputs

Long-read/direct analyses produce:

```text
pas_all.bed PAS before support filtering
pas.bed retained representative PAS
pas_counts.tsv PAS counts by sample
pau.tsv per-sample poly(A)-site usage
apa_sites.tsv PAS-level differential results
apa.tsv main gene-level APA results
apa_architecture.tsv
run_info.txt
```

The short-read engine produces analogous PAS/PAU/APA tables based on annotated terminal-exon coverage.

RBP analysis produces:

```text
rbp.tsv
rbp_hits.tsv
rbp_profile.tsv
rbp_windows.tsv
rbp_run_info.txt
```

See `docs/OUTPUT_SCHEMA.md` for column definitions.

## Core quantities

For poly(A) site `k` in gene `g`, sample `s`:

```text
PAU_gks = C_gks / sum_j C_gjs
```

Site-level change:

```text
delta_PAU = mean(PAU_case) - mean(PAU_control)
```

Sites are ordered proximal-to-distal in transcript orientation and assigned normalized positions `x_k` from 0 to 1:

```text
WPI = sum_k PAU_k * x_k
```

Therefore:

```text
delta_WPI > 0 -> distal movement
delta_WPI < 0 -> proximal movement
```

The overall PAS-usage difference is summarized by total-variation distance:

```text
APD = 0.5 * sum_k |PAU_case,k - PAU_control,k|
```

APAURA also reports architecture-specific metrics for ALE switching, TUTR shifts, IPA usage, and termination-compartment changes. See `docs/ARCHITECTURE_MODEL.md` and `docs/ALGORITHM.md`.

## Statistical methods

Depending on the analysis and available replication, APAURA uses:

- beta-binomial likelihood-ratio testing for PAS usage;
- Fisher's exact test when replicate-aware fitting is not available;
- Benjamini-Hochberg FDR correction;
- a tie-aware two-sample ECDF maximum-distance statistic at gene level;
- total-variation distance and Jensen-Shannon divergence for distributional summaries;
- weighted isotonic regression for short-read coverage deconvolution.

## Help

Global help:

```bash
java -jar APAURA.jar --help
```

Command-specific help:

```bash
java -jar APAURA.jar run-direct --help
java -jar APAURA.jar run-isoquant --help
java -jar APAURA.jar short-read --help
java -jar APAURA.jar rbp-enrich --help
java -jar APAURA.jar ends --help
```

## Building from source

APAURA targets Java 8-compatible bytecode. JDK 17 or 21 can be used for development.

With Maven:

```bash
mvn clean package
```

This runs the JUnit test suite and creates:

```text
target/APAURA.jar
```

The portable helper can also be used:

```bash
./build.sh
```

If Maven is available, it runs the Maven build. Otherwise it compiles the main source directly with `javac` and `jar`.

## Testing

### Unit tests

JUnit tests for the mathematical and low-level core live in:

```text
src/test/java/org/apaura/
```

Run them with:

```bash
mvn test
```

### End-to-end synthetic regression tests

```bash
./tests/run_tests.sh
```

The integration suite covers direct-RNA/SAM, BAM parsing, IsoQuant, short-read coverage analysis, TUTR, ALE, IPA, mixed mechanisms, and RBP motif enrichment.

## Alignment readers

APAURA contains small sequential SAM/BAM/BGZF readers so the core JAR can remain dependency-light. The BAM reader extracts only fields required by APAURA and does not implement random access or the complete BAM optional-tag API. See `docs/DEVELOPMENT.md` for scope and maintenance notes.

## Plotting

The companion R script is:

```text
scripts/plot_apaura.R
```

Example:

```bash
Rscript scripts/plot_apaura.R --input results --rbp rbp_results --out results/figures
```

See `docs/PLOTTING.md`.

## Current boundaries

Version 0.10.1 does not yet treat conventional fragmented RNA-seq as direct cleavage-site evidence, perform de novo short-read PAS breakpoint discovery, or claim complete short-read ALE/IPA quantification. Internal-priming filtering is also not yet implemented.

## Documentation

```text
docs/ALGORITHM.md
docs/ARCHITECTURE_MODEL.md
docs/LONG_READ_ALGORITHM.md
docs/SHORT_READ_ALGORITHM.md
docs/RBP_ENRICHMENT.md
docs/OUTPUT_SCHEMA.md
docs/PLOTTING.md
docs/VALIDATION.md
docs/DEVELOPMENT.md
```

Additional method-positioning and implementation notes are available under `docs/` for developers and reviewers.

## Citation

```bash
java -jar APAURA.jar --citation
```

Current software citation:

> Boateng-Sarfo G. APAURA: Alternative Polyadenylation Architecture, Usage, and Regulatory Analysis. Version 0.10.1, 2026.

## License

MIT License. See `LICENSE`.
