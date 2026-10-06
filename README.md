# APAURA

**APAURA — Alternative Polyadenylation Architecture, Usage, and Regulatory Analysis**

APAURA is a Java toolkit for detecting, quantifying, classifying, and visualizing alternative polyadenylation (APA) from multiple RNA-sequencing data types.

**Current version:** `0.10.1`

APAURA follows the biological hierarchy:

```text
PAS -> PAU -> APA architecture/change -> candidate RBP regulators
```

It uses different evidence models for conventional short-read RNA-seq, direct/true 3′-end data, and IsoQuant-corrected long reads, while reporting results through a common APA framework.

## Features

- conventional short-read RNA-seq using annotated terminal-exon coverage;
- direct RNA / true 3′-end aligned-read analysis;
- IsoQuant-corrected long-read analysis;
- tandem 3′UTR APA (TUTR);
- alternative terminal-exon architecture (ALE);
- intronic polyadenylation (IPA);
- replicate-aware differential testing;
- publication-style locus plotting;
- multi-gene locus plotting;
- persistent plotting normalization cache;
- RBP motif enrichment around regulated PAS and ALE acceptors.

## Main commands

```text
ends          Extract transcript-oriented aligned 3′ ends
run-direct    Direct RNA / true 3′-end APA analysis
run-isoquant  IsoQuant corrected-read APA analysis
short-read    Conventional RNA-seq annotated APA analysis
run-short     Alias for short-read
plot-region   Publication-style locus plots
rbp-enrich    RBP motif enrichment
citation      Print software citation
version       Print version
```

```bash
java -jar APAURA.jar --help
java -jar APAURA.jar short-read --help
java -jar APAURA.jar plot-region --help
```

## Installation

```bash
git clone <YOUR-APAURA-REPOSITORY-URL>
cd APAURA
mvn clean package
./build.sh
```

Verify the final JAR:

```bash
java -jar APAURA.jar --version
java -jar APAURA.jar short-read --help
java -jar APAURA.jar plot-region --help
jar tf APAURA.jar | grep 'org/apaura/stats/MultipleTesting.class'
```

See [docs/INSTALLATION.md](docs/INSTALLATION.md).

## Short-read quick start

Sample sheet:

```text
sample\tcondition\tbam
WT_1\tWT\t/path/WT_1.bam
WT_2\tWT\t/path/WT_2.bam
WT_3\tWT\t/path/WT_3.bam
MUT_1\tMUT\t/path/MUT_1.bam
MUT_2\tMUT\t/path/MUT_2.bam
MUT_3\tMUT\t/path/MUT_3.bam
```

```bash
java -jar APAURA.jar short-read \
  --samples samples.tsv \
  --gtf annotation.gtf \
  --control WT \
  --case MUT \
  --library reverse \
  --out apaura_results
```

For replicated short-read studies, APAURA evaluates biological replicate support before inferential testing. Sparse replicated events are not rescued by pooled Fisher testing.

See [docs/TUTORIAL_SHORT_READ.md](docs/TUTORIAL_SHORT_READ.md).

## Direct RNA / true 3′-end quick start

```bash
java -jar APAURA.jar run-direct \
  --samples samples.tsv \
  --gtf annotation.gtf \
  --control control \
  --case treated \
  --out direct_results
```

Use `--cluster-distance 0` to preserve exact single-nucleotide cleavage coordinates.

See [docs/TUTORIAL_DIRECT_RNA.md](docs/TUTORIAL_DIRECT_RNA.md).

## IsoQuant quick start

```bash
java -jar APAURA.jar run-isoquant \
  --samples isoquant_samples.tsv \
  --gtf annotation.gtf \
  --control control \
  --case treated \
  --out isoquant_results
```

See [docs/TUTORIAL_ISOQUANT.md](docs/TUTORIAL_ISOQUANT.md).

## Plot a locus

```bash
java -jar APAURA.jar plot-region \
  --results apaura_results \
  --gene CD44 \
  --mode short-read \
  --samples samples.tsv \
  --gtf annotation.gtf \
  --control WT \
  --case MUT \
  --library reverse \
  --format all
```

Multiple genes:

```bash
java -jar APAURA.jar plot-region \
  --results apaura_results \
  --gene CD44,SF3B1,TIMP4,TCL1A \
  --mode short-read \
  --samples samples.tsv \
  --gtf annotation.gtf \
  --control WT \
  --case MUT \
  --library reverse \
  --format pdf
```

Or:

```bash
--genes-file genes.txt
```

For stranded short-read data, coverage is transcript-strand aware. Genomic x-coordinates remain left-to-right. Negative-strand coverage can be mirrored vertically for display. APAURA caches normalization totals in `plot_normalization.tsv`, so later plots can avoid rescanning entire BAM files.

See [docs/TUTORIAL_PLOT_REGION.md](docs/TUTORIAL_PLOT_REGION.md).

## RBP motif enrichment

```bash
java -jar APAURA.jar rbp-enrich \
  --sites apaura_results/apa_sites.tsv \
  --genes apaura_results/apa.tsv \
  --fasta genome.fa \
  --motifs rbp_motifs.meme \
  --out rbp_results
```

Motif enrichment identifies candidate sequence-level regulatory associations; it does not prove direct binding or causality.

See [docs/TUTORIAL_RBP_ENRICHMENT.md](docs/TUTORIAL_RBP_ENRICHMENT.md).

## Core quantities

```text
PAU_gks = C_gks / sum_j C_gjs

delta_PAU = mean(PAU_case) - mean(PAU_control)

WPI = sum_k PAU_k * x_k

delta_WPI > 0 -> distal movement
delta_WPI < 0 -> proximal movement

APD = 0.5 * sum_k |PAU_case,k - PAU_control,k|
```

A gene can show substantial PAS redistribution while `delta_WPI` remains near zero. APAURA therefore separates overall change from net positional direction.

## Architecture

```text
TUTR   tandem PAS within a terminal-exon branch
ALE    alternative terminal-exon branches
mixed  both TUTR and ALE architecture
none   no terminal architecture assigned
```

IPA is tracked separately from terminal architecture.

## Replicate-aware inference

For replicated short-read experiments:

- supported events use beta-binomial testing;
- genuine 1-vs-1 experiments can use Fisher's exact test;
- sparse replicated experiments are not rescued by pooled Fisher testing;
- BH correction is applied to valid finite p-values;
- unavailable inferential values are written as `NA`.

Key statuses:

```text
ok
insufficient_replicate_support
insufficient_evidence
```

## Compact short-read outputs

```text
pas.tsv
apa_architecture.tsv
apa_sites.tsv
apa.tsv
ipa_annotation.tsv
ipa_sites.tsv
ipa.tsv
apaura_run.tsv        # when run metadata is enabled
```

Optional larger outputs:

```text
pau.tsv               # --write-per-sample
ipa_usage.tsv         # --write-per-sample
pas_coverage.tsv      # --write-diagnostics
```

See [docs/OUTPUT_SCHEMA.md](docs/OUTPUT_SCHEMA.md).

## Documentation

- [Installation](docs/INSTALLATION.md)
- [Quick start](docs/QUICKSTART.md)
- [Short-read tutorial](docs/TUTORIAL_SHORT_READ.md)
- [Direct RNA tutorial](docs/TUTORIAL_DIRECT_RNA.md)
- [IsoQuant tutorial](docs/TUTORIAL_ISOQUANT.md)
- [Plot-region tutorial](docs/TUTORIAL_PLOT_REGION.md)
- [RBP enrichment tutorial](docs/TUTORIAL_RBP_ENRICHMENT.md)
- [Batch comparisons](docs/TUTORIAL_BATCH_COMPARISONS.md)
- [Algorithm](docs/ALGORITHM.md)
- [Architecture model](docs/ARCHITECTURE_MODEL.md)
- [Output schema](docs/OUTPUT_SCHEMA.md)
- [Strandedness](docs/STRANDEDNESS.md)
- [Troubleshooting](docs/TROUBLESHOOTING.md)
- [Development](docs/DEVELOPMENT.md)
- [Validation](docs/VALIDATION.md)
- [Release checklist](docs/RELEASE_CHECKLIST.md)

## Testing

```bash
mvn clean package
./build.sh
./tests/run_tests.sh
```

## Current boundaries

For conventional fragmented short-read RNA-seq, aligned read ends are not treated as direct cleavage sites. Short-read PAS inference is annotation/coverage based.

RBP motif enrichment is associative and should not be interpreted as direct binding evidence without orthogonal validation.

## Citation

```bash
java -jar APAURA.jar --citation
```

> Boateng-Sarfo G. APAURA: Alternative Polyadenylation Architecture, Usage, and Regulatory Analysis. Version 0.10.1, 2026.

## License

MIT License. See `LICENSE`.
