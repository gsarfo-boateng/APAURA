# Quick start

## Short-read RNA-seq

```bash
apaura short-read \
  --samples samples.tsv \
  --gtf annotation.gtf \
  --control control \
  --case treated \
  --library reverse \
  --out results
```

## Direct RNA / true 3'-end

```bash
apaura run-direct \
  --samples samples.tsv \
  --gtf annotation.gtf \
  --control control \
  --case treated \
  --out direct_results
```

## IsoQuant

```bash
apaura run-isoquant \
  --samples isoquant_samples.tsv \
  --gtf annotation.gtf \
  --control control \
  --case treated \
  --out isoquant_results
```

## Plot one gene

```bash
apaura plot-region \
  --results results \
  --gene CD44 \
  --mode short-read \
  --samples samples.tsv \
  --gtf annotation.gtf \
  --control control \
  --case treated \
  --library reverse
```

## Plot many genes

```bash
apaura plot-region \
  --results results \
  --gene CD44,SF3B1,TIMP4,TCL1A \
  --mode short-read \
  --samples samples.tsv \
  --gtf annotation.gtf \
  --control control \
  --case treated \
  --library reverse \
  --format pdf
```

## RBP enrichment

```bash
apaura rbp-enrich \
  --sites results/apa_sites.tsv \
  --genes results/apa.tsv \
  --fasta genome.fa \
  --motifs rbp_motifs.meme \
  --out results/rbp
```

## Standalone locus visualization

A completed APAURA result directory is not required for coverage/annotation plotting:

```bash
apaura plot-region \
  --samples samples.tsv \
  --gtf annotation.gtf \
  --gene HNRNPA2B1 \
  --control control \
  --case treated \
  --library reverse \
  --show-overview \
  --format pdf \
  --out plots/HNRNPA2B1
```

Standalone mode visualizes BAM/SAM coverage and GTF-derived terminal structure. It does not perform differential APA testing and therefore does not make ΔPAU, TUTR, ALE, IPA, or formal APA-architecture calls without completed APAURA results.

See `docs/TUTORIAL_PLOT_REGION.md`.
