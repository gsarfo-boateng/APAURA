# Quick start

## Short-read RNA-seq

```bash
java -jar APAURA.jar short-read \
  --samples samples.tsv \
  --gtf annotation.gtf \
  --control control \
  --case treated \
  --library reverse \
  --out results
```

## Direct RNA / true 3′-end

```bash
java -jar APAURA.jar run-direct \
  --samples samples.tsv \
  --gtf annotation.gtf \
  --control control \
  --case treated \
  --out direct_results
```

## IsoQuant

```bash
java -jar APAURA.jar run-isoquant \
  --samples isoquant_samples.tsv \
  --gtf annotation.gtf \
  --control control \
  --case treated \
  --out isoquant_results
```

## Plot one gene

```bash
java -jar APAURA.jar plot-region \
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
java -jar APAURA.jar plot-region \
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
java -jar APAURA.jar rbp-enrich \
  --sites results/apa_sites.tsv \
  --genes results/apa.tsv \
  --fasta genome.fa \
  --motifs rbp_motifs.meme \
  --out results/rbp
```
