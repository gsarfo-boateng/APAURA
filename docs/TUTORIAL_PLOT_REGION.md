# Tutorial: `plot-region`

`plot-region` aligns RNA evidence, terminal PAS changes, IPA changes, and gene/transcript structure on one genomic coordinate system.

## Single gene

```bash
java -jar APAURA.jar plot-region \
  --results results \
  --gene CD44 \
  --mode short-read \
  --samples samples.tsv \
  --gtf annotation.gtf \
  --control WT \
  --case MUT \
  --library reverse \
  --format all
```

## Multiple genes

```bash
java -jar APAURA.jar plot-region \
  --results results \
  --gene CD44,SF3B1,TIMP4,TCL1A \
  --mode short-read \
  --samples samples.tsv \
  --gtf annotation.gtf \
  --control WT \
  --case MUT \
  --library reverse \
  --format pdf
```

Or use a text file:

```bash
--genes-file genes.txt
```

## Genomic orientation

Genomic x-coordinates remain left-to-right.

For negative-strand genes in short-read mode, coverage can be mirrored below zero for visual clarity. Mirroring changes only the display sign; it does not alter genomic coordinates and does not mean the biological coverage is negative.

Use:

```bash
--genomic-orientation
```

to keep negative-strand coverage upright.

## Strand-aware coverage

For `--library forward` or `--library reverse`, paired-end read orientations consistent with the target transcript strand are selected before locus coverage is calculated.

This prevents antisense/overlapping transcription from being mixed into the target track.

## Persistent normalization cache

Exact assigned-read normalization is cached in:

```text
plot_normalization.tsv
```

Typical behavior:

```text
first locus -> cache miss -> each BAM counted once
later loci -> cache hit -> only regional extraction/depth
```

The cache is tied to sample/BAM and relevant plotting filters so incompatible normalization values are not silently reused.

## Replicates

Use:

```bash
--show-replicates
```

to overlay faint individual replicate traces. Replicates are faceted with their own condition so control/case colors do not mix.

## Transcript models

```bash
--show-transcripts
--max-transcripts 20
```

## Reading the panels

Top to bottom:

1. normalized RNA coverage/evidence;
2. terminal PAS `delta_PAU`;
3. IPA `delta_ipa_usage`;
4. gene/transcript model.

```text
delta_PAU > 0 -> higher relative PAS use in case
delta_PAU < 0 -> lower relative PAS use in case

delta_ipa_usage > 0 -> IPA gain in case
delta_ipa_usage < 0 -> IPA loss in case
```

For minus-strand genes, proximal/distal interpretation must follow transcript orientation rather than page left/right.
