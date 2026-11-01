# APAURA `plot-region` tutorial

`plot-region` generates locus-level APA figures from either a completed APAURA analysis or directly from BAM/SAM alignments plus a GTF.

## Requirements

The plotting workflow requires `Rscript`, `samtools`, and the R/Bioconductor packages documented in the main README.

Check the complete CLI:

```bash
apaura plot-region --help
```

## Completed-results mode

```bash
apaura plot-region --results results --gene ABCDEF --samples samples.tsv --gtf annotation.gtf --control U1_WT --case U1_MUT --library reverse --show-replicates --show-overview --format pdf --out plots/HNRNPA2B1
```

Depending on the analysis and locus, APAURA can display mean normalized coverage, replicate coverage, whole-gene overview, terminal-region zoom, tested PAS/ΔPAU, TUTR direction, IPA usage, ALE branch usage, and terminal architecture.

## Standalone mode

`--results` is optional.

```bash
apaura plot-region --samples samples.tsv --gtf annotation.gtf --gene ABDCEF --control U1_WT --case U1_MUT --library reverse --show-replicates --show-overview --format pdf --out plots/HNRNPA2B1
```

Standalone mode displays coverage and annotation-derived terminal structure. It **does not run differential APA analysis** and therefore intentionally does not claim ΔPAU, significant TUTR shortening/lengthening, ALE switching, IPA gain/loss, or formal TUTR/ALE/mixed architecture when completed APAURA results are absent.

## Multiple genes

```bash
apaura plot-region --samples samples.tsv --gtf annotation.gtf --gene HNRNPA2B1,PLA2G12A,TMEM254 --control U1_WT --case U1_MUT --library reverse --show-overview --format pdf --out plots
```

Or use `--genes-file genes.txt`.

## Genomic region

```bash
apaura plot-region --samples samples.tsv --gtf annotation.gtf --region chr7:26170000-26197000 --control U1_WT --case U1_MUT --library reverse --format pdf --out plots/chr7_region
```

Standalone region plotting currently uses the supplied GTF for annotation/terminal-structure context.

## Strand orientation

The genomic x-axis remains left-to-right. For a minus-strand gene, transcript 3' direction points left. Use `--genomic-orientation` when genomic rather than transcript-oriented coverage display is preferred.

## Coverage options

- `--show-replicates` adds individual sample tracks.
- `--show-overview` adds a whole-gene overview above the terminal-region zoom.
- `--verbose` prints per-sample/cache messages.

APAURA caches normalization totals so repeated plotting does not require recounting the same alignment totals.

## Output formats

```text
--format pdf
--format png
--format svg
--format all
```

PDF or SVG is recommended for publication figures.

## Interpretation

A standalone plot answers: **What does the coverage and annotated terminal structure look like at this locus?**

A completed-results plot can additionally answer: **Which PAS, terminal branches, or intronic termination events changed, by how much, and in which direction?**

Do not infer differential APA significance from the standalone coverage plot alone.

## Citation

> Boateng-Sarfo G, E Niewold, Sigor S, Wang L. APAURA: Alternative Polyadenylation Architecture, Usage, and Regulatory Analysis. Version 0.10.1, 2026.
