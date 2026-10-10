# Tutorial: direct RNA / true 3'-end evidence

Use `run-direct` when the transcript-oriented aligned 3' end is biologically meaningful, including Nanopore direct RNA and appropriate 3'-end-enriched assays.

## Sample sheet

```text
sample\tcondition\tbam
C1\tcontrol\t/path/C1.bam
C2\tcontrol\t/path/C2.bam
T1\ttreated\t/path/T1.bam
T2\ttreated\t/path/T2.bam
```

## Run

```bash
java -jar APAURA.jar run-direct \
  --samples samples.tsv \
  --gtf annotation.gtf \
  --control control \
  --case treated \
  --out direct_results
```

Use:

```bash
--cluster-distance 0
```

to preserve exact single-nucleotide endpoint coordinates.

Common controls include `--cluster-distance`, `--cluster-representative`, `--min-gene-reads`, `--min-site-reads`, `--min-pau`, `--min-mapq`, `--keep-duplicates`, `--strand-mode`, `--annotation-slop`, `--min-delta-pau`, `--fdr`, `--wpi-threshold`, `--min-internal-exon-delta`, `--min-ale-switch`, `--min-ipa-delta`, and `--min-tutr-shift`.

Check exact help with:

```bash
java -jar APAURA.jar run-direct --help
```

Typical outputs include:

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
