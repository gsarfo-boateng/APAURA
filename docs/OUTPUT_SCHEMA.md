# Output schema

This document describes the compact short-read result tables used by APAURA v0.10.1.

Unavailable inferential values are written as `NA`.

## `apa_sites.tsv`

```text
group_id
gene_id
gene_name
biotype
scope
site_id
chrom
pas_position_1based
strand
site_rank
normalized_position
control_pau
case_pau
delta_pau
control_pau_sd
case_pau_sd
control_n_total
case_n_total
control_n_any_evidence
case_n_any_evidence
control_n_informative
case_n_informative
control_evidence
case_evidence
pvalue
fdr
method
status
```

`pas_position_1based` is the canonical user-facing coordinate.

## `apa.tsv`

```text
group_id
gene_id
gene_name
biotype
chrom
strand
scope
apa_architecture
n_sites
control_wpi
case_wpi
delta_wpi
control_wpi_sd
case_wpi_sd
apd
control_n_total
case_n_total
control_n_any_evidence
case_n_any_evidence
control_n_informative
case_n_informative
min_fdr
apa_change
observed_direction
apa_direction
apa_event
status
```

## `ipa_sites.tsv`

```text
ipa_id
gene_id
gene_name
biotype
chrom
strand
intron_start0
intron_end0
representative_pas_1based
n_annotated_ipa_pas
annotated_ipa_positions_1based
control_ipa_usage
case_ipa_usage
delta_ipa_usage
control_n_total
case_n_total
control_n_any_evidence
case_n_any_evidence
control_n_informative
case_n_informative
control_evidence
case_evidence
pvalue
fdr
dispersion_rho
ipa_direction
method
status
```

## Default compact short-read files

```text
pas.tsv
apa_architecture.tsv
apa_sites.tsv
apa.tsv
ipa_annotation.tsv
ipa_sites.tsv
ipa.tsv
```

## Optional large files

`--write-per-sample`:

```text
pau.tsv
ipa_usage.tsv
```

`--write-diagnostics`:

```text
pas_coverage.tsv
```

## Direct/long-read files

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
