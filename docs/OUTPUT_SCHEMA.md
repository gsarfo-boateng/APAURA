# Output schema

This document summarizes the primary APAURA v0.10.1 output tables. `NA`
represents an unavailable or non-evaluable quantity; it must not be interpreted
as a non-significant value.

## Conventional short-read RNA-seq

The primary compact files are:

```text
apa.tsv
apa_events.tsv
apa_sites.tsv
```

Supporting tables are written under `details/` when requested.

### `apa.tsv`

One row per gene:

```text
gene_id
gene_name
biotype
chrom
gene_start_1based
gene_end_1based
strand
apa_architecture
apa_change
apa_direction
changed_mechanisms
tutr_change
tutr_direction
delta_tutr_wpi
tutr_shift_magnitude
tutr_jsd
tutr_min_fdr
ale_change
ale_direction
delta_ale_wpi
ale_switch_index
ale_min_fdr
ipa_change
ipa_direction
delta_ipa_usage
ipa_min_fdr
min_mechanism_fdr
status
```

`apa_architecture` is one of `SINGLE_PAS`, `TUTR`, `ALE`, `mixed`, or
`UNRESOLVED` where that state is applicable to the analyzed structure.
`changed_mechanisms` is mechanism-specific and may contain combinations such as
`IPA+TUTR`.

### `apa_events.tsv`

One row per tested TUTR/ALE/IPA event:

```text
event_id
gene_id
gene_name
biotype
chrom
strand
event_type
start_1based
end_1based
control_usage
case_usage
delta_usage
effect_size
effect_metric
pvalue
fdr
direction
status
```

### `apa_sites.tsv`

Terminal PAS-level differential results:

```text
group_id
gene_id
gene_name
chrom
strand
site_id
pas_position_1based
site_rank
control_pau
case_pau
delta_pau
pvalue
fdr
status
```

Short-read fragment ends are not interpreted as cleavage sites. These tables are
derived from annotation-aware, CIGAR-aware coverage/evidence and the short-read
deconvolution model.

### Short-read status values

Common values include:

```text
ok
insufficient_evidence
insufficient_replicate_support
not_tested
```

Replicated events are not rescued with pooled Fisher testing when biological
replicate support is inadequate.

## Direct RNA and long-read / IsoQuant

Typical outputs are:

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

### `apa_sites.tsv`

```text
site_id
gene_id
gene_name
biotype
chrom
pos0
strand
region
architecture_role
terminal_exon
terminal_exon_anchor0
terminal_branch_rank
site_rank
control_pau
case_pau
delta_pau
control_n_evaluable
case_n_evaluable
pvalue
fdr
dispersion_rho
method
status
```

`pos0` is zero-based. `terminal_exon_anchor0` preserves the source terminal-exon
anchor; structural classification uses the reconciled terminal structure
internally. `terminal_branch_rank` is the v0.10.1 compatibility name for the
reconciled terminal exon group rank.

Site-level `method` values include:

```text
beta-binomial-LRT
Fisher-exact-unreplicated
not_tested
```

`Fisher-exact-unreplicated` is used only when the **experimental design itself**
is 1-vs-1 and both samples are evaluable. A replicated design that collapses to
one evaluable sample per condition is `not_tested` with
`insufficient_replicate_support`.

### `apa.tsv`

One row per gene:

```text
gene_id
gene_name
biotype
apa_architecture
architecture_components
apa_change
apa_direction
tutr_direction
ale_direction
ipa_direction
internal_exon_direction
apa_event
changed_mechanisms
mechanism_directions
n_sites
n_terminal_branches
n_tutr_branches
control_wpi
case_wpi
delta_wpi
apd
usage_jsd
control_terminal_usage
case_terminal_usage
delta_terminal_usage
control_ipa_usage
case_ipa_usage
delta_ipa_usage
termination_shift_index
termination_jsd
ale_switch_index
ale_jsd
control_ale_wpi
case_ale_wpi
delta_ale_wpi
control_tutr_wpi
case_tutr_wpi
delta_tutr_wpi
tutr_shift_magnitude
dominant_control_termination
dominant_case_termination
termination_transition
ks_d
ks_signed_d
ks_max_position_1based
ks_z
ks_pvalue
ks_fdr
min_site_fdr
```

The main usage quantities are:

```text
PAU       poly(A)-site usage
delta_PAU case mean PAU - control mean PAU
WPI       weighted position index
delta_WPI case WPI - control WPI
APD       total-variation distance across the PAS distribution
usage_jsd Jensen-Shannon divergence across the PAS distribution
```

Architecture-specific fields summarize TUTR, ALE, IPA, and
termination-compartment changes. `changed_mechanisms` requires significant
support from the corresponding mechanism; `min_site_fdr` from an unrelated
mechanism does not certify another mechanism.

### `pau.tsv`

Per-sample site usage:

```text
sample
condition
site_id
gene_id
gene_name
chrom
pos0
strand
region
site_rank
count
gene_end_count
pau
```

### `pas_counts.tsv`

The fixed leading columns are:

```text
site_id
chrom
pos0
strand
gene_id
gene_name
region
```

followed by one count column per sample.

## Single-cell

The integrated `sc` workflow preserves sc3UTR outputs and adds harmonized APAURA
usage tables:

```text
apaura_usage.tsv
apaura_site_usage.tsv
```

`apaura_site_usage.tsv` contains site-level PAU summaries. `apaura_usage.tsv`
contains gene/group-level PAU-distribution summaries including APAURA
distributional metrics such as TVD/JSD and weighted-position behavior.

Biological samples, not individual cells, are the replicate unit for formal
condition-level inference. Unreplicated single-cell comparisons remain
descriptive rather than being presented as replicate-supported population
inference.

Because single-cell grouping columns depend on the selected `--group-by`
hierarchy, the header written by each run is the canonical schema for those
group-specific fields.

## RBP enrichment

RBP analysis writes:

```text
rbp.tsv
rbp_hits.tsv
rbp_profile.tsv
rbp_windows.tsv
rbp_run_info.txt
```

See `docs/RBP_ENRICHMENT.md`.
