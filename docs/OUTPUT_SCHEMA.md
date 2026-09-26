# APAURA output schema

APAURA separates **terminal APA architecture**, **whether APA changed**, and **the direction/mechanism of that change**. These are intentionally different fields.


## PAS region labels

The `region` column is a site-level genomic/transcript annotation and is separate from `apa_architecture`. APAURA writes the following human-readable values:

- `3UTR` — PAS lies in an annotated or inferred 3'UTR within a terminal exon.
- `terminal_exon_noncoding` — PAS lies in the terminal exon of a noncoding transcript/gene.
- `coding_exon` — PAS lies in a coding exon outside the annotated 3'UTR.
- `noncoding_exon` — PAS lies in a noncoding exon that is not classified as a terminal noncoding exon.
- `intron` — PAS lies in an intron and can contribute to IPA analysis.
- `downstream` — PAS lies just beyond the annotated transcript 3' end within the allowed annotation window.
- `intergenic` — no compatible gene annotation was assigned.

These labels answer **where is this PAS?** By contrast, `apa_architecture` answers **how are the gene's terminal PAS/exons organized?** The same `3UTR` site can participate in TUTR, ALE, or mixed architecture depending on the other terminal sites and terminal-exon branches in that gene.

## Terminal APA architecture

`apa_architecture` describes the organization of poly(A) sites across terminal-exon branches. The allowed values are:

- `TUTR` — tandem 3'UTR APA: at least one terminal-exon branch contains multiple PAS.
- `ALE` — alternative last exon: the gene contains multiple terminal-exon branches.
- `mixed` — both TUTR and ALE architecture are present.
- `none` — no TUTR/ALE terminal architecture was identified from the retained sites.

`architecture_components` contains only terminal architecture components (`TUTR`, `ALE`, or `TUTR+ALE`). IPA is not an `apa_architecture` value.

## APA change

`apa_change` is `yes` or `no` and indicates whether the gene/group passes the configured statistical and effect-size criteria for at least one APA mechanism.

This field is useful because a PAS-usage pattern can change significantly without having a net proximal or distal shift.

## APA direction

`apa_direction` is a human-readable description of the supported direction:

- `shortening`
- `lengthening`
- `proximal_switch`
- `distal_switch`
- `switch`
- `gain`
- `loss`
- `complex`
- `no_directional_shift`
- `stable`
- `not_applicable`

`no_directional_shift` is not a separate APA architecture. It means that APA usage changed, but the net positional statistic did not support either shortening or lengthening.

## Mechanism-specific direction columns

- `tutr_direction`: `shortening`, `lengthening`, `no_directional_shift`, `stable`, or `not_applicable`.
- `ale_direction`: `proximal_switch`, `distal_switch`, `switch`, `stable`, or `not_applicable`.
- `ipa_direction`: `gain`, `loss`, `stable`, or `not_applicable`.
- `internal_exon_direction`: `gain`, `loss`, `stable`, or `not_applicable`.

IPA is therefore represented independently from terminal APA architecture. For example, an IPA-only gene may have:

```text
apa_architecture = none
apa_change       = yes
apa_direction    = gain
ipa_direction    = gain
apa_event        = IPA:gain
```

## Mixed/complex changes

A gene with both TUTR and ALE architecture can have:

```text
apa_architecture    = mixed
architecture_components = TUTR+ALE
apa_change          = yes
apa_direction       = complex
tutr_direction      = shortening
ale_direction       = distal_switch
changed_mechanisms  = ALE+TUTR
mechanism_directions = ALE=distal_switch;TUTR=shortening
apa_event           = mixed:complex
```

`changed_mechanisms` records which mechanisms passed the change criteria. `mechanism_directions` preserves their individual directions.

## Convenience event field

`apa_event` is a convenience label, not the primary biological representation. Examples include:

```text
TUTR:shortening
TUTR:lengthening
ALE:distal_switch
IPA:gain
mixed:complex
multi_mechanism:complex
TUTR:no_directional_shift
```

`multi_mechanism:complex` is used when multiple mechanisms change but the terminal architecture itself is not `mixed` (for example TUTR shortening plus IPA gain). Users should generally filter/model using `apa_architecture`, `apa_change`, `changed_mechanisms`, and the mechanism-specific direction fields.

## TUTR direction and usage-change magnitude

For a TUTR-capable terminal branch, APAURA reports two complementary quantities.

The signed positional statistic is:

```text
delta_tutr_wpi = tutr_wpi_case - tutr_wpi_control
```

Positive values indicate distal movement (lengthening), while negative values indicate proximal movement (shortening).

`tutr_shift_magnitude` is the branch-usage-weighted total-variation distance between control and case PAS-usage distributions within TUTR-capable branches:

```text
TVD = 0.5 * sum_k |PAU_case,k - PAU_control,k|
```

A gene can therefore have a large `tutr_shift_magnitude` while `delta_tutr_wpi` is approximately zero. In that case APAURA reports:

```text
apa_change     = yes
tutr_direction = no_directional_shift
```

rather than inventing a separate biological APA class.

## Termination-compartment metrics

The terminal-exon/IPA/internal-exon/downstream usage vector is summarized independently of `apa_architecture`.

- `termination_shift_index` — total-variation distance between condition-level termination-compartment usage vectors.
- `termination_jsd` — Jensen-Shannon divergence between those vectors.
- `dominant_control_termination`
- `dominant_case_termination`
- `termination_transition`

These fields quantify termination-site movement and are especially useful for IPA.

## Short-read output

The conventional RNA-seq `apa.tsv` uses the same naming convention:

```text
scope = TANDEM_3UTR
apa_architecture = TUTR
apa_change = yes
apa_direction = shortening
apa_event = TUTR:shortening
```

The short-read engine currently quantifies annotated tandem terminal-exon APA groups; ALE/IPA short-read inference remains a separate development target.


## Canonical output filenames

APAURA v0.10.1 uses biological names for its primary result files.

| File | Meaning |
|---|---|
| `pas.bed` | retained single-base PAS calls for direct/long-read input |
| `pas_all.bed` | all pooled PAS calls before PAS filtering |
| `pas.tsv` | annotated PAS catalog/groups for conventional short-read analysis |
| `pas_counts.tsv` | per-sample PAS read counts |
| `pas_coverage.tsv` | coverage evidence used for short-read PAS deconvolution |
| `pau.tsv` | per-sample polyadenylation-site usage |
| `apa_sites.tsv` | site-level differential APA results |
| `apa.tsv` | main gene/group-level APA results, including architecture and direction |
| `apa_architecture.tsv` | focused architecture/annotation summary |
| `rbp.tsv` | RBP motif-enrichment results |
| `rbp_hits.tsv` | site-level RBP motif hits |
| `rbp_profile.tsv` | positional RBP motif profile |
| `rbp_windows.tsv` | sequences scanned by the RBP module |

The intention is that a user can identify the analysis layer directly from the filename: PAS → PAU → APA → RBP.
