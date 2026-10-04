# APAURA long-read / direct-RNA algorithm

## 1. Alignment filtering

BAM/SAM processing removes unmapped, secondary, supplementary and QC-failed records. Duplicate-flagged records are excluded by default. `--min-mapq` is applied before endpoint extraction.

## 2. Transcript-oriented 3' endpoint

For a plus-strand alignment:

```text
p0 = alignmentEndExclusive - 1
```

For a minus-strand alignment:

```text
p0 = alignmentStart0
```

Each observation is therefore a single genomic nucleotide in 0-based coordinates.

## 3. PAS clustering

Endpoints are pooled by chromosome and strand. `--cluster-distance 0` retains exact coordinates. Positive distances cluster nearby cleavage microheterogeneity.

Representative coordinate choices:

```text
mode    pooled coordinate with maximum support
median  weighted median coordinate
```

Counts from all member coordinates are retained at the representative site.

## 4. Structural annotation

Sites are assigned to same-strand genes and categorized as:

```text
3UTR
terminal_exon_noncoding
coding_exon
noncoding_exon
intron
downstream
```

Terminal-exon classification is transcript-aware. In particular, coding context is assessed within the matched terminal exon rather than using a single gene-global CDS endpoint, which is necessary for alternative-last-exon structures.

## 5. PAS support filters

Defaults:

```text
min gene depth per condition = 10
min pooled site reads = 5
min pooled PAU = 0.01
```

The unfiltered PAS set is retained in `pas_all.bed`.

## 6. PAU

For gene `g`, site `k`, sample `s`:

```text
PAU_gks = C_gks / N_gs
```

where `N_gs` is the annotated gene-end depth before low-PAS filtering.

Site effect size:

```text
deltaPAU = mean(PAU_case) - mean(PAU_control)
```

## 7. Site-level inference

With at least two informative replicates per condition, APAURA fits beta-binomial null and alternative models and uses a likelihood-ratio statistic.

Without sufficient replication, a two-sided Fisher exact test is used on PAS-versus-other-gene-end counts.

P values are corrected with Benjamini-Hochberg.

## 8. Hierarchical APA architecture

After site filtering, each gene is decomposed into:

```text
Level 1: termination compartment
Level 2: terminal-exon branch (ALE)
Level 3: PAS within terminal branch (TUTR)
```

See `ARCHITECTURE_MODEL.md` for full equations.

## 9. Global PAS summaries

### WPI

Retained PAS are ordered proximal to distal:

```text
x_k = (k-1)/(K-1)
WPI = sum_k PAU*_k x_k
```

where `PAU*` is renormalized over retained PAS.

### APD

```text
APD = 0.5 * sum_k |mean_case(PAU*_k) - mean_control(PAU*_k)|
```

WPI and APD are retained as useful whole-gene summaries but are not treated as sufficient mechanistic descriptions for mixed architecture genes.

## 10. Gene-level distribution diagnostic

APAURA also compares pooled control and case PAS distributions using a discrete empirical-CDF maximum-distance statistic. The sign is oriented in transcript space, including reversal on minus-strand genes. A tie-aware variance term is used for the normal approximation and resulting gene P values are BH corrected.

This statistic is a diagnostic/global test; architecture-specific effect sizes remain separately reported.

## 11. Architecture and direction selection

APAURA stores architecture and direction separately. Examples:

```text
apa_architecture=TUTR   apa_direction=shortening
apa_architecture=ALE    apa_direction=distal_switch
apa_architecture=none   ipa_direction=gain   apa_direction=gain
```

For genes with more than one active mechanism, the overall direction is `complex`, and mechanism-specific columns retain the individual changes:

```text
apa_architecture=mixed
apa_direction=complex
ale_direction=distal_switch
tutr_direction=shortening
mechanism_directions=ALE=distal_switch;TUTR=shortening
apa_event=mixed:complex
```

This avoids conflating biological architecture with change direction.
