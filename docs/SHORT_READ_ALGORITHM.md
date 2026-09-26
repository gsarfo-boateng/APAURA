# APAURA conventional short-read algorithm

## Scope

This module is for conventionally fragmented RNA-seq BAM/SAM. It **does not interpret read endpoints as cleavage sites**.

Version 0.9 directly quantifies annotated **TUTR** usage and emits an annotation architecture catalog describing potential TUTR/ALE/mixed/IPA-compatible structures.

## 1. Transcript architecture catalog

The GTF is parsed at transcript level. For each gene, APAURA records:

- number of transcript models;
- distinct terminal-exon anchors;
- distinct annotated transcript-end PAS;
- number of terminal branches with multiple PAS (TUTR-capable branches);
- transcript ends that fall within an intron of another transcript model (`IPA-compatible`).

This catalog is written to:

```text
apa_architecture.tsv
```

It is annotation structure, not a claim that every potential event is expressed or differential.

## 2. TUTR group construction

Transcripts from the same gene are grouped when their terminal-exon 5' anchors match within `--terminal-anchor-tolerance`.

This keeps PAS belonging to different terminal-exon branches separate rather than incorrectly treating ALE as one tandem 3'UTR.

Within a branch, nearby annotated transcript ends are merged within `--pas-merge-distance` and ordered proximal to distal.

## 3. Coverage segments

For `K` PAS, `K` cumulative coverage segments are created from the common terminal-exon region toward successive distal PAS.

Only segments at least `--min-segment-length` are retained.

## 4. CIGAR-aware coverage

Only aligned CIGAR blocks (`M`, `=`, `X`) contribute covered bases. `N` intronic skips and `D` deletions advance reference coordinates without contributing aligned bases.

For segment `j`:

```text
C_j = aligned_bases_j / segment_length_j
```

## 5. Monotonic cumulative profile

For tandem APA, sequence beyond a proximal PAS can only occur in transcripts using a more distal site. The expected cumulative profile is therefore:

```text
F_0 >= F_1 >= ... >= F_(K-1) >= 0
```

APAURA fits this constraint using weighted pool-adjacent-violators isotonic regression.

## 6. Site abundance reconstruction

```text
A_1 = F_0 - F_1
A_k = F_(k-1) - F_k
A_K = F_(K-1)
```

The fitted monotonic profile guarantees non-negative reconstructed abundances.

```text
PAU_k = A_k / sum_j A_j
```

## 7. Differential TUTR effect sizes

```text
deltaPAU_k = mean_case(PAU_k) - mean_control(PAU_k)
WPI = sum_k PAU_k x_k
APD = 0.5 * sum_k |mean_case(PAU_k)-mean_control(PAU_k)|
```

The short-read engine reports `shortening`, `lengthening`, `no_directional_shift`, or `stable` for each tandem group. `no_directional_shift` means PAS usage changed significantly but the net WPI did not support proximal or distal movement.

## 8. Statistical caution

The short-read engine derives continuous isoform proportions from coverage rather than observing PAS molecule counts directly. v0.9 maps these estimates onto an evidence scale for replicate-aware beta-binomial testing. The effect-size calculations are structurally well defined, but P-value calibration requires broader simulation and real-data benchmarking before v1.0.

## 9. Deliberate v0.9 boundary

Direct differential ALE and IPA inference from conventional fragmented RNA-seq is not claimed as complete in v0.9. Those structures often require unique exon evidence, splice-aware modeling, known IPA annotations or breakpoint inference. APAURA exposes the annotation architecture now but keeps the differential TUTR engine conservative rather than over-interpreting fragment endpoints.
