# Algorithm overview

APAURA uses a common conceptual hierarchy:

```text
PAS -> PAU -> APA architecture -> differential APA -> candidate regulators
```

The evidence model depends on sequencing technology.

## Annotation and terminal families

APAURA loads transcript models from the GTF and organizes terminal structures in transcript orientation. Redundant/overlapping terminal annotations are reconciled into terminal families before quantification.

## Short-read evidence

Conventional fragmented RNA-seq does not provide direct cleavage evidence. APAURA therefore uses CIGAR-aware coverage/evidence across annotated terminal structures and annotated PAS groups rather than interpreting ordinary aligned read ends as cleavage sites.

## Direct/true 3′-end evidence

When transcript-oriented aligned 3′ ends are biologically meaningful, endpoint evidence can be extracted directly and nearby sites can be clustered.

## PAU

```text
PAU_gks = C_gks / sum_j C_gjs
```

```text
delta_PAU_k = mean(PAU_case,k) - mean(PAU_control,k)
```

## WPI

Sites are ordered proximal-to-distal in transcript orientation and mapped to normalized positions `x_k` from 0 to 1:

```text
WPI = sum_k PAU_k * x_k
```

```text
delta_WPI = WPI_case - WPI_control
```

```text
delta_WPI > 0 -> distal movement
delta_WPI < 0 -> proximal movement
```

## APD

```text
APD = 0.5 * sum_k |PAU_case,k - PAU_control,k|
```

APD is unsigned and captures overall redistribution. A group can have substantial APD while `delta_WPI` is near zero.

## Replicate-aware short-read statistics

APAURA separates observed effect from inferential support.

For replicated experiments:

1. count total biological replicates;
2. count replicates with any evidence;
3. count informative replicates;
4. require the configured minimum informative support;
5. use beta-binomial testing when supported;
6. do not rescue sparse replicated events with pooled Fisher testing.

Genuine 1-vs-1 experiments can use Fisher's exact test.

BH correction is applied over finite tested p-values.

Statuses include:

```text
ok
insufficient_replicate_support
insufficient_evidence
```

## IPA

```text
delta_ipa_usage = case_ipa_usage - control_ipa_usage
```

Positive values indicate relative IPA gain in the case condition; negative values indicate loss.

## RBP enrichment

APAURA extracts defined sequence windows around regulated features and scores supplied motifs against foreground/background sequence sets. These outputs are candidate regulatory associations, not direct binding evidence.
