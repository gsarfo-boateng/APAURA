# APAURA APA architecture model

## Overview

APAURA separates two biological questions that are often conflated:

1. **What terminal APA architecture does the gene have?**
2. **What polyadenylation/termination mechanism changed between conditions?**

The terminal architecture vocabulary is intentionally small:

- `TUTR` — multiple PAS within the same terminal-exon branch;
- `ALE` — multiple alternative terminal-exon branches;
- `mixed` — both TUTR and ALE are present;
- `none` — neither TUTR nor ALE is identified.

IPA is not an architecture value. Intronic polyadenylation is quantified as a separate termination mechanism.

## Hierarchical representation

```text
gene
 |
 +-- terminal architecture
 |    +-- terminal branch 1
 |    |    +-- PAS 1
 |    |    +-- PAS 2
 |    +-- terminal branch 2
 |         +-- PAS 3
 |
 +-- non-terminal termination
      +-- intronic PAS (IPA)
      +-- internal-exon PAS
      +-- downstream end
```

This allows APAURA to describe a gene with both an ALE switch and TUTR shortening without collapsing the two effects into a single proximal/distal score.

## TUTR

Within each terminal branch containing at least two PAS, PAS usage is normalized:

```text
PAU_k = count_k / sum_j count_j
```

Sites are ordered from proximal to distal in transcript orientation and assigned normalized positions `x_k` between 0 and 1.

```text
TUTR_WPI = sum_k PAU_k * x_k
```

The condition difference is:

```text
delta_tutr_wpi = TUTR_WPI_case - TUTR_WPI_control
```

- `delta_tutr_wpi > 0` supports `lengthening`.
- `delta_tutr_wpi < 0` supports `shortening`.

APAURA also calculates a within-branch PAS-usage change magnitude:

```text
TVD = 0.5 * sum_k |PAU_case,k - PAU_control,k|
```

This is reported as `tutr_shift_magnitude` after branch weighting. If usage changes strongly but the net WPI does not move, APAURA records `apa_change=yes` and `tutr_direction=no_directional_shift`. It does not create an additional APA architecture category for that pattern.

## ALE

Terminal reads are normalized across terminal-exon branches:

```text
B_b = terminal reads in branch b / total terminal reads
```

The magnitude of branch switching is:

```text
ALE_switch = 0.5 * sum_b |B_case,b - B_control,b|
```

A branch-level weighted position index summarizes movement toward proximal or distal terminal branches:

```text
ALE_WPI = sum_b B_b * x_b
```

The resulting direction is reported as `proximal_switch`, `distal_switch`, or `switch` when branch usage changes but no robust positional direction can be assigned.

## IPA

For each sample, termination evidence is partitioned into terminal and non-terminal compartments. Intronic PAS usage is:

```text
IPA_usage = intronic PAS evidence / total retained gene-end evidence
```

The condition change is:

```text
delta_ipa_usage = mean(IPA_usage_case) - mean(IPA_usage_control)
```

A positive supported change is `ipa_direction=gain`; a negative supported change is `ipa_direction=loss`.

An IPA-only event is represented as:

```text
apa_architecture = none
ipa_direction    = gain
apa_direction    = gain
```

This keeps terminal architecture separate from intronic termination.

## Mixed architecture and complex changes

A gene is `mixed` only when both TUTR and ALE terminal architectures are present.

If more than one APA mechanism changes, `apa_direction=complex` while mechanism-specific fields retain the details. For example:

```text
apa_architecture     = mixed
architecture_components = TUTR+ALE
tutr_direction       = shortening
ale_direction        = distal_switch
apa_direction        = complex
changed_mechanisms   = ALE+TUTR
mechanism_directions = ALE=distal_switch;TUTR=shortening
```

IPA can co-occur with either TUTR or ALE without becoming part of `apa_architecture`; it is added to `changed_mechanisms` when it changes.

For example, a gene with TUTR shortening and IPA gain remains `apa_architecture=TUTR`, with `changed_mechanisms=IPA+TUTR` and `apa_direction=complex`.

## Termination-compartment shift

APAURA builds a condition-level usage vector:

```text
u = [terminal, IPA, internal_exon, downstream]
```

The total-variation distance is:

```text
termination_shift_index = 0.5 * sum_a |u_case,a - u_control,a|
```

and the Jensen-Shannon divergence is reported as `termination_jsd`.

These metrics summarize where transcript termination is occurring. They are deliberately named as termination metrics rather than terminal APA architecture metrics.

## Overall result interpretation

The recommended interpretation order is:

1. `apa_architecture` — TUTR, ALE, mixed, or none;
2. `apa_change` — did at least one supported mechanism change?;
3. `tutr_direction`, `ale_direction`, `ipa_direction` — which mechanism changed and how?;
4. `apa_direction` — concise overall direction;
5. quantitative effect sizes such as `delta_tutr_wpi`, `ale_switch_index`, `delta_ipa_usage`, and `termination_shift_index`.
