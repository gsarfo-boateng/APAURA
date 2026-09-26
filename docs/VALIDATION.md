# APAURA validation

APAURA v0.10.1 includes unit tests and synthetic regression tests for statistical routines, evidence extraction, differential APA calculations, architecture classification, and RBP enrichment.

Run JUnit tests with:

```bash
mvn test
```

Run the end-to-end synthetic regression suite with:

```bash
./tests/run_tests.sh
```


## Core unit tests

JUnit tests under `src/test/java/org/apaura/` cover the mathematical and parsing primitives used by the workflows, including Fisher exact testing, Benjamini-Hochberg correction, isotonic regression, the tie-aware ECDF statistic, beta-binomial differential usage, PAS clustering, and CIGAR parsing.

These tests complement rather than replace the end-to-end synthetic workflows.

## Direct RNA / true 3' end test

The direct-end synthetic test contains a two-PAS TUTR gene in which control samples favor the distal PAS and treated samples favor the proximal PAS.

Expected result:

```text
apa_architecture = TUTR
apa_change       = yes
apa_direction    = shortening
tutr_direction   = shortening
```

The same molecule-end truth is tested through SAM, BAM, and IsoQuant corrected-read BED input.

## Conventional short-read test

A synthetic terminal-exon coverage dataset is analyzed through both SAM and BAM parsing. The expected result is TUTR shortening, and the two input formats must produce identical summary and site-level outputs.

## Architecture/mechanism truth set

The architecture test contains six genes:

| Gene | Synthetic truth | Expected result |
|---|---|---|
| TUTR1 | two PAS in one terminal branch, proximal shift in treatment | `TUTR / shortening` |
| ALE1 | two distinct terminal branches, treatment favors distal branch | `ALE / distal_switch` |
| IPA1 | canonical terminal PAS plus intronic PAS, treatment gains IPA | `apa_architecture=none`, `ipa_direction=gain` |
| MIX1 | TUTR branch plus an alternative terminal branch; both mechanisms change | `mixed / complex` |
| NODIR1 | three PAS in one TUTR branch; PAS usage changes but net WPI is unchanged | `apa_change=yes`, `no_directional_shift` |
| TUTRIPA1 | TUTR shortening plus intronic PAS gain | `apa_architecture=TUTR`, `ipa_direction=gain`, overall `complex` |

For MIX1, APAURA must retain both mechanism-specific directions:

```text
apa_architecture     = mixed
apa_direction        = complex
ale_direction        = distal_switch
tutr_direction       = shortening
changed_mechanisms   = ALE+TUTR
mechanism_directions = ALE=distal_switch;TUTR=shortening
```

NODIR1 tests an important boundary condition. Its PAS distribution changes strongly (`tutr_shift_magnitude` is large), but `delta_tutr_wpi` is zero. APAURA therefore records a supported APA change without inventing a shortening/lengthening direction:

```text
apa_change       = yes
tutr_direction   = no_directional_shift
apa_direction    = no_directional_shift
```

## What synthetic testing establishes

The suite verifies implementation behavior against known truth. It does not establish biological sensitivity/specificity in real datasets or fully calibrate all P values. Real-data benchmarking remains necessary before v1.0.

## Architecture-versus-IPA separation test

`TUTRIPA1` is specifically included to prevent IPA from leaking into the terminal architecture ontology. The expected result is:

```text
apa_architecture     = TUTR
architecture_components = TUTR
tutr_direction       = shortening
ipa_direction        = gain
apa_direction        = complex
changed_mechanisms   = IPA+TUTR
apa_event            = multi_mechanism:complex
```

If this gene is ever reported as `apa_architecture=IPA` or `apa_architecture=mixed`, the architecture separation has regressed.

## RBP enrichment validation (v0.10.1)

The RBP module includes a deterministic synthetic truth dataset in `examples/rbp/`.

Two independent regulatory patterns are planted in the reference sequence:

1. `RBP_A` motif at all four foreground PAS associated with `TUTR:shortening`, absent from eight stable GC-matched TUTR PAS.
2. `RBP_C` motif at all four ALE terminal-exon acceptor windows associated with `ALE:distal_switch`, absent from eight stable ALE acceptor windows. The same motif is deliberately absent from the corresponding ALE PAS windows.

The regression suite requires APAURA to recover both enrichments and their expected target types. A ubiquitous `RBP_G` control motif is present in both foreground and background and must not be called enriched.

This validates sequence extraction, target stratification, ALE acceptor handling, PWM scanning, stable-site background selection, Fisher enrichment and positional-hit output. It does not substitute for benchmarking with experimentally measured RBP binding.
