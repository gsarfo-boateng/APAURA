# Related APA tools and APAURA positioning

APAURA is not intended to reproduce an existing APA package. Existing software already solves important portions of the APA problem.

| Tool | Representative strength | Why APAURA should not duplicate it |
|---|---|---|
| DaPars/DaPars2 | coverage-based de novo 3'UTR APA from conventional RNA-seq | established coverage-breakpoint approach |
| QAPA | annotated 3'UTR isoform usage from transcript abundance | established annotation/quantification workflow |
| LABRAT | differential APA and TUTR/ALE/mixed structural classification | already demonstrates value of structure-aware APA labels |
| REPAC | compositional annotated-PAS analysis from conventional RNA-seq | strong conventional RNA-seq compositional framework |
| LAPA | PAS calling across long-read/3'-end data and differential testing | already shows cross-source PAS analysis is possible |
| APALORD | differential APA from long-read RNA-seq | strong long-read-specific comparator |
| APAlyzer | 3'UTR and intronic polyadenylation analysis | explicit 3'UTR/IPA analysis already exists |
| InPACT / IPAFinder | dedicated IPA inference | specialized IPA benchmarks |

## APAURA's intended distinction

The project focuses on **nested APA architecture within one gene**:

```text
termination compartment
 -> terminal-exon branch
    -> PAS within branch
```

This produces separate quantitative outputs for:

- terminal-versus-intronic remodeling;
- ALE branch switching;
- TUTR shifts within a branch;
- simultaneous multi-mechanism remodeling.

Cross-platform evidence is useful, but it is not by itself the novelty claim.

See `METHOD_POSITIONING.md` and `ARCHITECTURE_MODEL.md`.
