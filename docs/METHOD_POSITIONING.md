# Method positioning and intended gap

## Existing method space

The APA software ecosystem already contains strong methods for several individual problems:

- DaPars/DaPars2: de novo/coverage-based 3'UTR APA from conventional RNA-seq.
- QAPA: annotated 3'UTR isoform/PAS usage from transcript abundance.
- LABRAT: differential APA with TUTR/ALE/mixed structural labels.
- REPAC: compositional analysis of annotated PAS usage from conventional RNA-seq.
- LAPA: PAS calling from diverse 3'-end and long-read evidence with statistical testing.
- APALORD: differential APA from long-read RNA-seq.
- APAlyzer: 3'UTR and intronic polyadenylation analysis.
- InPACT/IPAFinder: dedicated intronic-polyadenylation analysis.

Therefore, APAURA should not claim novelty merely because it accepts BAM files, calculates PAU, uses a beta-binomial distribution, calculates a KS-like statistic, or labels TUTR/ALE genes.

## Intended gap

APAURA instead targets a narrower integration problem:

> Differential APA can occur at several nested structural levels within the same gene, but many workflows summarize the gene primarily as a site-level change, a single proximal/distal score, or one categorical structure label.

APAURA treats APA as a hierarchy:

```text
termination compartment
    -> terminal-exon branch
        -> PAS within branch
```

This allows independent quantification of:

1. terminal-to-intronic remodeling / IPA;
2. ALE branch switching;
3. TUTR shortening or lengthening within a branch;
4. combinations of these processes in the same gene.

## Why this should remain useful

The architecture model is downstream of sequencing technology. New evidence types only need an adapter that provides strand-aware terminal/PAS evidence. The biological representation remains unchanged.

This is intentionally different from building the tool around one sequencing platform or one specific external quantification package.

## Claims that should wait until benchmarking

Do not claim that APAURA is more accurate, faster or more sensitive than existing methods until real and simulated benchmarking supports those claims.

For v1.0 validation, compare at least:

- PAS coordinate concordance;
- PAU correlation;
- TUTR direction agreement;
- ALE branch-switch agreement;
- IPA event recovery;
- false-positive rate on null simulations;
- sensitivity across depth and replicate number;
- runtime and memory;
- robustness to annotation complexity.
