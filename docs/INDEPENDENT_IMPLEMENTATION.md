# Independent implementation statement

APAURA is intended to be an independent software implementation.

## What was not done

The project does not vendor, translate, copy or embed source files from APALORD, LABRAT, LAPA, DaPars/DaPars2, QAPA, REPAC, APAlyzer, InPACT or other APA analysis packages.

## What can legitimately be shared between methods

Many statistical and computational concepts used in bioinformatics are general mathematical methods rather than package-specific inventions, including:

- Fisher's exact test;
- beta-binomial models;
- negative-binomial models;
- empirical cumulative distributions;
- Kolmogorov-Smirnov-type distribution distances;
- Benjamini-Hochberg false-discovery-rate correction;
- isotonic regression;
- Jensen-Shannon divergence;
- total-variation distance.

APAURA implements the methods it uses directly in Java and documents their role in the analysis.

## What defines APAURA scientifically

APAURA is not intended to be novel because it uses a particular significance test. Its intended novelty is the **hierarchical architecture representation** and the separation of evidence extraction from architecture inference.

The same architecture engine can therefore be fed by future sequencing technologies without redefining the biological output model.

## Reproducibility policy

Tests use programmatically generated synthetic data with known truth. They are designed to verify behavior, not to reproduce proprietary or package-specific test fixtures.

## RBP module

The v0.10 RBP motif module is implemented directly in Java within APAURA. It parses generic MEME position-probability matrices, performs transcript-oriented sequence extraction, computes PWM log-odds scores, constructs architecture-aware foreground/background sets, and performs enrichment statistics internally. It does not wrap, translate, or vendor another RBP motif-scanning program.
