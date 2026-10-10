# Validation

APAURA v0.10.1 is validated with unit tests, synthetic end-to-end regression
tests, CLI/package checks, and real-data acceptance analyses across multiple RNA
sequencing modalities.

## Required build and regression suite

After every code change:

```bash
mvn clean package
./tests/run_tests.sh
```

The v0.10.1 finalization audit compiled 56 production source files and 20 test
source files and ran 46 JUnit tests with zero failures, errors, or skips. The
end-to-end synthetic suite also completed successfully.

Tests cover direct-RNA/SAM and BAM parsing, IsoQuant, short-read TUTR/ALE/IPA,
architecture, usage metrics, statistics, single-cell integration, and RBP motif
enrichment.

## Statistical safeguard

PAS-level inference follows the experimental design rather than sequencing depth
alone.

```text
true experimental 1-vs-1 + evaluable 1-vs-1
    -> Fisher-exact-unreplicated

replicated design + >=2 evaluable samples in each condition
    -> beta-binomial-LRT

replicated design + <2 evaluable samples in either condition
    -> not_tested
    -> insufficient_replicate_support
```

A replicated event must never fall back to Fisher simply because only one sample
per group is evaluable.

## Real-data acceptance analyses

### Conventional short-read RNA-seq

A CLL `U1_MUT` versus `U1_WT` 10-vs-10 analysis produced 30,268 gene-level
results and 122 APA changes.

The validation audit found:

```text
mechanism-support leakage = 0
APA-change inconsistencies = 0
```

This dataset exercises replicated short-read TUTR, ALE, IPA, and
multi-mechanism inference.

### IsoQuant long-read

Mouse METTL3 IsoQuant analyses were used to test the reconciled terminal
structure model.

Post-resolution-fix acceptance:

```text
CLL comparison:
  genes = 8,840
  sites = 21,096
  significant APA changes = 20
  structural violations = 0
  cross-mechanism support leakage = 0

NM comparison:
  genes = 10,233
  sites = 27,207
  significant APA changes = 51
  structural violations = 0
  cross-mechanism support leakage = 0
```

Previously ambiguous retained terminal PAS moved to `UNRESOLVED` rather than
being forced into TUTR/ALE classes.

### Oxford Nanopore direct RNA

A `w1118` versus `FR113N` 3-vs-3 direct-RNA analysis produced:

```text
genes = 9,531
sites = 19,590
significant APA changes = 46
beta-binomial-LRT sites = 19,570
not_tested sites = 20
Fisher-exact-unreplicated sites = 0
```

After the Fisher safeguard fix, the two sparse `dhd` sites that had only one
evaluable sample per condition changed from Fisher testing to:

```text
method = not_tested
status = insufficient_replicate_support
```

All 46 significant gene-level APA calls were preserved.

Final acceptance:

```text
strict structural violations = 0
cross-mechanism support leakage = 0
APA-change consistency flags = 0
```

### PacBio HiFi technical stress test

A genuine 1-vs-1 PacBio comparison exercised the opposite side of the Fisher
rule:

```text
genes = 5,720
sites = 16,253
Fisher-exact-unreplicated = 16,253
beta-binomial-LRT = 0
not_tested = 0

structural violations = 0
mechanism-support leakage = 0
APA-change consistency flags = 0
TUTR-direction flags = 0
```

This run is retained as a **technical/statistical stress test**, not as a
biological benchmark. The FLAIR BAMs were generated from parent CCS molecules
rather than the segmented/refined FLNC transcript molecules, so the resulting
IPA-heavy biological calls are not interpreted as validated PacBio biology.

### Single-cell

The integrated single-cell workflow preserves the sc3UTR core and adds
harmonized APAURA PAU/TVD/JSD/WPI summaries.

The underlying single-cell workflow had already been benchmarked against
SCAPE, Sierra, scUTRquant, and SCAPTURE on the same dataset. v0.10.1 regression
testing therefore focuses on preserving the established sc3UTR behavior while
testing the APAURA usage overlay.

Formal single-cell inference uses the **biological sample** as the replicate
unit. Cells from the same sample are not treated as independent biological
replicates.

## Structural acceptance criteria

For retained terminal PAS:

```text
one resolved terminal group + one PAS      -> SINGLE_PAS
one resolved terminal group + >=2 PAS      -> TUTR
multiple groups + one PAS per group        -> ALE
multiple groups + tandem PAS in any group  -> mixed
unassigned retained terminal PAS           -> UNRESOLVED
```

IPA remains a separate termination mechanism.

## Mechanism-specific support acceptance

A mechanism can enter `changed_mechanisms` only when both its effect threshold
and mechanism-specific statistical support are satisfied.

```text
TUTR          -> significant TUTR or ALE+TUTR site/event support
ALE           -> significant ALE or ALE+TUTR site/event support
IPA           -> significant IPA support
internal_exon -> significant internal_exon support
```

A gene-wide minimum FDR from another mechanism is not sufficient.

## CLI and packaging smoke tests

At minimum:

```bash
java -jar target/APAURA.jar --version
java -jar target/APAURA.jar --citation
java -jar target/APAURA.jar --help
java -jar target/APAURA.jar run-direct --help
java -jar target/APAURA.jar run-isoquant --help
java -jar target/APAURA.jar short-read --help
java -jar target/APAURA.jar sc --help
java -jar target/APAURA.jar plot-region --help
java -jar target/APAURA.jar rbp-enrich --help
```

## Plot regression

For locus plots, inspect at least one plus-strand and one minus-strand gene and
confirm:

- condition colors remain separate;
- genomic coordinates align across panels;
- minus-strand mirroring is display-only;
- PAS/IPA coordinates match result tables;
- cache hits are reused where expected;
- multi-gene plotting writes one plot per target.
