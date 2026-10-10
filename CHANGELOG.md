# Changelog

## [0.10.1] - 2026-10-06

- Added `SINGLE_PAS` and `UNRESOLVED` terminal-architecture states and reconciled terminal exon groups for long-read architecture assignment.
- Harmonized PAU, WPI, TVD/APD, and Jensen-Shannon usage summaries across bulk modalities and the integrated single-cell workflow.
- Restricted Fisher exact testing to genuine experimental 1-vs-1 designs; sparse replicated designs remain `not_tested` with `insufficient_replicate_support`.
- Enforced mechanism-specific significance support for TUTR, ALE, IPA, and internal-exon changes.
- Added real-data acceptance validation across short-read RNA-seq, IsoQuant, ONT direct RNA, and a PacBio statistical stress test, while retaining the existing single-cell cross-tool benchmark.


### Added

- Annotation-guided conventional short-read TUTR, ALE, mixed, and IPA analysis.
- Mechanism-specific APA direction reporting and multi-mechanism integration.
- RBP motif enrichment around PAS and ALE splice-acceptor windows.
- Publication-oriented `plot-region` command.
- Multi-gene, gene-file, region, whole-gene-overview, replicate-track, and strand-aware plotting.
- Standalone BAM/SAM + GTF locus visualization without a completed `--results` directory.
- Synthetic regression coverage for direct RNA, BAM/SAM, IsoQuant, short-read TUTR/ALE/IPA, mixed mechanisms, and RBP enrichment.

### Validation

- All APAURA v0.10.1 synthetic tests pass.
- Standalone plotting validated on a real HNRNPA2B1 short-read locus.

## Unreleased

### Plotting

- Strand-aware short-read locus coverage for forward/reverse paired-end libraries.
- Genomic x coordinates preserved for negative-strand genes.
- Optional negative-strand coverage mirroring as a display convention.
- Persistent `plot_normalization.tsv` cache.
- Multiple-gene plotting through comma-separated `--gene` and `--genes-file`.
- Condition-specific replicate faceting.
- Optional transcript-level gene models.

### Documentation

- Reworked README and tutorials.
- Added validation and release checklists.

## 0.10.1 - 2026

### Short-read APA

- Replicate-aware site testing.
- Minimum informative-replicate support.
- No pooled Fisher fallback for sparse replicated experiments.
- Fisher retained for genuine 1-vs-1 comparisons.
- Explicit `insufficient_replicate_support` and `insufficient_evidence` statuses.
- Separate observed and inferential APA direction.
- PAU/WPI SD fields.
- Compact default outputs with optional per-sample/diagnostic tables.

### IPA

- Replicate-support accounting in short-read IPA testing.
- Compact IPA event output with replicate-support fields.

### RBP

- Motif enrichment around regulated PAS and ALE splice acceptors.
