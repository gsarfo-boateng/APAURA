# Changelog

## 0.10.1

### Added

- `rbp-enrich` command for architecture-aware RBP motif enrichment.
- Generic MEME PWM parser and Java log2-odds motif scanner.
- Transcript-oriented FASTA extraction with reverse-complement handling on negative-strand genes.
- `.fai` random-access FASTA support for large uncompressed genomes, with unindexed FASTA fallback.
- Mechanism-specific RBP strata for TUTR shortening/lengthening, ALE proximal/distal switching, and IPA gain/loss.
- ALE analysis at both terminal PAS and terminal-exon splice-acceptor windows.
- GC-matched stable-site background selection with transparent non-foreground fallback.
- One-sided Fisher enrichment testing, BH FDR, label-permutation enrichment P values and permutation FDR.
- Site-level motif-hit output and transcript-relative positional motif profiles.
- Synthetic RBP-enrichment truth test with a planted motif.

### Preserved

- Long-read/direct-RNA PAS analysis.
- IsoQuant corrected-read support.
- Conventional STAR/BWA-style coverage-based TUTR analysis.
- TUTR/ALE/mixed terminal architecture model and separate IPA mechanism.

### Maintenance

- Reformatted the Java source into conventional readable style without changing analysis behavior.
- Added descriptive helper methods and technical Javadocs around statistical, BAM/BGZF, and APA-analysis code.
- Added JUnit tests for core statistical and parsing routines.
- Added development documentation and stricter ignore rules for generated/macOS metadata files.

### Interpretation

The RBP module prioritizes candidate regulators from sequence enrichment. Motif evidence is not presented as proof of RBP binding or causality.
