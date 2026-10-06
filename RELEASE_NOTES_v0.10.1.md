# APAURA v0.10.1

**Alternative Polyadenylation Architecture, Usage, and Regulatory Analysis**

APAURA v0.10.1 is the first public validated source release.

## Highlights

- Direct-RNA / true 3′-end APA analysis from SAM/BAM.
- IsoQuant corrected-read APA analysis.
- Conventional short-read RNA-seq APA analysis.
- Transcript-aware TUTR, ALE, mixed terminal architecture, and independent IPA analysis.
- Mechanism-specific direction calls for TUTR, ALE, and IPA.
- Replicate-aware beta-binomial differential testing with 1-vs-1 Fisher fallback.
- RBP motif enrichment around regulated PAS and ALE splice-acceptor windows.
- Publication-style `plot-region` visualization.
- Standalone `plot-region` mode using BAM/SAM + GTF without a completed APAURA result directory.
- Multi-gene plotting, whole-gene overview, replicate tracks, strand-aware terminal zoom, and PDF/PNG/SVG output.
- Synthetic SAM/BAM, IsoQuant, architecture, short-read, and RBP regression tests.

## Validation

```text
All APAURA v0.10.1 synthetic tests passed.
```

Synthetic tests verify expected implementation behavior against constructed truth. Real-data benchmarking remains an ongoing validation layer before v1.0.

## Citation

> Boateng-Sarfo G. APAURA: Alternative Polyadenylation Architecture, Usage, and Regulatory Analysis. Version 0.10.1, 2026.
