# APAURA algorithm map

APAURA separates **evidence extraction** from **APA architecture inference**.

## Long-read / true 3'-end evidence

```text
BAM/SAM or IsoQuant corrected BED
  -> transcript-oriented single-base 3' ends
  -> optional PAS microheterogeneity clustering
  -> GTF structural annotation
  -> PAS support filtering
  -> per-sample PAU
  -> site-level differential usage
  -> hierarchical APA architecture analysis
```

The architecture analysis separately measures:

```text
termination compartment
ALE branch choice
TUTR PAS choice within terminal branch
```

See `LONG_READ_ALGORITHM.md` and `ARCHITECTURE_MODEL.md`.

## Conventional short-read RNA-seq

Random-fragment RNA-seq endpoints are not treated as PAS. APAURA uses CIGAR-aware coverage across annotated tandem terminal-exon segments, fits a monotonic cumulative coverage profile and reconstructs annotated TUTR isoform usage.

The workflow also outputs an annotation-derived catalog of potential TUTR/ALE/mixed/IPA-compatible gene structures.

See `SHORT_READ_ALGORITHM.md`.

## Statistical layer

- beta-binomial LRT for replicated PAS usage;
- Fisher exact fallback;
- Benjamini-Hochberg FDR;
- gene-level discrete ECDF maximum-distance diagnostic;
- total-variation distance and Jensen-Shannon divergence for architecture changes;
- WPI/APD effect sizes.

The project uses independent Java implementations of standard mathematical methods. See `INDEPENDENT_IMPLEMENTATION.md`.
