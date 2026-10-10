# Library strandedness

APAURA accepts:

```text
--library unstranded
--library forward
--library reverse
```

For a typical paired-end library:

- forward-stranded: R1 follows transcript/sense orientation and R2 the opposite orientation;
- reverse-stranded: R1 is antisense to the transcript and R2 is sense.

Many dUTP/TruSeq stranded protocols are reverse-stranded, but protocol documentation or an empirical strandedness check should be used rather than guessing.

For stranded paired-end `plot-region`, APAURA selects R1/R2 orientation combinations consistent with the target transcript strand before calculating locus coverage.

Genomic coordinates remain left-to-right. Negative-strand coverage may be mirrored vertically as a visualization convention.
