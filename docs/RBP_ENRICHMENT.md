# RBP motif enrichment

`rbp-enrich` tests supplied RNA-binding-protein sequence motifs around
APAURA-regulated poly(A) sites and, for ALE events, terminal-exon splice
acceptors.

The analysis identifies **candidate regulatory associations**. Motif enrichment
does not by itself demonstrate physical RBP binding or causal regulation.

## Required inputs

```bash
apaura rbp-enrich   --sites results/apa_sites.tsv   --genes results/apa.tsv   --fasta genome.fa   --motifs rbp_motifs.meme   --out rbp_results
```

Required files:

- `apa_sites.tsv` from APAURA;
- `apa.tsv` from the same comparison;
- the matching reference genome FASTA;
- RBP motifs in MEME format.

The FASTA and APAURA results must use compatible chromosome/contig names.

## Default windows and thresholds

```text
PAS upstream window            100 nt
PAS downstream window           50 nt
ALE acceptor upstream window   200 nt
ALE acceptor downstream window 100 nt
site FDR threshold             0.05
minimum |delta_PAU|            0.10
motif score fraction           0.80
background ratio               3
permutations                   1000
random seed                    20260830
```

CLI options:

```text
--upstream INT
--downstream INT
--splice-upstream INT
--splice-downstream INT
--site-fdr FLOAT
--min-delta-pau FLOAT
--motif-score-frac FLOAT
--background-ratio INT
--permutations INT
--seed INT
--include-no-directional-shift
```

## Foreground interpretation

Foreground sites are selected from regulated APAURA results using the requested
site-FDR and effect thresholds. Directional APA classes can be used to separate
candidate sequence associations around different forms of remodeling.

`--include-no-directional-shift` allows supported redistribution events without
a net proximal/distal directional shift to enter the analysis.

For ALE, the analysis also evaluates sequence windows around the corresponding
terminal-exon splice acceptor context.

## Outputs

```text
rbp.tsv
rbp_hits.tsv
rbp_profile.tsv
rbp_windows.tsv
rbp_run_info.txt
```

- `rbp.tsv` — motif-level enrichment summary;
- `rbp_hits.tsv` — individual motif-hit records;
- `rbp_profile.tsv` — positional motif profile;
- `rbp_windows.tsv` — sequence-window metadata used for the analysis;
- `rbp_run_info.txt` — run parameters and provenance.

## Recommended interpretation

Use RBP enrichment to prioritize candidate regulators for follow-up. Strong
results are most useful when they are recurrent across comparisons, supported by
the relevant APA mechanism, and consistent with independent RBP expression,
binding, perturbation, or literature evidence.

Do not describe a motif-enrichment result alone as proof that an RBP binds the
site or causes the APA event.
