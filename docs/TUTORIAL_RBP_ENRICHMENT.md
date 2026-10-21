# Tutorial: RBP motif enrichment

## Inputs

```text
apa_sites.tsv
apa.tsv
reference genome FASTA
RBP motifs in MEME format
```

## Run

```bash
apaura rbp-enrich \
  --sites results/apa_sites.tsv \
  --genes results/apa.tsv \
  --fasta genome.fa \
  --motifs rbp_motifs.meme \
  --out rbp_results
```

Common options:

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

Documented defaults include 100 nt PAS-upstream, 50 nt PAS-downstream, 200 nt ALE-acceptor upstream, 100 nt downstream, site FDR 0.05, minimum |delta_PAU| 0.10, motif score fraction 0.80, background ratio 3, 1000 permutations, and seed 20260830.

Outputs:

```text
rbp.tsv
rbp_hits.tsv
rbp_profile.tsv
rbp_windows.tsv
rbp_run_info.txt
```

Motif enrichment is associative. It does not establish physical binding or causality.
