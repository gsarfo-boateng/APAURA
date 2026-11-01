# Tutorial: protocol-aware 3′-end sequencing with `run-3seq`

This tutorial describes the APAURA `run-3seq` workflow for aligned 3′-end sequencing datasets.

`run-3seq` is intended for protocols in which an aligned read endpoint represents or closely tracks a transcript cleavage/polyadenylation endpoint. Examples include 3SEQ, PAS-seq, PolyA-seq, 3′READS/3′READS+, PAT-seq, A-seq/A-seq2, SAPAS, QuantSeq, and related 3′-end-enriched assays.

The current implementation is protocol-aware without claiming that every published protocol has been independently validated. Protocol presets primarily control endpoint geometry. Protocol-specific UMI parsing and FASTA-based internal-priming filtering are planned separately.

## 1. Requirements

You need:

- APAURA v0.10.1 or later with `run-3seq`;
- aligned BAM or SAM files for each biological sample;
- a matching GTF/GTF.GZ annotation;
- at least two biological replicates per condition for replicated gene-level moderated Dirichlet-multinomial inference;
- at least two retained PAS per gene for gene-level compositional testing.

Check the command:

```bash
java -jar APAURA.jar run-3seq --help
```

## 2. Sample sheet

Create a tab-separated `samples.tsv`:

```text
sample	condition	bam
WT_1	WT	/absolute/path/WT_1.bam
WT_2	WT	/absolute/path/WT_2.bam
KO_1	KO	/absolute/path/KO_1.bam
KO_2	KO	/absolute/path/KO_2.bam
```

Absolute paths are recommended.

## 3. Protocol preset

A conservative starting point is:

```text
generic
```

Current recognized names include:

```text
3seq
3p-seq
pas-seq
polyaseq
3reads
3reads+
pat-seq
a-seq
a-seq2
sapas
3t-fill
paperclip
2p-seq
pa-seq
wtts-seq
maps
3pc
exprss
quantseq-fwd
quantseq-rev
quantseq-ffpe
generic
```

The presence of a name does not imply that every chemistry-specific rule has been fully modeled or independently benchmarked.

Implemented endpoint behavior includes:

- `generic`: transcript-oriented 3′ alignment endpoint;
- `quantseq-rev`: alignment 5′ endpoint with strand flipping to recover the biological transcript end.

Not yet implemented:

- FASTA-based internal-priming filtering;
- protocol-specific UMI extraction/deduplication.

## 4. Run APAURA

```bash
java -jar APAURA.jar run-3seq \
  --samples /absolute/path/samples.tsv \
  --gtf /absolute/path/annotation.gtf \
  --protocol generic \
  --control WT \
  --case KO \
  --out /absolute/path/APAURA_3SEQ_WT_vs_KO
```

Useful options:

```text
--cluster-distance INT
--cluster-representative mode|median
--min-gene-reads INT
--min-site-reads INT
--min-pau FLOAT
--min-mapq INT
--keep-duplicates
--strand-mode MODE
--annotation-slop INT
--min-delta-pau FLOAT
--min-ale-switch FLOAT
--fdr FLOAT
--wpi-threshold FLOAT
--min-internal-exon-delta FLOAT
--min-ipa-delta FLOAT
--min-tutr-shift FLOAT
```

Current defaults include a PAS clustering distance of 20 nt, minimum pooled gene support of 10 reads per condition, minimum pooled PAS support of 5 reads, minimum pooled PAU of 0.01, minimum MAPQ of 10, and FDR threshold of 0.05.

## 5. Workflow

```text
aligned BAM/SAM
    ↓
protocol-aware endpoint extraction
    ↓
PAS clustering
    ↓
GTF annotation
    ↓
PAS counting by biological replicate
    ↓
PAU estimation
    ↓
TUTR / ALE / IPA architecture and mechanism analysis
    ↓
PAS-level differential usage
    ↓
gene-level moderated Dirichlet-multinomial testing
    ↓
canonical APAURA output
```

## 6. Poly(A)-site usage

```text
PAU_gks = C_gks / sum_j C_gjs
```

and:

```text
delta_PAU_k = mean(PAU_case,k) - mean(PAU_control,k)
```

Positive ΔPAU indicates increased usage in the case condition; negative ΔPAU indicates decreased usage.

## 7. Gene-level moderated Dirichlet-multinomial test

For a gene with `K` retained PAS, each sample contributes:

```text
C_gs = (C_g1s, C_g2s, ..., C_gKs)
```

The test asks whether the whole PAS composition changes:

```text
H0: control and case share the same PAS composition
H1: control and case have different PAS compositions
```

The likelihood-ratio statistic is:

```text
LR = 2 × (logLik_full - logLik_null)
```

APAURA applies Benjamini-Hochberg correction across evaluable gene-level tests.

The moderated implementation estimates a common precision parameter across testable genes to stabilize inference.

Current requirements:

- at least two retained PAS;
- at least two evaluable control replicates;
- at least two evaluable case replicates.

A single-PAS gene is reported as:

```text
insufficient_sites
```

## 8. Effect size and direction

The canonical gene-level effect size is APD:

```text
APD = 0.5 × sum_k |PAU_case,k - PAU_control,k|
```

APD is unsigned.

For tandem 3′ UTR APA:

```text
delta_WPI > 0   distal shift / lengthening
delta_WPI < 0   proximal shift / shortening
```

Minus-strand genes are handled in transcript orientation.

## 9. Output files

```text
apa.tsv
apa_events.tsv
apa_sites.tsv
apa_architecture.tsv
pas_all.bed
pas.bed
pas_counts.tsv
pau.tsv
run_info.txt
details/
```

Additional `run-3seq` diagnostic:

```text
details/three_seq_dm.tsv
```

### `apa.tsv`

Primary gene-level table:

```text
pvalue      = moderated DM gene-level p-value
fdr         = BH-adjusted DM FDR
effect_size = APD
```

### `apa_sites.tsv`

PAS-level table with PAU, ΔPAU, p-value, FDR, and status. Replicated PAS-level inference uses beta-binomial testing.

### `apa_events.tsv`

Mechanism/event table for TUTR, ALE, IPA, and supported combinations.

### `details/three_seq_dm.tsv`

```text
gene_id
n_sites
n_control_evaluable
n_case_evaluable
pvalue
fdr
converged
common_precision
status
```

## 10. Recommended interpretation

Start with:

```text
apa_change == yes
fdr <= 0.05
```

Then inspect:

- `effect_size` for magnitude;
- `apa_direction` for direction;
- `changed_mechanisms` and `primary_mechanism`;
- `apa_sites.tsv` for individual PAS;
- `apa_events.tsv` for mechanism-level detail.

## 11. Current limitations

1. FASTA-based internal-priming filtering is not yet implemented.
2. Protocol-specific UMI extraction and UMI-aware deduplication are not yet implemented.
3. Not every protocol preset has been independently benchmarked.
4. `run-3seq` expects aligned BAM/SAM input, not FASTQ.
5. Moderated gene-level DM requires biological replication and at least two retained PAS.
6. `apa.tsv` is the primary gene-level inference table; PAS-level significance localizes contributing sites.

## 12. Reproducibility

Preserve:

- APAURA version;
- exact command line;
- protocol preset;
- genome build;
- GTF release;
- sample sheet;
- filtering thresholds;
- `run_info.txt`;
- complete result directory.

## 13. Help

```bash
java -jar APAURA.jar run-3seq --help
```

See also:

```text
README.md
docs/ALGORITHM.md
docs/ARCHITECTURE_MODEL.md
docs/OUTPUT_SCHEMA.md
```
