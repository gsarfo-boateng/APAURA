# RBP motif enrichment in APAURA

APAURA v0.10.1 adds an architecture-aware RBP motif-enrichment module. The purpose is to prioritize candidate RNA-binding proteins whose sequence motifs are enriched near PAS or ALE splice-acceptor regions that change between conditions.

This module is an independent implementation. It does not call or embed another RBP-analysis toolkit. It consumes generic MEME-format position-weight matrices (PWMs), so motif collections may come from any appropriately licensed source or from user-generated experiments.

## Biological strata

APAURA does not pool all regulated RNA ends into one set. Foreground sites are stratified by mechanism and direction:

- `TUTR:shortening`: significantly gained PAS in genes with TUTR shortening;
- `TUTR:lengthening`: significantly gained PAS in genes with TUTR lengthening;
- `ALE:proximal_switch`: gained PAS/terminal-exon branch associated with a proximal ALE switch;
- `ALE:distal_switch`: gained PAS/terminal-exon branch associated with a distal ALE switch;
- `IPA:gain`: significantly gained intronic PAS;
- `IPA:loss`: significantly lost intronic PAS.

Genes with mixed TUTR+ALE architecture can contribute to both mechanism-specific strata. This is intentional: regulatory interpretation follows the changed mechanism rather than forcing the whole gene into one pooled category.

A TUTR change with `no_directional_shift` is excluded by default and can be included with `--include-no-directional-shift`.

## Regulatory targets

### PAS windows

TUTR and IPA are scanned around the strand-aware PAS. ALE is also scanned around the PAS of the gaining terminal branch.

Default transcript-oriented window:

```text
-100 nt                          PAS                    +50 nt
  |------------------------------|-----------------------|
```

For a minus-strand gene the genomic sequence is reverse-complemented so all sequences are represented 5' to 3' in transcript orientation.

### ALE splice-acceptor windows

ALE additionally receives a terminal-exon 5' boundary window, which represents the splice-acceptor region for spliced alternative terminal exons.

Default window:

```text
-200 nt                    terminal-exon boundary        +100 nt
  |--------------------------------|------------------------|
```

This lets APAURA distinguish candidate factors enriched near terminal-exon splice choice from factors enriched near the terminal PAS.

## Motif scoring

MEME `letter-probability matrix` motifs are parsed as A/C/G/T probabilities. RNA U is normalized to T for genomic sequence scanning.

For motif position `i` and base `b`, APAURA uses the log2 odds contribution

\[
s_{ib}=\log_2\left(\frac{p_{ib}}{0.25}\right).
\]

The score of a sequence match is

\[
S=\sum_i s_{i,b_i}.
\]

To make thresholds comparable among motifs with different lengths and information contents, the raw score is normalized to the motif-specific theoretical score range:

\[
S_{norm}=\frac{S-S_{min}}{S_{max}-S_{min}}.
\]

The default hit threshold is

\[
S_{norm}\ge 0.80.
\]

Configure with `--motif-score-frac`.

APAURA scans only the transcript-oriented RNA sequence, not both DNA strands, because an RBP motif is interpreted in RNA orientation.

## Matched background

The preferred background comprises non-changing sites of the same mechanism/target type. Background sequences are selected for each foreground target using nearest GC content, excluding the same gene where possible.

`--background-ratio 3` requests up to three matched background targets per foreground target.

If no stable background sites exist for a stratum, APAURA falls back to non-foreground sites of the same mechanism while excluding foreground genes. The output records `background_mode` so this is transparent.

## Enrichment statistics

For each RBP motif, a site is considered motif-positive if at least one match passes the configured PWM threshold.

The enrichment table is

```text
                         motif +       motif -
regulated foreground        a             b
matched background          c             d
```

The effect size is an odds ratio with a 0.5 continuity correction for stability:

\[
OR=\frac{(a+0.5)(d+0.5)}{(b+0.5)(c+0.5)}.
\]

Because the biological question is enrichment, APAURA uses a one-sided Fisher exact test for increased motif prevalence in the foreground. Benjamini-Hochberg correction is applied across motif/stratum tests.

APAURA also performs label permutations while preserving foreground/background sample sizes. The empirical enrichment P value is

\[
p_{perm}=\frac{1+\#(D_{perm}\ge D_{obs})}{1+N_{perm}},
\]

where

\[
D=\Pr(motif+\mid foreground)-\Pr(motif+\mid background).
\]

Permutation P values receive a separate BH correction.

## Outputs

`rbp.tsv`
: One row per architecture/direction, target type and motif. Contains hit rates, odds ratio, Fisher P/FDR, permutation P/FDR and background mode.

`rbp_hits.tsv`
: Site-level motif matches, best normalized score and transcript-relative motif start positions.

`rbp_profile.tsv`
: Motif hit counts by transcript-relative start position. This supports positional enrichment plots around PAS or ALE acceptors.

`rbp_windows.tsv`
: Extracted foreground and candidate-background sequences with GC content.

`rbp_run_info.txt`
: Analysis parameters and test definitions.

## Example

```bash
java -jar APAURA.jar rbp-enrich \
  --sites apaura_results/apa_sites.tsv \
  --genes apaura_results/apa.tsv \
  --fasta genome.fa \
  --motifs rbp_motifs.meme \
  --upstream 100 \
  --downstream 50 \
  --splice-upstream 200 \
  --splice-downstream 100 \
  --site-fdr 0.05 \
  --min-delta-pau 0.10 \
  --motif-score-frac 0.80 \
  --background-ratio 3 \
  --permutations 1000 \
  --out rbp_results
```

An uncompressed FASTA with a samtools-compatible `.fai` index is recommended for large genomes. Without an index, APAURA can load an unindexed FASTA into memory; gzipped FASTA is therefore better reserved for small references/testing.

## Interpretation boundary

Motif enrichment identifies candidate regulatory RBPs. It does **not** establish binding or causality. Stronger support can come from CLIP/eCLIP overlap, RBP expression or activity changes, perturbation experiments, and replication in independent samples. These evidence layers can be integrated in future versions without changing the APA architecture model.
