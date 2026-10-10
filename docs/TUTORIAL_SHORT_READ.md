# Tutorial: conventional short-read RNA-seq

## Evidence model

Conventional fragmented RNA-seq does not provide direct biochemical cleavage-site evidence. APAURA therefore does **not** treat ordinary short-read alignment ends as cleavage sites.

The short-read engine instead uses annotation-aware, CIGAR-aware coverage/evidence across terminal structures and annotated PAS groups.

Typical processing:

1. load transcript models;
2. reconcile terminal-exon families;
3. construct annotated tandem APA groups;
4. construct annotated IPA candidates;
5. quantify per-sample evidence;
6. estimate PAS usage;
7. perform replicate-aware differential testing;
8. summarize positional movement with WPI and redistribution with APD.

## Sample sheet

```text
sample\tcondition\tbam
WT_1\tWT\t/path/WT_1.bam
WT_2\tWT\t/path/WT_2.bam
WT_3\tWT\t/path/WT_3.bam
MUT_1\tMUT\t/path/MUT_1.bam
MUT_2\tMUT\t/path/MUT_2.bam
MUT_3\tMUT\t/path/MUT_3.bam
```

Accepted path column names include `bam`, `file`, and `path`.

## Run

```bash
java -jar APAURA.jar short-read \
  --samples samples.tsv \
  --gtf annotation.gtf \
  --control WT \
  --case MUT \
  --library reverse \
  --out apaura_short
```

## Important options

```text
--library unstranded|forward|reverse
--min-mapq INT
--keep-duplicates
--terminal-anchor-tolerance INT
--pas-merge-distance INT
--ipa-window INT
--min-segment-length INT
--min-sample-evidence INT
--min-informative-replicates INT
--min-delta-pau FLOAT
--min-ipa-delta FLOAT
--fdr FLOAT
--wpi-threshold FLOAT
--write-per-sample
--write-diagnostics
--no-ipa
```

Use the installed JAR help for the authoritative list:

```bash
java -jar APAURA.jar short-read --help
```

## Replicate-aware inference

For replicated studies, APAURA tracks:

- total biological replicates;
- replicates with any evidence;
- informative replicates.

Supported replicated events use beta-binomial testing. Genuine 1-vs-1 experiments can use Fisher's exact test. Replicated experiments with too few informative replicates are reported as `insufficient_replicate_support` rather than being rescued by pooled Fisher testing.

BH correction is applied only across valid finite p-values.

## Interpreting `apa_sites.tsv`

Key fields:

```text
gene_name
pas_position_1based
site_rank
normalized_position
control_pau
case_pau
delta_pau
control_pau_sd
case_pau_sd
control_n_informative
case_n_informative
pvalue
fdr
method
status
```

```text
delta_pau > 0 -> higher relative use in case
delta_pau < 0 -> lower relative use in case
```

## Interpreting `apa.tsv`

Key fields:

```text
apa_architecture
n_sites
control_wpi
case_wpi
delta_wpi
apd
min_fdr
apa_change
observed_direction
apa_direction
apa_event
status
```

```text
delta_wpi > 0 -> distal movement
delta_wpi < 0 -> proximal movement
```

`observed_direction` preserves the observed positional effect. `apa_direction` is the inferential/classified direction and can be `no_directional_shift` when redistribution is supported but net WPI movement is small.

## Conservative filtering

A simple high-confidence terminal set is:

```text
status == ok
apa_change == yes
```

Then stratify by `apa_architecture`, `apa_direction`, and `apa_event`.

Do not treat `insufficient_replicate_support` as a biological negative; it means current biological replication is insufficient for the inferential call.
