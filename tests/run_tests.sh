#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

python3 tests/make_synthetic.py
python3 tests/make_bam.py
python3 tests/make_shortread_synthetic.py
python3 tests/make_shortread_bam.py
python3 tests/make_isoquant_synthetic.py
python3 tests/make_architecture_synthetic.py
./build.sh

rm -rf tests/out_direct tests/out_isoquant tests/out_short_sam tests/out_short_bam

# Direct-RNA / true 3'-end regression test
java -jar APAURA.jar run-direct \
  --samples examples/samples.tsv \
  --gtf examples/mini.gtf \
  --control control --case treated \
  --cluster-distance 0 --min-total 1 \
  --min-gene-reads 10 --min-site-reads 5 --min-pau 0.01 \
  --out tests/out_direct
grep -q $'TUTR\tTUTR\tyes\tshortening\tshortening' tests/out_direct/apa.tsv
grep -P -q '^GENE1\t.*\tTUTR\tTUTR' tests/out_direct/apa.tsv
grep -q $'0.550000\t-0.550000' tests/out_direct/apa.tsv
grep -q 'ks_fdr' tests/out_direct/apa.tsv
java -jar APAURA.jar ends --bam tests/synthetic.bam --out tests/synthetic.ends.tsv
grep -q $'chr1\t699\t700' tests/synthetic.ends.tsv
grep -q $'chr1\t899\t900' tests/synthetic.ends.tsv

# IsoQuant corrected_reads BED should reproduce the direct-end result.
java -jar APAURA.jar run-isoquant \
  --samples examples/isoquant/samples.tsv \
  --gtf examples/mini.gtf \
  --control control --case treated \
  --cluster-distance 0 \
  --out tests/out_isoquant
cmp tests/out_direct/apa.tsv tests/out_isoquant/apa.tsv
cmp tests/out_direct/apa_sites.tsv tests/out_isoquant/apa_sites.tsv

# Conventional short-read coverage engine: SAM
java -jar APAURA.jar short-read \
  --samples examples/shortread/samples.tsv \
  --gtf examples/shortread/mini_shortread.gtf \
  --control control --case treated \
  --min-sample-evidence 10 \
  --out tests/out_short_sam
grep -q $'TUTR\tyes\tshortening\tTUTR:shortening' tests/out_short_sam/apa.tsv

# Same engine through the internal BAM parser
java -jar APAURA.jar short-read \
  --samples examples/shortread/samples_bam.tsv \
  --gtf examples/shortread/mini_shortread.gtf \
  --control control --case treated \
  --min-sample-evidence 10 \
  --out tests/out_short_bam
grep -q $'TUTR\tyes\tshortening\tTUTR:shortening' tests/out_short_bam/apa.tsv
cmp tests/out_short_sam/apa.tsv tests/out_short_bam/apa.tsv
cmp tests/out_short_sam/apa_sites.tsv tests/out_short_bam/apa_sites.tsv


# Hierarchical APA architecture truth set: TUTR, ALE, IPA, and mixed TUTR+ALE.
rm -rf tests/out_architecture
java -jar APAURA.jar run-direct \
  --samples examples/architecture/samples.tsv \
  --gtf examples/architecture/architecture.gtf \
  --control control --case treated \
  --cluster-distance 0 --min-gene-reads 10 --min-site-reads 5 --min-pau 0.01 \
  --out tests/out_architecture
grep -P -q '^TUTR1\tTUTR1\tprotein_coding\tTUTR\tTUTR\tyes\tshortening\tshortening\tnot_applicable\tnot_applicable\tnot_applicable\tTUTR:shortening\tTUTR\tTUTR=shortening' tests/out_architecture/apa.tsv
grep -P -q '^ALE1\tALE1\tprotein_coding\tALE\tALE\tyes\tdistal_switch\tnot_applicable\tdistal_switch\tnot_applicable\tnot_applicable\tALE:distal_switch\tALE\tALE=distal_switch' tests/out_architecture/apa.tsv
grep -P -q '^IPA1\tIPA1\tprotein_coding\tnone\tnone\tyes\tgain\tnot_applicable\tnot_applicable\tgain\tnot_applicable\tIPA:gain\tIPA\tIPA=gain' tests/out_architecture/apa.tsv
grep -P -q '^MIX1\tMIX1\tprotein_coding\tmixed\tTUTR\+ALE\tyes\tcomplex\tshortening\tdistal_switch\tnot_applicable\tnot_applicable\tmixed:complex\tALE\+TUTR\tALE=distal_switch;TUTR=shortening' tests/out_architecture/apa.tsv
grep -P -q '^NODIR1\tNODIR1\tprotein_coding\tTUTR\tTUTR\tyes\tno_directional_shift\tno_directional_shift\tnot_applicable\tnot_applicable\tnot_applicable\tTUTR:no_directional_shift\tTUTR\tTUTR=no_directional_shift' tests/out_architecture/apa.tsv
grep -P -q '^TUTRIPA1\tTUTRIPA1\tprotein_coding\tTUTR\tTUTR\tyes\tcomplex\tshortening\tnot_applicable\tgain\tnot_applicable\tmulti_mechanism:complex\tIPA\+TUTR\tIPA=gain;TUTR=shortening' tests/out_architecture/apa.tsv
grep -q 'apa_direction' tests/out_architecture/apa_architecture.tsv
grep -q 'tutr_direction' tests/out_architecture/apa_architecture.tsv
grep -q 'ale_direction' tests/out_architecture/apa_architecture.tsv
grep -q 'ipa_direction' tests/out_architecture/apa_architecture.tsv
grep -q 'apa_event' tests/out_architecture/apa_architecture.tsv
grep -q 'termination_shift_index' tests/out_architecture/apa_architecture.tsv
grep -q 'ale_switch_index' tests/out_architecture/apa_architecture.tsv
grep -q 'delta_tutr_wpi' tests/out_architecture/apa_architecture.tsv

# RBP motif-enrichment regression test
python3 tests/make_rbp_synthetic.py
rm -rf tests/out_rbp
java -jar APAURA.jar rbp-enrich \
  --sites examples/rbp/apa_sites.tsv \
  --genes examples/rbp/apa.tsv \
  --fasta examples/rbp/synthetic.fa \
  --motifs examples/rbp/motifs.meme \
  --upstream 20 --downstream 10 \
  --background-ratio 2 --permutations 200 --seed 7 \
  --out tests/out_rbp
grep -P -q '^TUTR:shortening\tPAS\tRBP_A\tRBP_A\t4\t8\t4\t0\t' tests/out_rbp/rbp.tsv
grep -P -q '^ALE:distal_switch\tALE_acceptor\tRBP_C\tRBP_C\t4\t8\t4\t0\t' tests/out_rbp/rbp.tsv
grep -q 'relative_hit_starts' tests/out_rbp/rbp_hits.tsv
grep -q $'RBP_A\tRBP_A\t-10\t4' tests/out_rbp/rbp_profile.tsv

echo 'All APAURA v0.10.1 synthetic tests passed.'
