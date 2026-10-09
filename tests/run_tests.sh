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
  --min-informative-replicates 2 \
  --write-diagnostics \
  --out tests/out_short_sam
awk -F'\t' '
  NR==1 {for(i=1;i<=NF;i++) c[$i]=i; next}
  $(c["apa_architecture"])=="TUTR" &&
  $(c["tutr_change"])=="yes" &&
  $(c["tutr_direction"])=="shortening" &&
  $(c["changed_mechanisms"])=="TUTR" {ok=1}
  END {exit !ok}
' tests/out_short_sam/apa.tsv
test -s tests/out_short_sam/apa_events.tsv
test -s tests/out_short_sam/apa_sites.tsv
test -s tests/out_short_sam/details/apa_architecture.tsv
test -s tests/out_short_sam/details/apa_full.tsv
test -s tests/out_short_sam/details/apa_events_full.tsv
test -s tests/out_short_sam/details/apa_sites_full.tsv
! head -1 tests/out_short_sam/apa.tsv | grep -q 'n_tutr_groups'
head -1 tests/out_short_sam/apa.tsv | grep -q 'min_mechanism_fdr'

# Same engine through the internal BAM parser
java -jar APAURA.jar short-read \
  --samples examples/shortread/samples_bam.tsv \
  --gtf examples/shortread/mini_shortread.gtf \
  --control control --case treated \
  --min-sample-evidence 10 \
  --min-informative-replicates 2 \
  --write-diagnostics \
  --out tests/out_short_bam
awk -F'\t' '
  NR==1 {for(i=1;i<=NF;i++) c[$i]=i; next}
  $(c["apa_architecture"])=="TUTR" &&
  $(c["tutr_change"])=="yes" &&
  $(c["tutr_direction"])=="shortening" &&
  $(c["changed_mechanisms"])=="TUTR" {ok=1}
  END {exit !ok}
' tests/out_short_bam/apa.tsv
cmp tests/out_short_sam/apa.tsv tests/out_short_bam/apa.tsv
cmp tests/out_short_sam/apa_events.tsv tests/out_short_bam/apa_events.tsv
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
grep -P -q '^IPA1\tIPA1\tprotein_coding\tSINGLE_PAS\tnone\tyes\tgain\tnot_applicable\tnot_applicable\tgain\tnot_applicable\tIPA:gain\tIPA\tIPA=gain' tests/out_architecture/apa.tsv
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

# Single-cell APA regression tests
rm -rf tests/sc_synthetic tests/out_sc_reference tests/out_sc_compare
mkdir -p tests/sc_synthetic

# Public CLI help
java -jar APAURA.jar sc --help > tests/sc_synthetic/sc_help.txt

grep -q 'gtf-sites' tests/sc_synthetic/sc_help.txt
grep -q 'import-sites' tests/sc_synthetic/sc_help.txt
grep -q 'merge-sites' tests/sc_synthetic/sc_help.txt
grep -q 'quantify' tests/sc_synthetic/sc_help.txt
grep -q 'reference' tests/sc_synthetic/sc_help.txt
grep -q 'compare' tests/sc_synthetic/sc_help.txt

for cmd in gtf-sites import-sites merge-sites quantify reference compare; do
  java -jar APAURA.jar sc "$cmd" --help >/dev/null
done

# Small transcript annotation with two PAS per gene.
cat > tests/sc_synthetic/mini.gtf <<'EOF'
chr1	test	gene	401	600	.	+	.	gene_id "GENE1"; gene_name "GENE1";
chr1	test	transcript	401	500	.	+	.	gene_id "GENE1"; transcript_id "GENE1_TX1"; gene_name "GENE1";
chr1	test	exon	401	500	.	+	.	gene_id "GENE1"; transcript_id "GENE1_TX1"; gene_name "GENE1";
chr1	test	transcript	401	600	.	+	.	gene_id "GENE1"; transcript_id "GENE1_TX2"; gene_name "GENE1";
chr1	test	exon	401	600	.	+	.	gene_id "GENE1"; transcript_id "GENE1_TX2"; gene_name "GENE1";
chr2	test	gene	1100	1300	.	-	.	gene_id "GENE2"; gene_name "GENE2";
chr2	test	transcript	1200	1300	.	-	.	gene_id "GENE2"; transcript_id "GENE2_TX1"; gene_name "GENE2";
chr2	test	exon	1200	1300	.	-	.	gene_id "GENE2"; transcript_id "GENE2_TX1"; gene_name "GENE2";
chr2	test	transcript	1100	1300	.	-	.	gene_id "GENE2"; transcript_id "GENE2_TX2"; gene_name "GENE2";
chr2	test	exon	1100	1300	.	-	.	gene_id "GENE2"; transcript_id "GENE2_TX2"; gene_name "GENE2";
EOF

# PAS catalog plus architecture sidecars.
java -jar APAURA.jar sc gtf-sites \
  --annotation tests/sc_synthetic/mini.gtf \
  --out tests/sc_synthetic/sites.tsv

test -s tests/sc_synthetic/sites.tsv
test -s tests/sc_synthetic/sites.site_architecture.tsv
test -s tests/sc_synthetic/sites.gene_architecture.tsv
test -s tests/sc_synthetic/sites.architecture_summary.tsv

# UMI-deduplication test: read1/read2 are the same molecule.
cat > tests/sc_synthetic/reads.sam <<'EOF'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:5000
@SQ	SN:chr2	LN:5000
read1	0	chr1	451	60	50M	*	0	0	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	*	CB:Z:CELL1	UB:Z:UMI1	GX:Z:GENE1
read2	0	chr1	451	60	50M	*	0	0	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	*	CB:Z:CELL1	UB:Z:UMI1	GX:Z:GENE1
read3	0	chr1	551	60	50M	*	0	0	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	*	CB:Z:CELL1	UB:Z:UMI2	GX:Z:GENE1
read4	0	chr1	551	60	50M	*	0	0	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	*	CB:Z:CELL2	UB:Z:UMI3	GX:Z:GENE1
read5	16	chr2	1200	60	50M	*	0	0	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	*	CB:Z:CELL1	UB:Z:UMI4	GX:Z:GENE2
read6	16	chr2	1100	60	50M	*	0	0	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	*	CB:Z:CELL2	UB:Z:UMI5	GX:Z:GENE2
EOF

java -jar APAURA.jar sc quantify \
  --sites tests/sc_synthetic/sites.tsv \
  --alignment tests/sc_synthetic/reads.sam \
  --sample-id TEST_SAMPLE \
  --gene-tag GX \
  --three-prime-by gene \
  --site-window 50 \
  --out tests/sc_synthetic/counts.tsv

test -s tests/sc_synthetic/counts.tsv

# Six alignments contain five unique UMIs.
awk -F'\t' '
  NR > 1 {total += $4}
  END {exit !(total == 5)}
' tests/sc_synthetic/counts.tsv

# Grouped reference generation.
cat > tests/sc_synthetic/cells.tsv <<'EOF'
barcode	sample_id	cell_type
CELL1	TEST_SAMPLE	TestCells
CELL2	TEST_SAMPLE	TestCells
EOF

java -jar APAURA.jar sc reference \
  --sites tests/sc_synthetic/sites.tsv \
  --counts tests/sc_synthetic/counts.tsv \
  --cells tests/sc_synthetic/cells.tsv \
  --group-by cell_type \
  --out tests/out_sc_reference

test -n "$(find tests/out_sc_reference -type f -print -quit)"

# Descriptive control-versus-treated comparison.
cat > tests/sc_synthetic/samples.tsv <<'EOF'
sample_id	condition
CTRL	control
TRT	treated
EOF

cat > tests/sc_synthetic/compare_counts.tsv <<'EOF'
sample_id	group	site_id	count
CTRL	ALL	GENE1_PAS1	30
CTRL	ALL	GENE1_PAS2	5
CTRL	ALL	GENE2_PAS1	5
CTRL	ALL	GENE2_PAS2	30
TRT	ALL	GENE1_PAS1	5
TRT	ALL	GENE1_PAS2	30
TRT	ALL	GENE2_PAS1	30
TRT	ALL	GENE2_PAS2	5
EOF

java -jar APAURA.jar sc compare \
  --sites tests/sc_synthetic/sites.tsv \
  --samples tests/sc_synthetic/samples.tsv \
  --counts tests/sc_synthetic/compare_counts.tsv \
  --control control \
  --treated treated \
  --no-interactions \
  --out tests/out_sc_compare

test -s tests/out_sc_compare/all_results.tsv
test -s tests/out_sc_compare/all_results.with_architecture.tsv
test -s tests/out_sc_compare/observed_unreplicated_events.tsv
test -s tests/out_sc_compare/observed_unreplicated_events.with_architecture.tsv
test -s tests/out_sc_compare/apaura_usage.tsv
test -s tests/out_sc_compare/apaura_site_usage.tsv

grep -q 'OBSERVED_' tests/out_sc_compare/all_results.tsv
grep -q 'apa_architecture' tests/out_sc_compare/all_results.with_architecture.tsv
grep -q $'\tapd\tusage_jsd\t' tests/out_sc_compare/apaura_usage.tsv
grep -q $'\tcontrol_pau\ttreated_pau\tdelta_pau' tests/out_sc_compare/apaura_site_usage.tsv
grep -P -q '^ALL\tGENE1\tGENE1\tchr1\t\+\t2\tUNREPLICATED_OBSERVED_SHIFT\t1\t1\t35\t35\t0\.142857\t0\.857143\t0\.714286\t0\.714286\t.*\tyes\tno\tno\tstable\tlengthening\tstable\tobserved_lengthening\tdescriptive_only\t' tests/out_sc_compare/apaura_usage.tsv
grep -P -q '^ALL\tGENE2\tGENE2\tchr2\t-\t2\tUNREPLICATED_OBSERVED_SHIFT\t1\t1\t35\t35\t0\.857143\t0\.142857\t-0\.714286\t0\.714286\t.*\tyes\tno\tno\tstable\tshortening\tstable\tobserved_shortening\tdescriptive_only\t' tests/out_sc_compare/apaura_usage.tsv

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
