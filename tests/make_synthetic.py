#!/usr/bin/env python3

from pathlib import Path

root = Path(__file__).resolve().parents[1]
ex = root / "examples"
ex.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------
# Synthetic direct-RNA / true 3' end dataset
#
# GENE1 contains two transcript ends in the same terminal-exon
# branch:
#
#   proximal PAS = 700
#   distal PAS   = 900
#
# Control favors distal usage.
# Treated favors proximal usage.
#
# Expected:
#   architecture = TUTR
#   direction    = shortening
#   delta WPI    = -0.55
# ------------------------------------------------------------

header = (
    "@HD\tVN:1.6\tSO:coordinate\n"
    "@SQ\tSN:chr1\tLN:2000\n"
)

def make(path, prox, dist):
    rows = []
    n = 0

    for end, count in [(700, prox), (900, dist)]:
        start = end - 49  # 1-based SAM start; 50M ends at 'end'

        for _ in range(count):
            n += 1
            rows.append(
                f"r{n}\t0\tchr1\t{start}\t60\t50M\t*\t0\t0\t"
                + ("A" * 50)
                + "\t"
                + ("I" * 50)
            )

    path.write_text(
        header + "\n".join(rows) + "\n"
    )


# Control: mostly distal PAS
make(ex / "control1.sam", 8, 32)
make(ex / "control2.sam", 10, 30)

# Treated: mostly proximal PAS
make(ex / "treated1.sam", 32, 8)
make(ex / "treated2.sam", 30, 10)


# ------------------------------------------------------------
# Sample sheet
# Paths are relative to examples/samples.tsv
# ------------------------------------------------------------

(ex / "samples.tsv").write_text(
    "sample\tcondition\tbam\n"
    "C1\tcontrol\tcontrol1.sam\n"
    "C2\tcontrol\tcontrol2.sam\n"
    "T1\ttreated\ttreated1.sam\n"
    "T2\ttreated\ttreated2.sam\n"
)


# ------------------------------------------------------------
# Minimal GTF
#
# Two transcripts share the same terminal-exon start (600),
# therefore they belong to the same terminal branch, while
# terminating at 700 and 900 respectively.
#
# This gives one TUTR branch containing two PAS.
# ------------------------------------------------------------

gtf = """chr1\tAPAURA\tgene\t100\t900\t.\t+\t.\tgene_id "GENE1"; gene_name "GENE1"; gene_type "protein_coding";
chr1\tAPAURA\ttranscript\t100\t700\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; gene_name "GENE1"; gene_type "protein_coding"; transcript_type "protein_coding";
chr1\tAPAURA\texon\t100\t199\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; gene_name "GENE1"; exon_number "1";
chr1\tAPAURA\texon\t600\t700\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; gene_name "GENE1"; exon_number "2";
chr1\tAPAURA\ttranscript\t100\t900\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX2"; gene_name "GENE1"; gene_type "protein_coding"; transcript_type "protein_coding";
chr1\tAPAURA\texon\t100\t199\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX2"; gene_name "GENE1"; exon_number "1";
chr1\tAPAURA\texon\t600\t900\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX2"; gene_name "GENE1"; exon_number "2";
"""

(ex / "mini.gtf").write_text(gtf)
