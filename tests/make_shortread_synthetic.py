#!/usr/bin/env python3
from pathlib import Path
root=Path(__file__).resolve().parents[1]
ex=root/'examples'/'shortread'
ex.mkdir(parents=True,exist_ok=True)

gtf='''chr1\ttest\texon\t101\t200\t.\t+\t.\tgene_id "GAPA1"; transcript_id "GAPA1_short"; gene_name "ApaDemo"; gene_type "protein_coding"; transcript_type "protein_coding";\nchr1\ttest\tCDS\t101\t120\t.\t+\t0\tgene_id "GAPA1"; transcript_id "GAPA1_short"; gene_name "ApaDemo"; gene_type "protein_coding"; transcript_type "protein_coding";\nchr1\ttest\texon\t101\t300\t.\t+\t.\tgene_id "GAPA1"; transcript_id "GAPA1_long"; gene_name "ApaDemo"; gene_type "protein_coding"; transcript_type "protein_coding";\nchr1\ttest\tCDS\t101\t120\t.\t+\t0\tgene_id "GAPA1"; transcript_id "GAPA1_long"; gene_name "ApaDemo"; gene_type "protein_coding"; transcript_type "protein_coding";\n'''
(ex/'mini_shortread.gtf').write_text(gtf)
header='@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n'

def make(path, common_reads, distal_reads):
    rows=[]
    # Segment 0 is [100,200); segment 1 is [200,300). Each alignment is 10M
    # and lies wholly inside one segment. Number of aligned bases therefore
    # gives exact mean depth: common_reads*10/100 versus distal_reads*10/100.
    idx=0
    for i in range(common_reads):
        pos0=100+(i%91); idx+=1
        rows.append((pos0,f'c{idx}\t0\tchr1\t{pos0+1}\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII'))
    for i in range(distal_reads):
        pos0=200+(i%91); idx+=1
        rows.append((pos0,f'd{idx}\t0\tchr1\t{pos0+1}\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII'))
    rows.sort(key=lambda x:x[0])
    path.write_text(header+'\n'.join(x[1] for x in rows)+'\n')

make(ex/'control1.sam',100,75)
make(ex/'control2.sam',100,72)
make(ex/'treated1.sam',100,25)
make(ex/'treated2.sam',100,28)
(ex/'samples.tsv').write_text('sample\tcondition\tbam\ncontrol1\tcontrol\tcontrol1.sam\ncontrol2\tcontrol\tcontrol2.sam\ntreated1\ttreated\ttreated1.sam\ntreated2\ttreated\ttreated2.sam\n')
print(ex)
