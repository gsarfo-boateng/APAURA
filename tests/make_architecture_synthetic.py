#!/usr/bin/env python3
from pathlib import Path
root=Path(__file__).resolve().parents[1]
out=root/'examples'/'architecture'
out.mkdir(parents=True,exist_ok=True)

gtf=[]
def add(chrom,feature,s,e,strand,gid,tid,gname):
    attrs=f'gene_id "{gid}"; transcript_id "{tid}"; gene_name "{gname}"; gene_type "protein_coding";'
    gtf.append(f'{chrom}\ttest\t{feature}\t{s}\t{e}\t.\t{strand}\t.\t{attrs}')
# TUTR1: same terminal-exon 5' anchor, two annotated 3' ends
for tid,end in [('TUTR1_t1',700),('TUTR1_t2',900)]:
    add('chr1','exon',100,300,'+','TUTR1',tid,'TUTR1'); add('chr1','CDS',120,280,'+','TUTR1',tid,'TUTR1')
    add('chr1','exon',500,end,'+','TUTR1',tid,'TUTR1'); add('chr1','CDS',500,600,'+','TUTR1',tid,'TUTR1')
# ALE1: distinct terminal-exon branches
for tid,ts,te in [('ALE1_t1',2200,2400),('ALE1_t2',2600,2800)]:
    add('chr1','exon',2000,2100,'+','ALE1',tid,'ALE1'); add('chr1','CDS',2020,2080,'+','ALE1',tid,'ALE1')
    add('chr1','exon',ts,te,'+','ALE1',tid,'ALE1'); add('chr1','CDS',ts,min(ts+40,te),'+','ALE1',tid,'ALE1')
# IPA1: canonical terminal end plus an observed intronic cleavage site at 4300
add('chr1','exon',4000,4100,'+','IPA1','IPA1_t1','IPA1'); add('chr1','CDS',4020,4080,'+','IPA1','IPA1_t1','IPA1')
add('chr1','exon',4500,5000,'+','IPA1','IPA1_t1','IPA1'); add('chr1','CDS',4500,4600,'+','IPA1','IPA1_t1','IPA1')
# NODIR1: three PAS in one terminal branch; PAS-usage change with unchanged net WPI
for tid,end in [('NODIR1_t1',8500),('NODIR1_t2',8700),('NODIR1_t3',8900)]:
    add('chr1','exon',8000,8100,'+','NODIR1',tid,'NODIR1'); add('chr1','CDS',8020,8080,'+','NODIR1',tid,'NODIR1')
    add('chr1','exon',8200,end,'+','NODIR1',tid,'NODIR1'); add('chr1','CDS',8200,8300,'+','NODIR1',tid,'NODIR1')
# TUTRIPA1: TUTR terminal architecture plus an intronic PAS; both change in treatment
for tid,end in [('TUTRIPA1_t1',9700),('TUTRIPA1_t2',9900)]:
    add('chr1','exon',9000,9100,'+','TUTRIPA1',tid,'TUTRIPA1'); add('chr1','CDS',9020,9080,'+','TUTRIPA1',tid,'TUTRIPA1')
    add('chr1','exon',9400,end,'+','TUTRIPA1',tid,'TUTRIPA1'); add('chr1','CDS',9400,9500,'+','TUTRIPA1',tid,'TUTRIPA1')
# MIX1: TUTR within one branch + distinct ALE branch
for tid,end in [('MIX1_t1',6500),('MIX1_t2',6700)]:
    add('chr1','exon',6000,6100,'+','MIX1',tid,'MIX1'); add('chr1','CDS',6020,6080,'+','MIX1',tid,'MIX1')
    add('chr1','exon',6200,end,'+','MIX1',tid,'MIX1'); add('chr1','CDS',6200,6300,'+','MIX1',tid,'MIX1')
add('chr1','exon',6000,6100,'+','MIX1','MIX1_t3','MIX1'); add('chr1','CDS',6020,6080,'+','MIX1','MIX1_t3','MIX1')
add('chr1','exon',7000,7300,'+','MIX1','MIX1_t3','MIX1'); add('chr1','CDS',7000,7100,'+','MIX1','MIX1_t3','MIX1')
(out/'architecture.gtf').write_text('\n'.join(gtf)+'\n')

header='@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:10000\n'
def make(path,counts):
    rows=[];n=0
    for end,count in sorted(counts.items()):
        start=end-49
        for _ in range(count):
            n+=1;seq='A'*50;qual='I'*50
            rows.append(f'r{n}\t0\tchr1\t{start}\t60\t50M\t*\t0\t0\t{seq}\t{qual}')
    path.write_text(header+'\n'.join(rows)+'\n')
# 1-based aligned end coordinates
samples={
 'C1':{700:10,900:30,2400:30,2800:10,4300:5,5000:35,6500:5,6700:25,7300:10,8500:20,8700:0,8900:20,9250:3,9700:10,9900:30},
 'C2':{700:12,900:28,2400:28,2800:12,4300:6,5000:34,6500:6,6700:24,7300:10,8500:20,8700:0,8900:20,9250:3,9700:10,9900:30},
 'T1':{700:30,900:10,2400:10,2800:30,4300:30,5000:10,6500:20,6700:5,7300:25,8500:0,8700:40,8900:0,9250:25,9700:30,9900:10},
 'T2':{700:28,900:12,2400:12,2800:28,4300:28,5000:12,6500:18,6700:7,7300:25,8500:0,8700:40,8900:0,9250:25,9700:30,9900:10},
}
for name,c in samples.items():make(out/f'{name}.sam',c)
(out/'samples.tsv').write_text('sample\tcondition\tfile\nC1\tcontrol\tC1.sam\nC2\tcontrol\tC2.sam\nT1\ttreated\tT1.sam\nT2\ttreated\tT2.sam\n')
print(out)
