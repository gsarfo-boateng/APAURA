#!/usr/bin/env python3
from pathlib import Path
root=Path(__file__).resolve().parents[1]
ex=root/'examples'
header='@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:2000\n'
# Two PAS: last aligned base 700 (proximal) or 900 (distal), 50M on + strand.
def make(path, prox, dist):
    rows=[]; n=0
    for end,count in [(700,prox),(900,dist)]:
        start=end-49 # 1-based SAM start so 50M ends at 'end'
        for i in range(count):
            n+=1; rows.append(f'r{n}\t0\tchr1\t{start}\t60\t50M\t*\t0\t0\t'+('A'*50)+'\t'+('I'*50))
    path.write_text(header+'\n'.join(rows)+'\n')
make(ex/'control1.sam',8,32)
make(ex/'control2.sam',10,30)
make(ex/'treated1.sam',32,8)
make(ex/'treated2.sam',30,10)
