#!/usr/bin/env python3
from pathlib import Path
root=Path(__file__).resolve().parents[1]
base=root/'examples'/'isoquant'
base.mkdir(parents=True,exist_ok=True)
counts={
    'control1':(8,32),
    'control2':(10,30),
    'treated1':(32,8),
    'treated2':(30,10),
}
for sample,(prox,dist) in counts.items():
    d=base/sample
    d.mkdir(parents=True,exist_ok=True)
    rows=[]; n=0
    for end,c in [(700,prox),(900,dist)]:
        start=end-50
        for _ in range(c):
            n+=1
            rows.append(f'chr1\t{start}\t{end}\t{sample}_r{n}\t0\t+')
    (d/'OUT.corrected_reads.bed').write_text('\n'.join(rows)+'\n')
(base/'samples.tsv').write_text(
    'sample\tcondition\tisoquant_dir\n'
    'control1\tcontrol\tcontrol1\n'
    'control2\tcontrol\tcontrol2\n'
    'treated1\ttreated\ttreated1\n'
    'treated2\ttreated\ttreated2\n'
)
