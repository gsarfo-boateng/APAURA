from pathlib import Path
root=Path(__file__).resolve().parents[1]
out=root/'examples'/'rbp'; out.mkdir(parents=True, exist_ok=True)
seq=['G']*12000
# TUTR shortening truth
fg=[500,700,900,1100]
bg=[2000,2200,2400,2600,2800,3000,3200,3400]
for p in fg: seq[p-10:p-5]=list('AAAAA')
# ALE lengthening truth: motif is planted around the terminal-exon acceptor, not PAS.
ale_fg=[(5000,5100),(5400,5500),(5800,5900),(6200,6300)] # (acceptor, PAS)
ale_bg=[(7000+i*300,7100+i*300) for i in range(8)]
for a,p in ale_fg: seq[a-10:a-5]=list('CCCCC')
with open(out/'synthetic.fa','w') as f:
    f.write('>chrR\n'); s=''.join(seq)
    for i in range(0,len(s),80): f.write(s[i:i+80]+'\n')
with open(out/'motifs.meme','w') as f:
    f.write('MEME version 4\n\nALPHABET= ACGT\n\n')
    for motif,base in [('RBP_A',0),('RBP_C',1),('RBP_G',2)]:
        f.write(f'MOTIF {motif} {motif}\nletter-probability matrix: alength= 4 w= 5\n')
        for _ in range(5):
            row=[0.001]*4; row[base]=0.997
            f.write(' '.join(f'{v:.3f}' for v in row)+'\n')
        f.write('\n')
with open(out/'apa.tsv','w') as f:
    f.write('gene_id\tapa_architecture\tapa_change\tapa_direction\ttutr_direction\tale_direction\tipa_direction\n')
    for i,p in enumerate(fg,1): f.write(f'FG{i}\tTUTR\tyes\tshortening\tshortening\tnot_applicable\tnot_applicable\n')
    for i,p in enumerate(bg,1): f.write(f'BG{i}\tTUTR\tno\tstable\tstable\tnot_applicable\tnot_applicable\n')
    for i,(a,p) in enumerate(ale_fg,1): f.write(f'AF{i}\tALE\tyes\tlengthening\tnot_applicable\tlengthening\tnot_applicable\n')
    for i,(a,p) in enumerate(ale_bg,1): f.write(f'AB{i}\tALE\tno\tstable\tnot_applicable\tstable\tnot_applicable\n')
with open(out/'apa_sites.tsv','w') as f:
    f.write('site_id\tgene_id\tchrom\tpos0\tstrand\tregion\tarchitecture_role\tterminal_exon_anchor0\tsite_rank\tdelta_pau\tpvalue\tfdr\n')
    for i,p in enumerate(fg,1): f.write(f'FGS{i}\tFG{i}\tchrR\t{p}\t+\t3UTR\tTUTR\t400\t1\t0.50\t0.0001\t0.001\n')
    for i,p in enumerate(bg,1): f.write(f'BGS{i}\tBG{i}\tchrR\t{p}\t+\t3UTR\tTUTR\t1900\t1\t0.00\t1\t1\n')
    for i,(a,p) in enumerate(ale_fg,1): f.write(f'AFS{i}\tAF{i}\tchrR\t{p}\t+\t3UTR\tALE\t{a}\t2\t0.50\t0.0001\t0.001\n')
    for i,(a,p) in enumerate(ale_bg,1): f.write(f'ABS{i}\tAB{i}\tchrR\t{p}\t+\t3UTR\tALE\t{a}\t2\t0.00\t1\t1\n')
