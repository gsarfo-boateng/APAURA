#!/usr/bin/env python3
from pathlib import Path
import struct,zlib,binascii,re
root=Path(__file__).resolve().parents[1]
ex=root/'examples'/'shortread'

CIGAR_OP={'M':0,'I':1,'D':2,'N':3,'S':4,'H':5,'P':6,'=':7,'X':8}

def bgzf_block(data):
    co=zlib.compressobj(level=6,wbits=-15); comp=co.compress(data)+co.flush()
    total=18+len(comp)+8; bsize=total-1
    hdr=b'\x1f\x8b\x08\x04'+b'\x00\x00\x00\x00'+b'\x00\xff'+struct.pack('<H',6)+b'BC'+struct.pack('<H',2)+struct.pack('<H',bsize)
    footer=struct.pack('<II',binascii.crc32(data)&0xffffffff,len(data)&0xffffffff)
    return hdr+comp+footer

EOF_BLOCK=bytes.fromhex('1f8b08040000000000ff0600424302001b0003000000000000000000')

def cigar_words(cigar):
    return [((int(n)<<4)|CIGAR_OP[op]) for n,op in re.findall(r'(\d+)([MIDNSHP=X])',cigar)]

def encode_seq(seq):
    code={'=':0,'A':1,'C':2,'M':3,'G':4,'R':5,'S':6,'V':7,'T':8,'W':9,'Y':10,'H':11,'K':12,'D':13,'B':14,'N':15}
    out=bytearray()
    for i in range(0,len(seq),2):
        a=code.get(seq[i].upper(),15); b=code.get(seq[i+1].upper(),15) if i+1<len(seq) else 0
        out.append((a<<4)|b)
    return bytes(out)

def make_bam(sam_path,bam_path):
    refs=[('chr1',1000)]
    text='@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n'.encode()
    raw=b'BAM\x01'+struct.pack('<i',len(text))+text+struct.pack('<i',len(refs))
    for name,ln in refs:
        nb=(name+'\0').encode(); raw+=struct.pack('<i',len(nb))+nb+struct.pack('<i',ln)
    for line in sam_path.read_text().splitlines():
        if not line or line.startswith('@'): continue
        f=line.split('\t'); qname,flag,rname,pos,mapq,cigar=f[:6]; seq=f[9]; qual=f[10]
        flag=int(flag); pos0=int(pos)-1; mapq=int(mapq); words=cigar_words(cigar); rn=(qname+'\0').encode(); lseq=len(seq)
        bin_mq_nl=(0<<16)|(mapq<<8)|len(rn); flag_nc=(flag<<16)|len(words)
        body=struct.pack('<iiIIiiii',0,pos0,bin_mq_nl,flag_nc,lseq,-1,-1,0)+rn
        body+=b''.join(struct.pack('<I',x) for x in words)+encode_seq(seq)+bytes(max(0,ord(ch)-33) for ch in qual)
        raw+=struct.pack('<i',len(body))+body
    bam_path.write_bytes(bgzf_block(raw)+EOF_BLOCK)

for stem in ['control1','control2','treated1','treated2']:
    make_bam(ex/f'{stem}.sam',ex/f'{stem}.bam')
(ex/'samples_bam.tsv').write_text('sample\tcondition\tbam\ncontrol1\tcontrol\tcontrol1.bam\ncontrol2\tcontrol\tcontrol2.bam\ntreated1\ttreated\ttreated1.bam\ntreated2\ttreated\ttreated2.bam\n')
print(ex/'samples_bam.tsv')
