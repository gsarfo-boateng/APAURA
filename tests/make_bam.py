#!/usr/bin/env python3
from pathlib import Path
import struct,zlib,binascii
root=Path(__file__).resolve().parents[1]

def bam_record(name,pos0,flag=0,mapq=60,length=50):
    rn=(name+'\0').encode()
    cigar=struct.pack('<I',(length<<4)|0)
    seq=bytes([0x11]*((length+1)//2))
    qual=bytes([30]*length)
    bin_mq_nl=(0<<16)|(mapq<<8)|len(rn)
    flag_nc=(flag<<16)|1
    body=struct.pack('<iiIIiiii',0,pos0,bin_mq_nl,flag_nc,length,-1,-1,0)+rn+cigar+seq+qual
    return struct.pack('<i',len(body))+body

def bgzf_block(data):
    co=zlib.compressobj(level=6,wbits=-15)
    comp=co.compress(data)+co.flush()
    # gzip fixed header with FEXTRA + XLEN=6 + BC subfield
    total=18+len(comp)+8
    bsize=total-1
    hdr=b'\x1f\x8b\x08\x04'+b'\x00\x00\x00\x00'+b'\x00\xff'+struct.pack('<H',6)+b'BC'+struct.pack('<H',2)+struct.pack('<H',bsize)
    footer=struct.pack('<II',binascii.crc32(data)&0xffffffff,len(data)&0xffffffff)
    return hdr+comp+footer

header=b'@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:2000\n'
raw=b'BAM\x01'+struct.pack('<i',len(header))+header+struct.pack('<i',1)+struct.pack('<i',5)+b'chr1\0'+struct.pack('<i',2000)
# read ending at 700 (SAM last aligned base): pos0 = 650, 50M => last base0=699
raw+=bam_record('r1',650)
# read ending at 900: pos0=850 => last base0=899
raw+=bam_record('r2',850)
# standard BGZF EOF block
EOF_BLOCK=bytes.fromhex('1f8b08040000000000ff0600424302001b0003000000000000000000')
out=root/'tests'/'synthetic.bam'
out.write_bytes(bgzf_block(raw)+EOF_BLOCK)
print(out)
