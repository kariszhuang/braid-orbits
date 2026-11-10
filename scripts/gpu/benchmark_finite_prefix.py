#!/usr/bin/env python3
"""Bounded shared-buffer CPU4/Metal screening comparison from a real finite case.
No production switch and no scientific negative certificate is generated here.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import struct
import subprocess

from verify_genus_three_finite import mul,inv,braid,smooth
PROJECT=Path(__file__).resolve().parents[2]


def packed(p):
    return sum(x<<(4*i) for i,x in enumerate(p))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('folder',type=Path);p.add_argument('out',type=Path)
    args=p.parse_args();root=args.out.resolve();root.mkdir();(root/'source').mkdir()
    inputs=['finite_prefix_cpu.rs','finite_prefix.swift','finite_prefix.metal','benchmark_finite_prefix.py','verify_genus_three_finite.py']
    for name in inputs:shutil.copyfile(_scripts_path.find(name),root/'source'/name)
    source=root/'source';frames=[]
    with (args.folder/'FRAMES.tsv').open()as f:
        for row in csv.DictReader(f,delimiter='\t'):
            a=tuple(map(int,row['a'].split(',')));b=tuple(map(int,row['b'].split(',')));ab=mul(a,b)
            factors=[mul(mul(inv(b),a),b),a,b,mul(mul(a,b),inv(a)),mul(mul(inv(ab),b),ab)]
            frames.append((a,b,factors))
    with (root/'FACTORS.u64le').open('xb')as out:
        for _,_,factors in frames:
            for f in factors:out.write(struct.pack('<Q',packed(f)))
    count=0
    with (args.folder/'COMPLETIONS.tsv').open()as f,(root/'COMPLETIONS.u64le').open('xb')as out:
        for row in csv.DictReader(f,delimiter='\t'):
            c=tuple(map(int,row['c'].split(',')));out.write(struct.pack('<Q',packed(c)))
            # Literal independent short-word checks, strided over actual input.
            if count%997==0:
                a,b,factors=frames[(count//997)%len(frames)];d=inv(mul(mul(a,b),c));t=(a,b,c,d)
                for factor,word in zip(factors,[(1,),(4,),(0,1),(0,4),(3,1)]):
                    u=t
                    for m in word:u=braid(u,m)
                    seen=set();x=0
                    while x not in seen:seen.add(x);x=factor[c[x]]
                    assert (len(seen)!=15)==(not smooth(u))
            count+=1
    subprocess.run(['rustc','--edition=2021','--crate-type','cdylib','-C','opt-level=3','-C','target-cpu=native',str(source/'finite_prefix_cpu.rs'),'-o',str(root/'cpu.dylib')],check=True)
    subprocess.run(['swiftc','-O',str(source/'finite_prefix.swift'),'-o',str(root/'gpu')],check=True)
    subprocess.run([str(root/'gpu'),str(root/'cpu.dylib'),str(source/'finite_prefix.metal'),str(root/'COMPLETIONS.u64le'),str(root/'FACTORS.u64le'),str(root/'RESULTS.json'),'5','4'],check=True)
    bindings={str(f.relative_to(root)):hashlib.sha256(f.read_bytes()).hexdigest() for f in root.rglob('*') if f.is_file()}
    (root/'BINDINGS.json').write_text(json.dumps(bindings,indent=2)+'\n')


if __name__=='__main__':main()
