#!/usr/bin/env python3
"""Finite exact GPU-screen benchmark. Does not modify or stop production jobs."""
import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path
import random
import statistics
import struct
import subprocess
import tempfile

PROJECT=Path(__file__).resolve().parents[2]


def oracle():
    d=4;input_bytes=bytearray();outputs=[]
    def index(p):
        seen=set();cycles=0
        for x in range(d):
            if x in seen:continue
            cycles+=1
            while x not in seen:seen.add(x);x=p[x]
        return d-cycles
    def transitive(a,b):
        seen={0};todo=[0]
        while todo:
            x=todo.pop()
            for y in (a[x],b[x]):
                if y not in seen:seen.add(y);todo.append(y)
        return len(seen)==d
    permutations=list(itertools.permutations(range(d)))
    for a,b,c in itertools.product(permutations,repeat=3):
        fourth=[0]*d
        for x in range(d):fourth[a[b[c[x]]]]=x
        t=[a,b,c,fourth];mask=0
        for bit,(i,j) in enumerate(itertools.combinations(range(4),2)):
            mask|=int(transitive(t[i],t[j]))<<bit
        outputs.append(sum(index(p)<<(6*j) for j,p in enumerate(t)) | mask<<24)
        input_bytes.extend(a+b+c)
    return bytes(input_bytes),struct.pack(f'<{len(outputs)}I',*outputs)


def summary(rows):
    groups={}
    for row in rows:groups.setdefault((row['degree'],row['family'],row['batch'],row['mode']),[]).append(row)
    results=[]
    for key,samples in groups.items():
        ratios=[float(r['cpu_seconds'])/float(r['gpu_wall_seconds']) for r in samples]
        rng=random.Random(20260929);bootstrap=sorted(statistics.median(rng.choices(ratios,k=len(ratios))) for _ in range(3000))
        item=dict(zip(('degree','family','batch','mode'),map(int,key)),screen_speedup_median=statistics.median(ratios),
                  screen_ratio_bootstrap_low=bootstrap[75],screen_ratio_bootstrap_high=bootstrap[2924],
                  cpu_seconds_median=statistics.median(float(r['cpu_seconds']) for r in samples),
                  gpu_wall_seconds_median=statistics.median(float(r['gpu_wall_seconds']) for r in samples),
                  gpu_device_seconds_median=statistics.median(float(r['gpu_device_seconds']) for r in samples))
        if key[3]=='1':
            item.update(pipeline_speedup_median=statistics.median(float(r['cpu_pipeline_seconds'])/float(r['gpu_pipeline_seconds']) for r in samples),
                        cpu_pipeline_seconds_median=statistics.median(float(r['cpu_pipeline_seconds']) for r in samples),
                        gpu_pipeline_seconds_median=statistics.median(float(r['gpu_pipeline_seconds']) for r in samples),
                        survivors_min=min(int(r['survivors']) for r in samples),survivors_max=max(int(r['survivors']) for r in samples))
        results.append(item)
    return results


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repeats',type=int,default=7);parser.add_argument('--out',type=Path)
    args=parser.parse_args()
    if not 3<=args.repeats<=31:parser.error('repeats must be 3..31')
    out=(args.out or Path(tempfile.mkdtemp(prefix='gpu_screen_',dir=PROJECT/'runs/benchmarks'))).resolve()
    out.mkdir(parents=True,exist_ok=True)
    source_paths=[PROJECT/f'scripts/{f}' for f in ('gpu_screen_cpu.rs','gpu_screen.swift','gpu_screen.metal','benchmark_gpu.py')]+[PROJECT/'src/main.rs',PROJECT/'src/smart.rs',PROJECT/'src/hybrid.rs']
    metadata=dict(source_hashes={str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths},
                  cpu_workers=1,repeats=args.repeats,families={'0':'uniform triples','1':'matching triples','2':'known-failure neighborhoods'},
                  modes={'0':'four indices and all six pair flags','1':'necessary weak filter, short-circuiting'},
                  caveat='Technical repeats under concurrent production load; not independent experiments or full-search speedups.')
    (out/'METADATA.json').write_text(json.dumps(metadata,indent=2)+'\n')
    data,expected=oracle();(out/'ORACLE_INPUT.bin').write_bytes(data);(out/'ORACLE_OUTPUT.bin').write_bytes(expected)
    commands=[['rustc','--crate-type','cdylib','--edition','2021','-C','opt-level=3','-C','target-cpu=native','-C','codegen-units=1','-C','panic=abort',PROJECT/'scripts/gpu/gpu_screen_cpu.rs','-o',out/'libscreen.dylib'],
              ['xcrun','swiftc','-O','-whole-module-optimization',PROJECT/'scripts/gpu/gpu_screen.swift','-o',out/'gpu_screen']]
    with (out/'BUILD.log').open('w') as log:
        for command in commands:subprocess.run(list(map(str,command)),cwd=PROJECT,stdout=log,stderr=log,check=True)
    command=[out/'gpu_screen',out/'libscreen.dylib',PROJECT/'scripts/gpu/gpu_screen.metal',out/'MEASUREMENTS.tsv',out/'ORACLE_INPUT.bin',out/'ORACLE_OUTPUT.bin',args.repeats]
    (out/'COMMAND.json').write_text(json.dumps(list(map(str,command)),indent=2)+'\n')
    print(out,flush=True)
    with (out/'RUN.log').open('w') as log:subprocess.run(list(map(str,command)),cwd=PROJECT,stdout=log,stderr=log,check=True)
    with (out/'MEASUREMENTS.tsv').open() as stream:results=summary(list(csv.DictReader(stream,delimiter='\t')))
    (out/'SUMMARY.json').write_text(json.dumps(results,indent=2)+'\n')
    for r in results:
        if r['mode']==1 and r['batch'] in (4096,262144):print(r,flush=True)


if __name__=='__main__':main()
