#!/usr/bin/env python3
"""Bounded paired benchmark and independent controls; never launches discovery."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import random
import shutil
import statistics
import subprocess
import time

from run_smooth_frontier import audit_guided
from odd_order_spin import parse_seed
from run_s18_linked_roots import state_key

PROJECT=Path(__file__).resolve().parents[2]

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def save(path,value): path.write_text(json.dumps(value,indent=2)+'\n')
def run(command,log):
    result=subprocess.run(list(map(str,command)),check=True,capture_output=True,text=True)
    log.write_text(result.stdout+result.stderr)
    return result.stdout

def aggregate(text):
    lines=[line for line in text.splitlines() if line and not line.startswith('#')]
    rows=list(csv.DictReader(lines,delimiter='\t'));groups={}
    for row in rows: groups.setdefault((row['stage'],int(row['degree'])),[]).append(row)
    result=[]
    rng=random.Random(20261012)
    for (stage,n),group in groups.items():
        variants={v:{int(r['repeat']):r for r in group if r['variant']==v} for v in {r['variant'] for r in group}}
        comparisons=[('old','fast')] if 'old' in variants else [('cpu','gpu'),('cpu','gpu_shadow')]
        for baseline,new in comparisons:
            ratios=[]
            for repeat in sorted(variants[baseline]):
                a,b=variants[baseline][repeat],variants[new][repeat]
                if (a['checksum'],a['units'])!=(b['checksum'],b['units']): raise ValueError('benchmark output mismatch')
                ratios.append(float(a['seconds'])/float(b['seconds']))
            draws=sorted(statistics.median(rng.choices(ratios,k=len(ratios))) for _ in range(10000))
            result.append(dict(stage=stage,degree=n,baseline=baseline,new=new,
                               repetitions=len(ratios),paired_median_speedup=statistics.median(ratios),
                               bootstrap_median_interval95=[draws[249],draws[9749]],
                               old_median_seconds=statistics.median(float(r['seconds']) for r in variants[baseline].values()),
                               new_median_seconds=statistics.median(float(r['seconds']) for r in variants[new].values()),
                               units=int(group[0]['units']),paired_ratios=ratios))
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',type=Path,required=True);p.add_argument('--old-binary',type=Path)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--size',type=int,default=100000)
    p.add_argument('--attempts',type=int,default=10000)
    p.add_argument('--end-to-end-only',action='store_true')
    p.add_argument('--skip-metal',action='store_true');args=p.parse_args()
    args.out.mkdir();binary=args.binary.resolve();start=time.time()
    catalog_path=PROJECT/'evidence/smooth_positive_catalog_20261012/CATALOG.json'
    catalog=json.loads(catalog_path.read_text())['components']
    controls=args.out/'controls.tsv'
    controls.write_text(''.join(f"{c['degree']}\t{c['states']}\t{c['representative']}\n" for c in catalog))
    run([binary,'--check-controls',controls],args.out/'CONTROLS.log')
    results=[]
    for n in ([]if args.end_to_end_only else range(11,15)):
        log=args.out/f'CPU_{n}.tsv'
        results+=aggregate(run([binary,'--bench','--d',n,'--size',args.size],log))
        print('CPU_BENCH',n,flush=True)
    if not args.end_to_end_only:results+=aggregate(run([binary,'--io-bench','--out',args.out/'io'],args.out/'IO.tsv'))
    # Finite same-seed streams, not equal-duration unequal-population throughput.
    matched=[]
    for n in [11,12,14]:
        reports=[]
        for variant in ['reference','fast']:
            folder=args.out/f'matched_{n}_{variant}'
            command=[binary,'--d',n,'--seconds',120,'--attempts',10000,'--seed',81,
                     '--legacy-proposals','--out',folder]
            if variant=='reference':command.append('--reference')
            run(command,args.out/f'matched_{n}_{variant}.log')
            report=audit_guided(folder);save(folder/'CHECKS.json',report);(folder/'VERIFIED').write_text('benchmark\n')
            reports.append(report)
        fields=['constructed','scope','transitive','passport','short_reject','duplicate','closed','unresolved']
        if any(reports[0]['construction'][f]!=reports[1]['construction'][f] for f in fields):raise ValueError('finite stream outcomes differ')
        leads=[]
        for variant in ['reference','fast']:
            leads.append((args.out/f'matched_{n}_{variant}/LEADS.tsv').read_text())
        if leads[0]!=leads[1]:raise ValueError('promoted finite streams differ')
        matched.append(dict(degree=n,counts_equal=True,leads_equal=True,reports=reports))
    smokes=[]
    for n,lane in [(11,'forest'),(12,'forest'),(14,'forest'),(12,'blocks')]:
        folder=args.out/f'smoke_{n}_{lane}'
        run([binary,'--d',n,'--lane',lane,'--seconds',2,'--out',folder],args.out/f'smoke_{n}_{lane}.log')
        report=audit_guided(folder);save(folder/'CHECKS.json',report);(folder/'VERIFIED').write_text('benchmark\n');smokes.append(report)
    metal=None
    if not args.end_to_end_only and not args.skip_metal and shutil.which('swiftc'):
        executable=args.out/'forest_screen'
        run(['swiftc','-O',PROJECT/'scripts/gpu/forest_screen.swift','-o',executable],args.out/'METAL_BUILD.log')
        trace=args.out/'GPU_TRACE.bin'
        run([binary,'--bench','--d',14,'--size',262144,'--gpu-trace',trace],args.out/'GPU_TRACE.log')
        metal=aggregate(run([executable,trace,PROJECT/'scripts/gpu/forest_screen.metal'],args.out/'METAL.tsv'))
        results+=metal
    end_to_end=[]
    if args.old_binary:
        for n in [11,12,13,14]:
            timings=[]
            for repeat in range(1,8):
                pair={}
                for variant in (['old','fast'] if repeat%2 else ['fast','old']):
                    folder=args.out/f'end_to_end_{n}_{repeat}_{variant}'
                    executable=args.old_binary.resolve()if variant=='old'else binary
                    command=[executable,'--d',n,'--seconds',120,'--attempts',args.attempts,
                             '--seed',8100+repeat,'--out',folder]
                    if variant=='fast':command.append('--legacy-proposals')
                    run(command,args.out/f'end_to_end_{n}_{repeat}_{variant}.log')
                    report=audit_guided(folder);save(folder/'CHECKS.json',report)
                    (folder/'VERIFIED').write_text('benchmark\n');pair[variant]=report
                fields=['constructed','scope','transitive','passport','short_reject','duplicate','closed','unresolved']
                if any(pair['old']['construction'][f]!=pair['fast']['construction'][f]for f in fields):raise ValueError('old/new end-to-end populations differ')
                if (args.out/f'end_to_end_{n}_{repeat}_old/LEADS.tsv').read_text()!=(args.out/f'end_to_end_{n}_{repeat}_fast/LEADS.tsv').read_text():raise ValueError('end-to-end promoted streams differ')
                for variant in ['old','fast']:
                    timings.append(f"end_to_end\t{n}\t{repeat}\t{variant}\t{pair[variant]['construction']['seconds']}\t{pair[variant]['construction']['constructed']}\t{args.attempts}")
            text='stage\tdegree\trepeat\tvariant\tseconds\tchecksum\tunits\n'+'\n'.join(timings)+'\n'
            (args.out/f'END_TO_END_{n}.tsv').write_text(text)
            end_to_end+=aggregate(text);print('END_TO_END',n,flush=True)
        results+=end_to_end
    save(args.out/'RESULTS.json',dict(seconds=time.time()-start,binary_sha256=sha(binary),
        catalog_sha256=sha(catalog_path),controls=len(catalog),control_states=sum(c['states'] for c in catalog),
        stages=results,matched=matched,smokes=smokes,metal=metal,end_to_end=end_to_end,
        old_binary_sha256=sha(args.old_binary.resolve())if args.old_binary else None,
        limitations='Stage ratios are not discovery speedups; shared busy Mac; bootstrap intervals describe paired timing repeats.'))
    print(json.dumps(results,indent=2),flush=True)

if __name__=='__main__':main()
