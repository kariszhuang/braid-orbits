#!/usr/bin/env python3
"""Bounded policy-eleven validation; writes new evidence, never launches a hunt."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from run_campaign import stability_passport_plan, inherited_passport_exclusions


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--predecessor',type=Path,required=True)
    parser.add_argument('--binary',type=Path,default=ROOT/'target/smart/release/hurwitz_breaker')
    args=parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=False)
    results=[]
    def run(name,command):
        start=time.monotonic()
        with (args.out/f'{name}.log').open('w') as log:
            process=subprocess.run(list(map(str,command)),cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
        record=dict(name=name,command=list(map(str,command)),returncode=process.returncode,seconds=time.monotonic()-start)
        results.append(record)
        (args.out/'CHECKS.json').write_text(json.dumps(dict(checks=results),indent=2)+'\n')
        if process.returncode:raise RuntimeError(f'Check failed: {name}; see {args.out}')
        print(f'{name}: passed ({record["seconds"]:.2f}s)',flush=True)
    run('make_check',['make','check'])
    run('make_check_spectra',['make','check-spectra'])
    certs=sorted((ROOT/'evidence/certificates').glob('CERTIFIED_*_weak.txt'))
    if len(certs)!=9:raise RuntimeError('Expected nine weak controls')
    for i,cert in enumerate(certs):run(f'control_{i:02}',[args.binary,'verify','--certificate',cert])
    for group,d in [('S',10),('A',11),('A',14),('A',25)]:
        for rounds in (0,8):
            name=f'joint_{group}{d}_r{rounds}'
            run(name,[args.binary,'smart-hunt','--group',group,'--d',d,'--family','joint-pure',
                '--weak-only','--all-pairs-intransitive','--threads',1,'--seconds',2,
                '--candidate-seconds',1,'--cap',5000,'--steps',16,'--screen-steps',16,
                '--seed',20261010,'--local-block-max',min(d,11),'--local-frame-cap',128,
                '--local-combination-cap',10000,'--point-orbit-max',d,
                '--joint-repair-rounds',rounds,'--out',args.out/name])
            summary=(args.out/name/'SUMMARY.txt').read_text()
            if 'joint pure constraint =' not in summary:raise RuntimeError('New constructor was not exercised')
    for cycle_type in ('5.4.1','6.3.1'):
        name=f'fiber_{cycle_type}'
        run(name,[args.binary,'braid-library','--fiber-actions','--d',10,'--fiber-type',cycle_type,
            '--seconds',60,'--out',args.out/name])
        summary=(args.out/name/'SUMMARY.txt').read_text()
        if 'search complete = true' not in summary:raise RuntimeError('Finite fiber pilot did not finish within budget')
    for group,d,name in [('A',25,'CERTIFIED_A25_XXYY_genus4_weak.txt'),('S',18,'CERTIFIED_S18_genus1_weak.txt')]:
        label=f'joint_source_{group}{d}'
        run(label,[args.binary,'smart-hunt','--group',group,'--d',d,'--family','joint-pure',
            '--joint-source-only','--certificate',ROOT/'evidence/certificates'/name,
            '--weak-only','--all-pairs-intransitive','--threads',1,'--seconds',2,
            '--candidate-seconds',1,'--cap',5000,'--steps',16,'--joint-repair-rounds',8,
            '--seed',20261010,'--skip-known-starts','--known-certificates',ROOT/'evidence/certificates',
            '--out',args.out/label])
    exclusions=inherited_passport_exclusions(args.predecessor)
    plan=stability_passport_plan(exclusions)
    (args.out/'PASSPORT_PLAN.json').write_text(json.dumps(plan,indent=2)+'\n')
    for group,d in [('S',9),('A',10),('A',11),('S',11)]:
        row=next(r for r in plan['rows'] if r['group']==group and r['d']==d)
        name=f'exact_{group}{d}'
        run(name,[args.binary,'passport-exact','--group',group,'--d',d,'--passport',row['passport'],
            '--weak-only','--all-pairs-intransitive','--centralizer-rows','--short-weak-witnesses',
            '--threads',1,'--seconds',2,'--candidate-seconds',1,'--cap',5000,'--steps',16,
            '--count',1024,'--out',args.out/name])
    sources=[ROOT/'scripts/weak/run_campaign.py',ROOT/'src/smart.rs',ROOT/'src/smart/braid_frames.rs',
        ROOT/'src/smart/weak_discovery.rs',ROOT/'src/smart/weak_discovery/joint_pure.rs',Path(__file__)]
    (args.out/'CHECKS.json').write_text(json.dumps(dict(checks=results,
        binary=str(args.binary),binary_sha256=hashlib.sha256(args.binary.read_bytes()).hexdigest(),
        sources={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        workers=os.cpu_count(),predecessor=str(args.predecessor),exact_jobs=len(plan['rows']),
        interpretation='bounded correctness/operational checks, not a discovery-rate benchmark'),indent=2)+'\n')
    print(f'All {len(results)} checks passed. Exact queue: {len(plan["rows"])} passports.',flush=True)


if __name__=='__main__':main()
