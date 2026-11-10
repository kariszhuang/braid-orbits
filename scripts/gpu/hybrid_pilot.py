#!/usr/bin/env python3
"""Frozen native hybrid runtime and paired experiments for run_campaign.py."""
import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path
import random
import resource
import shutil
import signal
import struct
import subprocess
import sys
import time

PROJECT=Path(__file__).resolve().parents[2]
CHILD=None
STOP=False
CASES=[(d,f) for d in (9,10,11) for f in (0,1)]+[(18,2),(18,3),(25,4)]
STAT_NAMES=['screen_survivors','point_transitive','duplicate_starts','known_bad_samples','proper_groups','new_generating_starts','probe_witnesses','exact_witnesses','new_verified_failures','unresolved_starts','explored_states','cached_good_samples']

def walk_oracle():
    """Independent product-one S4 oracle; all four exact states after 64 moves."""
    mask=(1<<64)-1;perms=list(itertools.permutations(range(4)));ranks={p:i for i,p in enumerate(perms)}
    raw=bytearray();packed=[]
    def compose(a,b):return tuple(a[b[x]] for x in range(4))
    def inverse(a):return tuple(a.index(x) for x in range(4))
    for a,b,c in itertools.product(perms,repeat=3):
        raw.extend(a+b+c);t=[a,b,c,inverse(compose(a,compose(b,c)))];state=0xcbf29ce484222325
        for x in a+b+c:state=((state^x)*0x100000001b3)&mask
        for _ in range(64):
            while True:
                state=(state+0x9e3779b97f4a7c15)&mask;z=state
                z=((z^(z>>30))*0xbf58476d1ce4e5b9)&mask;z=((z^(z>>27))*0x94d049bb133111eb)&mask;z^=z>>31
                if z<mask-mask%6:break
            move=z%6;i=move%3;x,y=t[i:i+2]
            t[i:i+2]=[compose(x,compose(y,inverse(x))),x] if move<3 else [y,compose(inverse(y),compose(x,y))]
        packed.append(sum(ranks[p]<<(5*j) for j,p in enumerate(t)))
    return bytes(raw),struct.pack(f'<{len(packed)}I',*packed)

def atomic(path,value):
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value,indent=2)+'\n');temporary.replace(path)

def stop(signum,frame):
    global STOP
    STOP=True
    if CHILD is not None and CHILD.poll() is None:CHILD.terminate()

def schedule(workers,batches,repeats,smoke=False):
    rng=random.Random(20260929)
    blocks=list(itertools.product(CASES,workers,batches,range(repeats)))
    if smoke:blocks=[((d,f),6,4096,0) for d,f in [(9,0),(10,1),(11,1),(18,2),(18,3),(25,4)]]
    rng.shuffle(blocks);rows=[]
    for block,((d,f),w,n,rep) in enumerate(blocks):
        variants=['cpu','gpu'];rng.shuffle(variants)
        for backend in variants:
            rows.append(dict(name=f'b{block:03d}_{backend}',backend=backend,degree=d,family=f,workers=w,batch=n,repeat=rep,
                             seed=20260929+rep*104729+d*8191+f*3571,block=block))
    return rows

def summarize(out,plan):
    rows=[];components={};paired={}
    for spec in plan['schedule']:
        done=out/'jobs'/spec['name']/'DONE.json'
        if not done.exists():continue
        entry=json.loads(done.read_text());attempt=out/entry['attempt'];summary=json.loads((attempt/'SUMMARY.json').read_text())
        stats=[sum(r['stats'][i] for r in summary['rows']) for i in range(12)]
        row=dict(spec,pipeline_seconds=summary['pipeline_seconds'],audit_seconds=summary['shadow_audit_seconds'],
                 cpu_seconds=entry['cpu_seconds'],input_sha256=entry['input_sha256'],**dict(zip(STAT_NAMES,stats)))
        rows.append(row);paired.setdefault(spec['block'],{})[spec['backend']]=row
        for path in attempt.glob('CERTIFIED_*.txt'):
            h=dict(line.split('=',1) for line in path.read_text().splitlines() if '=' in line and not line.startswith('STATE'))
            identity=(h['d'],h['representative'])
            components.setdefault(identity,dict(path=str(path.relative_to(out)),degree=int(h['d']),states=int(h['states']),group=h['group']))
    pairs=[]
    for block,arms in paired.items():
        if set(arms)!= {'cpu','gpu'}:continue
        c,g=arms['cpu'],arms['gpu']
        if c['input_sha256']!=g['input_sha256']:raise RuntimeError(f'Paired input mismatch: block {block}')
        pairs.append(dict(block=block,degree=c['degree'],family=c['family'],workers=c['workers'],batch=c['batch'],repeat=c['repeat'],
                          cpu_seconds=c['pipeline_seconds'],gpu_seconds=g['pipeline_seconds'],pipeline_ratio=c['pipeline_seconds']/g['pipeline_seconds'],
                          cpu_new_generating=c['new_generating_starts'],gpu_new_generating=g['new_generating_starts'],
                          cpu_unresolved=c['unresolved_starts'],gpu_unresolved=g['unresolved_starts']))
    for name,records in [('METRICS.tsv',rows),('PAIRED.tsv',pairs)]:
        if records:
            with (out/name).open('w') as stream:
                writer=csv.DictWriter(stream,fieldnames=list(records[0]),delimiter='\t');writer.writeheader();writer.writerows(records)
    atomic(out/'RESULTS.json',dict(completed_arms=len(rows),planned_arms=len(plan['schedule']),paired_blocks=len(pairs),
           verified_new_components=list(components.values()),caveat='Sample/start counts are not distinct components. Technical repetitions are not independent discoveries. Time-limited closure is unresolved, not a counterexample.'))

def main():
    global CHILD
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--max-workers',type=int,default=10)
    p.add_argument('--prepare-only',action='store_true',help='freeze/build the shared discovery runtime, without launching candidate work')
    p.add_argument('--known',type=Path,default=PROJECT/'evidence/certificates',help='immutable starting certificate catalog')
    p.add_argument('--smoke',action='store_true');p.add_argument('--repeats',type=int,default=2)
    p.add_argument('--compact-inputs',action='store_true',help='frozen generation recipes plus SHA256; packed exact flags; full regeneration/audit')
    p.add_argument('--early-witness',action='store_true',help='GPU/CPU 64-move witness rejection BEFORE full-group checks')
    p.add_argument('--batch-sizes',default='65536,262144')
    p.add_argument('--batches',type=int,default=3);p.add_argument('--arm-seconds',type=int,default=600)
    p.add_argument('--min-free-gb',type=float,default=30);p.add_argument('--max-seconds',type=float,default=0)
    a=p.parse_args();out=a.out.resolve();out.mkdir(parents=True,exist_ok=True)
    try:batch_sizes=list(dict.fromkeys(map(int,a.batch_sizes.split(','))))
    except ValueError:p.error('batch-sizes must be comma-separated integers')
    if not 6<=a.max_workers<=10 or not 1<=a.repeats<=7 or not 1<=a.batches<=128 or a.arm_seconds<1 or not batch_sizes or any(not 1<=n<=1048576 for n in batch_sizes):p.error('invalid finite limits')
    signal.signal(signal.SIGTERM,stop);signal.signal(signal.SIGINT,stop)
    plan_path=out/'PLAN.json'
    if not plan_path.exists():
        snapshot=out/'source';snapshot.mkdir()
        paths=sorted(str(p.relative_to(PROJECT)) for p in (PROJECT/'src').rglob('*.rs'))
        paths+=['scripts/gpu/gpu_screen_cpu.rs','scripts/gpu/gpu_screen.metal','scripts/gpu/hybrid_pilot.swift','scripts/gpu/hybrid_pilot.py']
        hashes={}
        for name in paths:
            path=PROJECT/name;target=snapshot/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(path,target);hashes[name]=hashlib.sha256(path.read_bytes()).hexdigest()
        known=snapshot/'evidence/certificates';known.mkdir(parents=True)
        for path in a.known.glob('CERTIFIED_*weak.txt'):
            shutil.copy2(path,known/path.name);hashes[f'evidence/certificates/{path.name}']=hashlib.sha256(path.read_bytes()).hexdigest()
        workers=[w for w in (6,8,10) if w<=a.max_workers]
        plan=dict(schedule=[] if a.prepare_only else schedule(workers,batch_sizes,a.repeats,a.smoke),batches=a.batches,
                  cap=50000,millis=250,source_hashes=hashes,arm_seconds=a.arm_seconds,
                  screen_mode=2 if a.early_witness else 1,runtime_only=a.prepare_only,compact_inputs=a.compact_inputs,
                  comparison='Matched candidates and frozen known-component cache; randomized block and backend order; complete pipeline includes batch persistence, CPU full-group checks and capped exact closure; shadow audit separately timed.')
        atomic(plan_path,plan)
    plan=json.loads(plan_path.read_text());snapshot=out/'source'
    for name,digest in plan['source_hashes'].items():
        if hashlib.sha256((snapshot/name).read_bytes()).hexdigest()!=digest:raise RuntimeError(f'Snapshot modified: {name}')
    if plan.get('screen_mode',1)==2 and not (out/'ORACLE_WALK_OUTPUT.bin').exists():
        raw,expected=walk_oracle();(out/'ORACLE_WALK_INPUT.bin').write_bytes(raw);(out/'ORACLE_WALK_OUTPUT.bin').write_bytes(expected)
    for name,cmd in [('libscreen.dylib',['rustc','--crate-type','cdylib','--edition','2021','-C','opt-level=3','-C','target-cpu=native','-C','codegen-units=1',snapshot/'scripts/gpu/gpu_screen_cpu.rs']),
                     ('hybrid_pilot',['xcrun','swiftc','-O','-whole-module-optimization',snapshot/'scripts/gpu/hybrid_pilot.swift'])]:
        if not (out/name).exists():
            with (out/'BUILD.log').open('a') as log:subprocess.run(list(map(str,cmd+['-o',out/name])),stdout=log,stderr=log,check=True)
    if a.prepare_only:
        if not plan.get('runtime_only') or plan.get('screen_mode')!=2:raise RuntimeError('Discovery requires a prepared early-witness runtime')
        built={name:hashlib.sha256((out/name).read_bytes()).hexdigest() for name in ('libscreen.dylib','hybrid_pilot')}
        metadata=out/'NATIVE_HASHES.json'
        if metadata.exists() and json.loads(metadata.read_text())!=built:raise RuntimeError('Prepared native runtime modified')
        atomic(metadata,built)
        print(f'PREPARED {out}',flush=True);return
    deadline=time.monotonic()+a.max_seconds if a.max_seconds else None
    duration_limited=False
    for spec in plan['schedule']:
        if STOP:break
        if deadline is not None and time.monotonic()>=deadline:duration_limited=True;break
        job=out/'jobs'/spec['name']
        if (job/'DONE.json').exists():continue
        if shutil.disk_usage(out).free<a.min_free_gb*1024**3:raise RuntimeError('Disk reserve reached; all work retained')
        job.mkdir(parents=True,exist_ok=True);attempt=job/f'attempt_{len(list(job.glob("attempt_*"))):03d}';attempt.mkdir()
        command=[out/'hybrid_pilot',out/'libscreen.dylib',snapshot/'scripts/gpu/gpu_screen.metal',attempt,snapshot/'evidence/certificates',spec['backend'],spec['degree'],spec['family'],spec['batch'],plan['batches'],spec['workers'],spec['seed'],plan['cap'],plan['millis']]
        if plan.get('screen_mode',1)==2:command.append(2)
        if plan.get('compact_inputs'):
            if plan.get('screen_mode',1)!=2:command.append(1)
            command.append(1)
        atomic(attempt/'COMMAND.json',list(map(str,command)));atomic(out/'STATUS.json',dict(status='running',job=spec,attempt=str(attempt),time=time.time()))
        print(f'RUN {spec}',flush=True);before=resource.getrusage(resource.RUSAGE_CHILDREN)
        with (attempt/'RUN.log').open('w') as log:
            CHILD=subprocess.Popen(list(map(str,command)),stdout=log,stderr=log)
            try:rc=CHILD.wait(timeout=min(plan['arm_seconds'],max(0.01,deadline-time.monotonic())) if deadline is not None else plan['arm_seconds'])
            except subprocess.TimeoutExpired:
                CHILD.terminate()
                try:CHILD.wait(timeout=10)
                except subprocess.TimeoutExpired:CHILD.kill();CHILD.wait()
                atomic(attempt/'CAPPED.json',dict(reason='arm wall-time ceiling; saved inputs and validation ledger retained'))
                rc=-1
        CHILD=None
        if STOP:break
        if rc==0:
            after=resource.getrusage(resource.RUSAGE_CHILDREN)
            digest=hashlib.sha256()
            if plan.get('compact_inputs'):
                for path in sorted(attempt.glob('batch_*.input.json')):digest.update(bytes.fromhex(json.loads(path.read_text())['sha256']))
            else:
                for path in sorted(attempt.glob('batch_*.input.bin')):digest.update(path.read_bytes())
            atomic(job/'DONE.json',dict(attempt=str(attempt.relative_to(out)),input_sha256=digest.hexdigest(),
                   digest_scheme='sha256_of_batch_sha256_bytes' if plan.get('compact_inputs') else 'sha256_of_raw_batch_bytes',
                   cpu_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime))
        elif rc!=-1:raise RuntimeError(f'Failed arm retained: {attempt}; see RUN.log')
        summarize(out,plan)
    summarize(out,plan);atomic(out/'STATUS.json',dict(status='stopped' if STOP else 'duration-limited' if duration_limited else 'complete',time=time.time()))
    if not STOP and not duration_limited and any(not (out/'jobs'/r['name']/'DONE.json').exists() for r in plan['schedule']):
        raise RuntimeError('Some arms capped; inputs retained, no complete comparison claimed')

if __name__=='__main__':main()
