"""Matched fixed-proposal third-pass timings, with independent omitted-seed checks."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import csv
import json
from pathlib import Path
import random
import statistics
import subprocess
import time

from audit_guided_native import audit_guided_native
from benchmark_forest_second import interval
from run_smooth_degree9 import sha

PROJECT=Path(__file__).resolve().parents[2]
CASES=[(11,'forest',256),(12,'forest',256),(13,'forest',256),(14,'forest',256),(12,'blocks',256),
       (11,'forest',1),(14,'forest',1024)]


def table(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))


def omitted_negatives(old,new,binary,auditor,n,attempts):
    ids={row['candidate'] for row in table(new/'RESOLUTIONS.tsv')}
    omitted=[row for row in table(old/'RESOLUTIONS.tsv') if row['candidate'] not in ids]
    if any(row['status']!='nonsmooth' for row in omitted):raise ValueError('new sieve omitted a nonnegative old outcome')
    root=new/'omitted_negatives';root.mkdir()
    seeds=root/'INPUT.txt';seeds.write_text(''.join(row['seed']+'\n' for row in omitted))
    words=subprocess.check_output([str(binary),'--d',str(n),'--screen-seeds',str(seeds)],text=True).splitlines()
    if len(words)!=len(omitted) or 'survives' in words:raise ValueError('unexplained candidate omission')
    (root/'LEADS.tsv').write_text(''.join(row['candidate']+'\t'+row['seed']+'\n' for row in omitted))
    (root/'RESOLUTIONS.tsv').write_text('candidate\tstatus\tword\tfolder\tseed\n'+''.join(
        row['candidate']+'\tnonsmooth\t'+word+'\t\t'+row['seed']+'\n' for row,word in zip(omitted,words)))
    header='degree\tlane\tseed\tconstructed\tscope\ttransitive\tpassport\tshort_reject\tduplicate\tclosed\tunresolved\tseconds\n'
    k=len(omitted)
    (root/'COUNTS.tsv').write_text(header+f'{n}\tforest\t0\t{attempts}\t{k}\t{k}\t{k}\t0\t0\t0\t0\t0\n')
    result=audit_guided_native(root,auditor)
    (root/'CHECKS.json').write_text(json.dumps(result,indent=2)+'\n')
    return k


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--old',type=Path,required=True);p.add_argument('--new',type=Path,required=True)
    p.add_argument('--auditor',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--attempts',type=int,default=500000);p.add_argument('--repeats',type=int,default=7)
    p.add_argument('--pilot',action='store_true');p.add_argument('--deep-level',type=int,default=3)
    args=p.parse_args();root=args.out.resolve();root.mkdir();old=args.old.resolve();new=args.new.resolve();auditor=args.auditor.resolve()
    catalog=json.loads((PROJECT/'evidence/smooth_positive_catalog_20261012/CATALOG.json').read_text())['components']
    controls=root/'controls.tsv';controls.write_text(''.join(f"{c['degree']}\t{c['states']}\t{c['representative']}\n" for c in catalog))
    text=subprocess.check_output([str(new),'--check-controls',str(controls)],text=True)
    (root/'CONTROLS.log').write_text(text);print('CONTROLS',text.strip(),flush=True)
    rng=random.Random(202610123);results=[];omitted=0
    cases=[(11,'forest',256),(14,'forest',256),(12,'blocks',256),(11,'forest',1)] if args.pilot else CASES
    levels=list(range(4)) if args.pilot else [args.deep_level]
    for n,lane,reuse in cases:
        for level in levels:
            timings=[]
            for repeat in range(args.repeats):
                order=['old','new'];rng.shuffle(order);pair={};folders={};reports={}
                seed=8400+repeat if args.pilot else 9500+repeat
                for variant in order:
                    folder=root/f'{lane}_{n}_reuse{reuse}_level{level}_{repeat}_{variant}';folders[variant]=folder
                    cmd=[str(old if variant=='old' else new),'--d',str(n),'--lane',lane,'--frame-reuse',str(reuse),
                         '--seconds','120','--attempts',str(args.attempts),'--seed',str(seed),'--out',str(folder)]
                    if variant=='new':
                        cmd+=['--deep-level',str(level if lane=='forest' or args.pilot else 0)]
                        if lane=='blocks' and not args.pilot:cmd+=['--no-orbit-cache']
                    began=time.monotonic();subprocess.run(cmd,check=True,stdout=subprocess.DEVNULL)
                    producer=time.monotonic()-began;began=time.monotonic()
                    reports[variant]=audit_guided_native(folder,auditor);audit=time.monotonic()-began
                    (folder/'CHECKS.json').write_text(json.dumps(reports[variant],indent=2)+'\n')
                    pair[variant]=dict(producer=producer,audit=audit,total=producer+audit,order=order.index(variant))
                for field in ['constructed','scope','passport','closed','unresolved']:
                    if reports['old']['construction'][field]!=reports['new']['construction'][field]:raise ValueError('scope/proposal/positive result changed: '+field)
                if int(reports['new']['construction']['constructed'])!=args.attempts:raise ValueError('timing cap hit')
                if sha(folders['old']/'FRAMES.tsv')!=sha(folders['new']/'FRAMES.tsv'):raise ValueError('RNG/frame stream changed')
                omitted+=omitted_negatives(folders['old'],folders['new'],new,auditor,n,args.attempts)
                pair['new']['deep_reject']=int(table(folders['new']/'POLICY.tsv')[0]['deep_reject'])
                timings.append(pair)
            result=dict(degree=n,lane=lane,reuse=reuse,
                        deep_level=level if lane=='forest' or args.pilot else 0,timings=timings,
                        producer=interval([r['old']['producer']/r['new']['producer'] for r in timings]),
                        pipeline=interval([r['old']['total']/r['new']['total'] for r in timings]))
            results.append(result);print('PIPELINE',n,lane,reuse,level,json.dumps(result['pipeline']),flush=True)
    result=dict(passed=True,pilot=args.pilot,controls=len(catalog),control_states=sum(c['states'] for c in catalog),
                binary_sha256=sha(new),audit_binary_sha256=sha(auditor),old_binary_sha256=sha(old),
                repeats=args.repeats,attempts=args.attempts,stages=results,independently_checked_omissions=omitted,
                identical_proposal_and_frame_streams=True,
                limitations='Fixed biased proposal streams, independent seed replicates; no discovery-rate promise. Deep filters change saved-lead/counter distribution. Startup/T7 logging/audit included; omitted-seed audit is separate validation, not production runtime.')
    (root/'RESULTS.json').write_text(json.dumps(result,indent=2)+'\n');print('RESULTS='+str(root/'RESULTS.json'),flush=True)


if __name__=='__main__':main()
