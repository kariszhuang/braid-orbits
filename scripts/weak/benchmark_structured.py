#!/usr/bin/env python3
"""Sequential before/after complete-passport measurements, including class setup."""
import argparse
import hashlib
import json
from pathlib import Path
import statistics
import subprocess
import time

p=argparse.ArgumentParser(__doc__)
p.add_argument('--baseline',type=Path,required=True)
p.add_argument('--current',type=Path,required=True)
p.add_argument('--out',type=Path,required=True)
args=p.parse_args()
args.out.mkdir(parents=True,exist_ok=False)
rows=[]
for d,passport in [(9,'2.2.2.2/2.2.2.2/2.2.2.2/3.3.3'),
                   (10,'3.3.3/2.2.2.2/2.2.2.2/3.3.3'),
                   (11,'3.3.3/2.2.2.2/2.2.2.2/3.3.3')]:
    for repetition in range(3):
        # Alternate order to reduce warm-cache and thermal-order bias.
        for name in (('baseline','current') if repetition%2==0 else ('current','baseline')):
            binary=getattr(args,name).resolve()
            folder=args.out/f'{d}_{repetition}_{name}'
            cmd=[str(binary),'passport-exact','--group','A','--d',str(d),'--passport',passport,
                 '--weak-only','--threads','5','--seconds','60','--seed','20260929','--out',str(folder)]
            if name=='current':cmd+=['--centralizer-rows']
            started=time.perf_counter()
            result=subprocess.run(cmd,capture_output=True,text=True)
            elapsed=time.perf_counter()-started
            (args.out/f'{folder.name}.log').write_text(result.stdout+result.stderr)
            result.check_returncode()
            v=dict(l.split(' = ',1) for l in (folder/'SUMMARY.txt').read_text().splitlines() if ' = ' in l)
            if v['scan range resolved']!='true' or v['scan next']!=v['scan total candidate pairs']:
                raise RuntimeError(f'Incomplete benchmark: {folder}; no speedup claim permitted')
            for cert in folder.glob('CERTIFIED_*.txt'):
                subprocess.run([str(binary),'verify','--certificate',str(cert)],check=True)
            rows.append(dict(degree=d,passport=passport,implementation=name,repetition=repetition,
                             wall_seconds=elapsed,binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                             command=cmd,summary=str(folder/'SUMMARY.txt'),candidate_pairs=int(v['scan total candidate pairs'])))
            (args.out/'MEASUREMENTS.json').write_text(json.dumps(rows,indent=2)+'\n')
            print(d,name,f'{elapsed:.4f}s',flush=True)
for d in (9,10,11):
    times={name:statistics.median(r['wall_seconds'] for r in rows if r['degree']==d and r['implementation']==name)
           for name in ('baseline','current')}
    print(d,'median speedup',times['baseline']/times['current'])
