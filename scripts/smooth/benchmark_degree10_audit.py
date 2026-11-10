#!/usr/bin/env python3
"""Matched prefix audit baseline; full-frame validation is separate."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse, itertools, json, statistics, subprocess, time
from pathlib import Path
from run_smooth_degree9 import binary_records, rows
from run_smooth_mission import values
from odd_order_spin import parse_seed, inverse, cycle_lengths, hurwitz
from analyze_fixed_braid_constraints import product, point_orbits

def reference(folder, a, b, limit):
    n=len(a); identity=tuple(range(n)); p=product(a,b); ip=n-len(cycle_lengths(p))
    left=len(point_orbits((a,b))); positive=(folder.parents[1]/'POSITIVE_ONLY').exists()
    dictionary=rows(folder/'WITNESSES.tsv')
    words=[]
    for row in dictionary:
        if row['status']!='nonsmooth': raise ValueError('negative-only benchmark scope')
        words.append(tuple(map(int,row['word_or_component'].split(','))))
    records=iter(binary_records(folder/'OUTCOMES.bin')); count=[0]*8; omitted=[0,0]
    for rank,c in enumerate(itertools.islice(itertools.permutations(range(n)),limit)):
        ic=n-len(cycle_lengths(c)); count[0]+=1; d=inverse(product(p,c)); t=(a,b,c,d)
        if positive:
            if any(q==identity for q in t): omitted[0]+=1; continue
            if ip+ic+n-len(cycle_lengths(d))<=2*n-2: omitted[1]+=1; continue
        if len(point_orbits(t[:3]))!=1: continue
        count[1]+=1; right=len(point_orbits(t[2:]))
        if n-ip-left-right+1!=0 or right!=1: raise ValueError('implication')
        if ic+n-len(cycle_lengths(d))==ip: count[2]+=1
        if any(q==identity for q in t): count[4]+=1; continue
        count[3]+=1; record_rank,token=next(records)
        if record_rank!=rank or token>=len(words): raise ValueError('rank/token')
        for move in words[token]: t=hurwitz(t,move)
        if values(t)['smooth']: raise ValueError('false word')
        count[5]+=1
    return dict(count=count,omitted=omitted)

def main():
    p=argparse.ArgumentParser();p.add_argument('folder',type=Path);p.add_argument('--binary',type=Path)
    p.add_argument('--limit',type=int,default=10000);p.add_argument('--repeat',type=int,default=3)
    p.add_argument('--frames',default='0,1000,2000,4000');p.add_argument('--out',type=Path,required=True)
    args=p.parse_args();frames=rows(args.folder/'FRAMES.tsv');results=[]
    for fi in map(int,args.frames.split(',')):
        row=frames[fi];a,b,_,_=parse_seed(f"[{row['a']}, {row['b']}, (), ()]",int(rows(args.folder/'COVERAGE.tsv')[0]['degree']))
        folder=args.folder/'frames'/f'frame_{fi:05}';times={'python':[],'native':[]}
        command='\t'.join([str(folder),','.join(map(str,a)),','.join(map(str,b)),str(int((args.folder/'POSITIVE_ONLY').exists())),str(args.limit)])+'\n'
        native=None
        if args.binary:
            native=subprocess.Popen([str(args.binary),'--server',str(len(a))],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True)
            if native.stdout.readline().strip()!='READY':raise ValueError('native startup')
        try:
            for repeat in range(args.repeat):
                answers={}
                for kind in (['python','native'] if repeat%2==0 else ['native','python']):
                    if kind=='native' and native is None:continue
                    start=time.perf_counter()
                    if kind=='python':answers[kind]=reference(folder,a,b,args.limit)
                    else:
                        native.stdin.write(command);native.stdin.flush();answers[kind]=json.loads(native.stdout.readline())
                    times[kind].append(time.perf_counter()-start)
                if native and answers['python']!={k:answers['native'][k] for k in ['count','omitted']}:raise ValueError('differential mismatch')
            results.append(dict(frame=fi,limit=args.limit,timings=times,counts=answers['python'],
                speedup=statistics.median(times['python'])/statistics.median(times['native']) if native else None))
        finally:
            if native:native.stdin.close();native.wait(timeout=10)
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps(results,indent=2))

if __name__=='__main__':main()
