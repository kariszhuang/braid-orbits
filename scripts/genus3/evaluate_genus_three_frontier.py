#!/usr/bin/env python3
"""Bounded, blocked producer comparison; never launches the long search pool.

Measure distinct late leads/second, not raw-draw speed across different laws.
Every saved negative/exception is independently audited. Same frozen word bank.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import csv
import json
from pathlib import Path
import random
import subprocess
import time

from audit_guided_native import audit_guided_native
from run_genus_three import write, sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['binary','auditor','words','out']:
        p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--seconds',type=int,default=3)
    p.add_argument('--degrees',type=int,nargs='+',default=[15,18,21,24])
    args=p.parse_args()
    if args.seconds<=0 or any(not 9<=n<=25 for n in args.degrees):
        p.error('positive seconds and degrees 9..25 are required')
    args.out.mkdir();order=random.Random(20261022);rows=[]
    for n in args.degrees:
        for rep in range(3):
            lanes=['collision','joint-genus3'];order.shuffle(lanes)
            block={}
            for lane in lanes:
                folder=args.out/f'{n}_{rep}_{lane}';cmd=[str(args.binary),'--d',str(n),'--lane',lane,
                    '--seed',str(98651+rep),'--seconds',str(args.seconds),'--cap','1000000',
                    '--closure-seconds','10','--words',str(args.words),'--out',str(folder)]
                began=time.perf_counter();subprocess.run(cmd,capture_output=True,text=True,check=True)
                wall=time.perf_counter()-began
                checked=audit_guided_native(folder,args.auditor,native_components=True);write(folder/'CHECKS.json',checked)
                with (folder/'COUNTS.tsv').open() as f:c=next(csv.DictReader(f,delimiter='\t'))
                with (folder/'STAGES.tsv').open() as f:s=next(csv.DictReader(f,delimiter='\t'))
                # Without an unbounded orbit walk, distinct labelled late leads
                # are only an upper bound on distinct S_n states/components.
                leads={line.split('\t',1)[1] for line in (folder/'LEADS.tsv').read_text().splitlines()}
                block[lane]=dict(wall_seconds=wall,counts=c,stages=s,distinct_labelled_leads=len(leads),
                    labelled_leads_per_second=len(leads)/wall,independent=checked['independent_native'])
            rows.append(dict(degree=n,replicate=rep,results=block))
            print(n,rep,{k:round(v['labelled_leads_per_second'],3) for k,v in block.items()},flush=True)
    write(args.out/'RESULT.json',dict(passed=True,binary_sha256=sha(args.binary),auditor_sha256=sha(args.auditor),
        words_sha256=sha(args.words),blocks=rows,scope='different biased populations,shared word bank; labelled distinct leads are not distinct components; timings under continuous main-search load'))


if __name__=='__main__':main()
