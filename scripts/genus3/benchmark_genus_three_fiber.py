#!/usr/bin/env python3
"""Randomized paired performance blocks of raw returns versus conjugacy solving.

The two proposal laws differ. Report qualified-candidate throughput, never a
same-stream speedup or multiplier in the probability of discovery.
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
import statistics
import subprocess
import time

from audit_guided_native import audit_guided_native
from run_genus_three import sha, write


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['raw','solved','auditor','out']:
        p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--attempts',type=int,default=5_000_000)
    p.add_argument('--raw-rejection-flag',action='store_true',help='use the diagnostic raw path in the same final binary')
    args=p.parse_args();args.out.mkdir()
    order=random.Random(20261022);rows=[]
    for degree in [18,24]:
        for replicate in range(7):
            modes=['raw','solved'];order.shuffle(modes);block={}
            for mode in modes:
                dest=args.out/f'{degree}_{replicate}_{mode}'
                binary=getattr(args,mode)
                cmd=[str(binary),'--d',str(degree),'--lane','full-fiber','--seed',str(931+replicate),
                    '--attempts',str(args.attempts),'--seconds','120','--cap','1000000',
                    '--closure-seconds','20','--out',str(dest)]
                if mode=='raw' and args.raw_rejection_flag:cmd.append('--fiber-rejection')
                began=time.perf_counter();subprocess.run(cmd,check=True,capture_output=True,text=True)
                wall=time.perf_counter()-began
                with (dest/'COUNTS.tsv').open() as stream:counts=next(csv.DictReader(stream,delimiter='\t'))
                assert int(counts['constructed'])==args.attempts
                checks=audit_guided_native(dest,args.auditor,native_components=True)
                write(dest/'CHECKS.json',checks)
                block[mode]=dict(wall_seconds=wall,counts=counts,
                    qualified_per_second=int(counts['scope'])/wall,
                    independent=checks['independent_native'])
            rows.append(dict(degree=degree,replicate=replicate,results=block))
            print(f'checked degree {degree} block {replicate+1}/7',flush=True)
    summaries={}
    for degree in [18,24]:
        ratios=[r['results']['solved']['qualified_per_second']/r['results']['raw']['qualified_per_second']
                for r in rows if r['degree']==degree]
        bootstrap=random.Random(199);means=sorted(statistics.mean(bootstrap.choices(ratios,k=7)) for _ in range(10000))
        summaries[str(degree)]=dict(mean_paired_throughput_ratio=statistics.mean(ratios),
            bootstrap95=[means[249],means[9749]],
            raw_median_qualified_per_second=statistics.median(r['results']['raw']['qualified_per_second'] for r in rows if r['degree']==degree),
            solved_median_qualified_per_second=statistics.median(r['results']['solved']['qualified_per_second'] for r in rows if r['degree']==degree))
    report=dict(passed=True,raw_sha256=sha(args.raw),solved_sha256=sha(args.solved),
        auditor_sha256=sha(args.auditor),blocks=rows,summary=summaries,
        scope='same marked quotient templates but DIFFERENT proposal laws; qualified-genus3 throughput, not a discovery multiplier')
    write(args.out/'RESULT.json',report);print(json.dumps(summaries,indent=2))


if __name__=='__main__':main()
