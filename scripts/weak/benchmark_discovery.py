#!/usr/bin/env python3
"""Finite paired discovery benchmarks. Technical repeats, not a speed guarantee."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import statistics
import subprocess
import tempfile
import time


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline',type=Path,required=True)
    parser.add_argument('--optimized',type=Path,required=True)
    parser.add_argument('--repeats',type=int,default=5)
    parser.add_argument('--out',type=Path)
    args=parser.parse_args()
    if args.repeats<3:parser.error('at least three alternating technical repeats required')
    out=args.out or Path(tempfile.mkdtemp(prefix='discovery_',dir='runs/benchmarks'))
    out.mkdir(parents=True,exist_ok=True)
    cases=[('d12',[4,6],[1,1]),('d13',[4,6],[3]),('d14',[6,6],[1,1])]
    measurements=[]
    for label,blocks,paths in cases:
        for rep in range(args.repeats):
            variants=[('baseline',args.baseline,[]),('optimized',args.optimized,['--quotient-h'])]
            if rep%2:variants.reverse()
            for variant,binary,extra in variants:
                folder=out/f'{label}_{rep}_{variant}';folder.mkdir()
                cmd=list(map(str,[binary,'hexagon-search','--pair-blocks',','.join(map(str,blocks)),
                    '--pair-paths',','.join(map(str,paths)),'--seconds',30,'--out',folder,*extra]))
                (folder/'COMMAND.json').write_text(json.dumps(cmd))
                start=time.perf_counter()
                with (folder/'RUN.log').open('w') as log:subprocess.run(cmd,stdout=log,stderr=log,check=True)
                elapsed=time.perf_counter()-start
                v=dict(line.split(' = ',1) for line in (folder/'SUMMARY.txt').read_text().splitlines() if ' = ' in line)
                assert v['traversal exhausted']=='true' and v['unresolved candidate starting states']=='0',folder
                assert v['certified weak-failure components']=='0',folder
                measurements.append(dict(case=label,rep=rep,variant=variant,seconds=elapsed,
                    thirds=int(v['compatible third involutions visited']),binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest()))
    with (out/'MEASUREMENTS.tsv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(measurements[0]),delimiter='\t');writer.writeheader();writer.writerows(measurements)
    print(out)
    for label,_,_ in cases:
        medians={variant:statistics.median(r['seconds'] for r in measurements if r['case']==label and r['variant']==variant) for variant in ('baseline','optimized')}
        print(label,medians,'ratio',medians['baseline']/medians['optimized'])


if __name__=='__main__':main()
