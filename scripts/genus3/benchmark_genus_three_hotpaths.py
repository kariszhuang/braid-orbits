#!/usr/bin/env python3
"""Fixed-stream v1/v2 CPU comparison; static isotrivial exclusions measured apart."""
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
from run_genus_three import sha,write
from verify_genus_three_exclusions import audit_exclusions


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['baseline','optimized','auditor','mutations','source','out']:
        p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--attempts',type=int,default=500000)
    p.add_argument('--words',type=Path)
    a=p.parse_args();a.out.mkdir();rng=random.Random(19841);measurements=[]
    for n,lane in [(12,'collision'),(20,'collision'),(15,'elliptic'),(24,'elliptic'),(12,'rational'),(21,'mutation')]:
        for replicate in range(5):
            modes=['baseline','optimized'];rng.shuffle(modes);pair={}
            for mode in modes:
                folder=a.out/f'{lane}_{n}_{replicate}_{mode}'; binary=getattr(a,mode)
                cmd=[str(binary),'--d',str(n),'--lane',lane,'--seed',str(73841+replicate),
                    '--attempts',str(a.attempts),'--seconds','600','--closure-seconds','15','--cap','1000000',
                    '--keep-isotrivial','--out',str(folder)]
                if lane=='mutation':cmd+=['--mutation-seeds',str(a.mutations)]
                if a.words:cmd+=['--words',str(a.words)]
                started=time.perf_counter();subprocess.run(cmd,check=True,capture_output=True,text=True);wall=time.perf_counter()-started
                checked=audit_guided_native(folder,a.auditor,native_components=True)
                with (folder/'COUNTS.tsv').open() as f:counts=next(csv.DictReader(f,delimiter='\t'))
                with (folder/'STAGES.tsv').open() as f:stages=next(csv.DictReader(f,delimiter='\t'))
                assert int(counts['constructed'])==a.attempts
                pair[mode]=dict(wall=wall,counts=counts,stages=stages,leads_sha256=sha(folder/'LEADS.tsv'),
                    resolutions_sha256=sha(folder/'RESOLUTIONS.tsv'),independent=checked['independent_native'])
            for field in ['constructed','scope','transitive','passport','short_reject','duplicate','closed','unresolved']:
                assert pair['baseline']['counts'][field]==pair['optimized']['counts'][field],(n,lane,field)
            assert pair['baseline']['stages']['rng_final']==pair['optimized']['stages']['rng_final']
            assert pair['baseline']['leads_sha256']==pair['optimized']['leads_sha256']
            assert pair['baseline']['resolutions_sha256']==pair['optimized']['resolutions_sha256']
            measurements.append(dict(degree=n,lane=lane,replicate=replicate,results=pair))
        rows=[v for v in measurements if v['degree']==n and v['lane']==lane]
        print(n,lane,statistics.mean(v['results']['baseline']['wall']/v['results']['optimized']['wall'] for v in rows),flush=True)
    # Different intended target: exclude proved isotrivial families before BFS.
    # Do not call the resulting count difference a same-law kernel speedup.
    filtered=[]
    for rep in range(3):
        dest=a.out/f'filtered_mutation21_{rep}'
        start=time.perf_counter()
        cmd=[str(a.optimized),'--d','21','--lane','mutation','--seed',str(73841+rep),
            '--attempts',str(a.attempts),'--seconds','600','--closure-seconds','15','--cap','1000000',
            '--mutation-seeds',str(a.mutations),'--out',str(dest)]
        if a.words:cmd+=['--words',str(a.words)]
        subprocess.run(cmd,check=True,capture_output=True,text=True)
        wall=time.perf_counter()-start
        checks=audit_guided_native(dest,a.auditor,native_components=True)
        with (dest/'STAGES.tsv').open() as f:stages=next(csv.DictReader(f,delimiter='\t'))
        exclusions=audit_exclusions(dest/'ISOTRIVIAL.tsv',int(stages['isotrivial_outside_target']))
        filtered.append(dict(replicate=rep,wall=wall,independent=checks['independent_native'],exclusions=exclusions))
    summary=[];boot=random.Random(981)
    for n,lane in [(v['degree'],v['lane']) for v in measurements if v['replicate']==0]:
        ratios=[v['results']['baseline']['wall']/v['results']['optimized']['wall'] for v in measurements if v['degree']==n and v['lane']==lane]
        ci=sorted(statistics.mean(boot.choices(ratios,k=len(ratios))) for _ in range(10000))
        summary.append(dict(degree=n,lane=lane,mean_speedup=statistics.mean(ratios),bootstrap95=[ci[249],ci[9749]]))
    write(a.out/'CHECKS.json',dict(passed=True,baseline_sha256=sha(a.baseline),binary_sha256=sha(a.optimized),audit_binary_sha256=sha(a.auditor),
        compiled_source_files={str(v.relative_to(a.source)):sha(v) for v in (a.source/'src').rglob('*.rs')},
        matched=measurements,summaries=summary,filtered_target=filtered,
        scope='Busy-Mac fixed streams, not speed of finding a non-isotrivial example'))
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':main()
