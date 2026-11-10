#!/usr/bin/env python3
"""Matched outer-block/ordinary BFS measurements; identical seeded populations."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import csv
import json
from pathlib import Path
import random
import re
import statistics
import subprocess
import time

from audit_guided_native import audit_guided_native, audit_component_native
from run_genus_three import sha, write, mutation_seed_text


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['binary','auditor','controls','tests','source','out']:
        p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--catalog',type=Path)
    a=p.parse_args(); a.out.mkdir()
    tests=[]
    for name in ['genus_three_tests','forest_fast_tests','structural_novelty_tests','discovery_design_tests']:
        r=subprocess.run([str(a.tests),name],capture_output=True,text=True,check=True)
        tests.append(dict(name=name,output=r.stdout))
    control_dir=a.out/'controls'
    control=subprocess.check_output([str(a.binary),'--controls',str(a.controls),'--compare','--out',str(control_dir)],text=True)
    state_total=0; control_count=0
    for folder in sorted(control_dir.glob('control_*')):
        check=audit_component_native(folder,a.auditor)
        assert check['independent_six_move_closure']
        state_total+=int(check['summary']['states']);control_count+=1
    print(control.strip(),flush=True)
    order=random.Random(20261021); blocks=[]
    for block in range(7):
        modes=['ordinary','macro'];order.shuffle(modes);row={}
        for mode in modes:
            text=subprocess.check_output([str(a.binary),'--controls',str(a.controls),'--verifier',mode],text=True)
            row[mode]=float(re.search(r'closure_seconds=(\S+)',text).group(1))
        blocks.append(row)
    # Fixed work, not different numbers of time-limited random draws. Certifier
    # order may change the literal negative word, never the candidate stream.
    pilots=[]
    populations=[(12,'collision'),(15,'elliptic'),(12,'rational')]
    if a.catalog:
        (a.out/'MUTATION_SEEDS.tsv').write_text(mutation_seed_text(json.loads(a.catalog.read_text())))
        populations.append((15,'mutation'))
    for n,lane in populations:
        for rep in range(3):
            modes=['ordinary','macro'];order.shuffle(modes);pair={}
            for mode in modes:
                dest=a.out/f'pilot_{n}_{lane}_{rep}_{mode}'
                start=time.perf_counter()
                cmd=[str(a.binary),'--d',str(n),'--lane',lane,'--seed',str(9731+rep),
                    '--attempts','2000000','--seconds','600','--cap','100000','--closure-seconds','10',
                    '--verifier',mode,'--out',str(dest),'--keep-isotrivial']
                if lane=='mutation':cmd.extend(['--mutation-seeds',str(a.out/'MUTATION_SEEDS.tsv')])
                subprocess.run(cmd,check=True,capture_output=True,text=True)
                wall=time.perf_counter()-start
                checks=audit_guided_native(dest,a.auditor,native_components=True)
                write(dest/'CHECKS.json',checks);(dest/'VERIFIED').write_text('All persisted outcomes audited.\n')
                with (dest/'COUNTS.tsv').open() as stream:c=next(csv.DictReader(stream,delimiter='\t'))
                with (dest/'STAGES.tsv').open() as stream:s=next(csv.DictReader(stream,delimiter='\t'))
                assert int(c['constructed'])==2000000
                pair[mode]=dict(wall_seconds=wall,producer_seconds=float(c['seconds']),counts=c,stages=s,
                    leads_sha256=sha(dest/'LEADS.tsv'),independent=checks['independent_native'])
            assert pair['ordinary']['leads_sha256']==pair['macro']['leads_sha256']
            assert pair['ordinary']['stages']['rng_final']==pair['macro']['stages']['rng_final']
            for field in ['constructed','scope','transitive','passport','short_reject','duplicate','closed','unresolved']:
                assert pair['ordinary']['counts'][field]==pair['macro']['counts'][field]
            pilots.append(dict(degree=n,lane=lane,replicate=rep,results=pair))
            print(f'matched {n}/{lane}/{rep} verified',flush=True)
    ratios=[b['ordinary']/b['macro'] for b in blocks]
    boot=random.Random(812); means=[]
    for _ in range(10000):means.append(statistics.mean(boot.choices(ratios,k=len(ratios))))
    means.sort()
    receipt=dict(passed=True,binary_sha256=sha(a.binary),audit_binary_sha256=sha(a.auditor),
        compiled_source_files={str(v.relative_to(a.source)):sha(v) for v in (a.source/'src').rglob('*.rs')},
        controls=control_count,control_states=state_total,independent_all_six_moves=True,
        rust_tests=tests,controls_sha256=sha(a.controls),closure_blocks=blocks,
        closure_speedup=dict(mean_paired_ratio=statistics.mean(ratios),bootstrap95=[means[249],means[9749]],
            scope='known positive-control closure kernels only; not a discovery multiplier'),
        matched_search_pilots=pilots,limitations='biased proposals, includes repeats; caps/timeouts undecided; no non-isotriviality proof')
    write(a.out/'CHECKS.json',receipt)
    print(json.dumps(receipt['closure_speedup'],indent=2),flush=True)


if __name__=='__main__':main()
