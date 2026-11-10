#!/usr/bin/env python3
"""Paired fixed-proposal sieve comparison, within one reserved CPU slot.

Randomized AB/BA order; construction RNG/counts and positive/capped outcomes
must agree. Negative witnesses can change and are independently replayed.
This measures throughput, never the probability of finding the target.
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
from verify_genus_three_exclusions import audit_exclusions


def row(folder, name):
    with (folder/name).open() as stream:
        return next(csv.DictReader(stream, delimiter='\t'))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('build',type=Path)
    p.add_argument('base',type=Path)
    p.add_argument('out',type=Path)
    p.add_argument('--attempts',type=int,default=1000000)
    p.add_argument('--repeats',type=int,default=5)
    p.add_argument('--placement',choices=['early','late'],default='late')
    p.add_argument('--screen',choices=['cyclic','forest'],default='cyclic',
        help='forest compares opt-in whole-subset cuts against existing cyclic sieve')
    a=p.parse_args()
    if min(a.attempts,a.repeats)<=0: p.error('positive bounds required')
    a.out.mkdir()
    checkpoint=json.loads((a.base/'CHECKPOINT.json').read_text())
    words=Path(checkpoint['word_bank'])
    if sha(words)!=checkpoint['word_bank_sha256']: raise ValueError('word binding changed')
    rng=random.Random(213); results=[]; summaries=[]
    cases=[(12,'collision'),(16,'collision'),(25,'collision'),(33,'general-fiber'),(21,'mutation')]
    for n,lane in cases:
        for rep in range(a.repeats):
            order=['reference','cyclic']; rng.shuffle(order); pair={}
            for mode in order:
                folder=a.out/f'{n}_{lane}_{rep}_{mode}'
                cmd=[a.build/'bin'/('compact' if n<=25 else 'wide'),'--d',str(n),
                     '--lane',lane,'--seed',str(2026102213+rep),'--attempts',str(a.attempts),
                     '--seconds','600','--closure-seconds','15','--cap','1000000',
                     '--words',words,'--trace-early','1000','--out',folder]
                if a.screen=='forest':
                    if mode=='cyclic':cmd+=['--cyclic-forest-cuts']
                elif mode=='reference': cmd+=['--no-cyclic-screen']
                elif a.placement=='early':cmd+=['--cyclic-before-deep']
                if lane!='mutation': cmd+=['--literal-trie']
                if lane=='mutation': cmd+=['--mutation-seeds',a.base/'MUTATION_SEEDS.tsv']
                start=time.perf_counter()
                subprocess.run(list(map(str,cmd)),check=True,capture_output=True)
                wall=time.perf_counter()-start
                checked=audit_guided_native(folder,a.build/'bin/auditor',native_components=True)
                stages=row(folder,'STAGES.tsv')
                checked['isotrivial_exclusions']=audit_exclusions(folder/'ISOTRIVIAL.tsv',int(stages['isotrivial_outside_target']))
                write(folder/'CHECKS.json',checked)
                counts=row(folder,'COUNTS.tsv')
                if int(counts['constructed'])!=a.attempts: raise ValueError('wall safety limit hit; not a matched comparison')
                identities=[dict(key=x['key'],states=x['summary']['states']) for x in checked['closed_identities']]
                pair[mode]=dict(wall=wall,counts=counts,stages=stages,cyclic=row(folder,'CYCLIC.tsv'),
                    forest=row(folder,'FOREST_CUTS.tsv') if (folder/'FOREST_CUTS.tsv').exists() else None,
                    independent=checked['independent_native'],positive_identities=identities,
                    checks_sha256=sha(folder/'CHECKS.json'))
            for key in ('constructed','scope','passport','closed','unresolved'):
                if pair['reference']['counts'][key]!=pair['cyclic']['counts'][key]:raise ValueError('proposal/outcome mismatch: '+key)
            if pair['reference']['positive_identities']!=pair['cyclic']['positive_identities']:raise ValueError('positive component mismatch')
            if pair['reference']['stages']['rng_final']!=pair['cyclic']['stages']['rng_final']:raise ValueError('proposal RNG changed')
            results.append(dict(degree=n,lane=lane,replicate=rep,results=pair))
            print(n,lane,rep,round(pair['reference']['wall']/pair['cyclic']['wall'],3),flush=True)
        ratios=[r['results']['reference']['wall']/r['results']['cyclic']['wall'] for r in results if r['degree']==n and r['lane']==lane]
        bootstrap=sorted(statistics.mean(rng.choices(ratios,k=len(ratios))) for _ in range(10000))
        summaries.append(dict(degree=n,lane=lane,median_speedup=statistics.median(ratios),mean_speedup=statistics.mean(ratios),bootstrap95=[bootstrap[249],bootstrap[9749]],
            qualification='paired busy-Mac full-job throughput, not discovery gain'))
    write(a.out/'CHECKS.json',dict(passed=True,paired=results,summaries=summaries,
        binary_sha256=sha(a.build/'bin/wide'),compact_sha256=sha(a.build/'bin/compact'),auditor_sha256=sha(a.build/'bin/auditor'),
        word_bank=dict(path=str(words),sha256=sha(words)),attempts=a.attempts,repeats=a.repeats,placement=a.placement,screen=a.screen))
    print(json.dumps(summaries,indent=2))


if __name__=='__main__': main()
