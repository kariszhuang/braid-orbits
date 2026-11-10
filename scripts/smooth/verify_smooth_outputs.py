#!/usr/bin/env python3
"""Independent literal closure and SymPy group orders for a frozen smooth run."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import csv
import hashlib
import json
from pathlib import Path

from odd_order_spin import parse_seed
from run_smooth_mission import audit_component, audit_job, sha


def verify(root):
    import sympy
    from sympy.combinatorics import Permutation, PermutationGroup
    plan=json.loads((root/'PLAN.json').read_text())
    for name,digest in plan['source_hashes'].items():
        if sha(root/'source'/name)!=digest:raise ValueError('frozen source mismatch: '+name)
    if sha(root/'smooth_orbits')!=plan['binary_sha256']:raise ValueError('binary changed')
    jobs=[audit_job(root/f"A{3*j['m']+j['fixed']}_{j['shape']}") for j in plan['jobs']]
    folders=[root/'reported_A9',*sorted(root.glob('A*/component_*')),*sorted(root.glob('block_family_*'))]
    checked=[];components={}
    for folder in folders:
        if not (folder/'SUMMARY.tsv').is_file():continue
        audit=audit_component(folder);s=audit['summary'];n=int(s['degree'])
        t=parse_seed(s['seed'],n)
        g=PermutationGroup([Permutation(list(p)) for p in t]);order=int(g.order())
        if order!=int(s['group_order']):raise ValueError('independent group order mismatch')
        digest=hashlib.sha256(repr(audit['state_keys']).encode()).hexdigest()
        item=dict(folder=str(folder.relative_to(root)),degree=n,genus=int(s['genus']),
                  states=int(s['states']),status=s['status'],full_group=s['full_group']=='true',
                  independent_group_order=order,component_id=digest,seed=s['seed'],
                  states_sha256=sha(folder/'STATES.tsv'),summary_sha256=sha(folder/'SUMMARY.tsv'))
        checked.append(item)
        if audit['independent_six_move_closure']:components.setdefault(digest,item)
    controls=json.loads((root/'CONTROLS.json').read_text())
    if sum(x['states'] for x in controls['weak_catalog'])!=141 or any(x['smooth_states'] for x in controls['weak_catalog']):raise ValueError('weak controls')
    return dict(root=str(root),sympy_version=sympy.__version__,binary_sha256=plan['binary_sha256'],
                plan_sha256=sha(root/'PLAN.json'),checks_sha256=sha(root/'CHECKS.json'),
                completion=all(j['coverage']['fully_resolved']=='true' for j in jobs),
                jobs=jobs,checks=checked,distinct_closed_components=list(components.values()),
                additional_components_beyond_reported=len(components)-1,
                known_weak_states_checked=141,known_weak_smooth_states=0,
                warning='Different generated groups are explicitly distinguished from ambient A_d/S_d; no weak discovery or literature novelty claim.')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',required=True,type=Path)
    p.add_argument('--out',required=True,type=Path);args=p.parse_args()
    result=verify(args.root.resolve());result['verifier_sha256']=sha(Path(__file__))
    with args.out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(dict(complete=result['completion'],distinct_components=len(result['distinct_closed_components']),
                         additional=result['additional_components_beyond_reported']),indent=2))


if __name__=='__main__':main()
