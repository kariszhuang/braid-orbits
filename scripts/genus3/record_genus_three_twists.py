#!/usr/bin/env python3
"""Seal C209 finite evidence, negative-control tests and controller reductions."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
from fractions import Fraction
import json
from math import factorial
import os
from pathlib import Path
import subprocess
import sys

from run_genus_three import sha,write
from audit_guided_native import audit_guided_native


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('chords',type=Path);p.add_argument('saturated',type=Path);p.add_argument('fixedpoint',type=Path)
    p.add_argument('--producer',type=Path);p.add_argument('--pilot-out',type=Path)
    a=p.parse_args();project=Path(__file__).resolve().parents[2];out=project/'evidence/genus_three_twists_20261022'
    assert bool(a.producer)==bool(a.pilot_out)
    chords=json.loads((a.chords/'RESULT.json').read_text());saturated=json.loads((a.saturated/'RESULT.json').read_text());fixed=json.loads((a.fixedpoint/'CHECKS.json').read_text())
    assert chords['passed']and saturated['passed']and fixed['passed']
    eleven=[c for c in chords['cases']if c['degree']==11]
    assert len(eleven)==2 and sum(c['coverage']['qualified']for c in eleven)==115390
    assert all(c['coverage']['pending']==0 for c in chords['cases'])
    assert all(c['pending']==0 for c in saturated['cases'])and sum(c['candidate_states']for c in saturated['cases'])==44800
    for c in chords['cases']:
        folder=a.chords/f'{c["mode"]}_{c["degree"]}'
        assert sha(folder/'EXACT_AUDIT.json')==c['audit_sha256']
        audit=json.loads((folder/'EXACT_AUDIT.json').read_text())
        assert audit['native']['independent_native']['nonsmooth']==c['coverage']['qualified']
    for c in saturated['cases']:
        folder=a.saturated/('special_'+'_'.join(map(str,c['special'])))
        for name,h in c['files_sha256'].items():assert sha(folder/name)==h
        assert c['native']['independent_native']['nonsmooth']==c['candidate_states']
    env=dict(os.environ,G3_CHORD_BINARY=str(a.chords/'chords'),G3_CHORD_AUDITOR=str(a.chords/'independent_auditor'),G3_SATURATED_ROOT=str(a.saturated))
    tests=[]
    for name in ('test_genus_three_chords.py','test_genus_three_saturated.py','test_genus_three_twists.py','test_genus_three.py','test_watch_genus_three.py'):
        cmd=[sys.executable,'-m','unittest','discover','-s',str(project/'scripts/genus3'),'-p',name,'-v']
        r=subprocess.run(cmd,check=True,capture_output=True,text=True,env=env)
        assert 'skipped'not in r.stderr,'required native fixtures must not be skipped'
        tests.append(dict(command=cmd,stdout=r.stdout,stderr=r.stderr,passed=True))
    # Unsigned Stirling recurrence counts ALL third permutations with3 cycles.
    counts=[1]+[0]*11
    for n in range(1,12):counts=[0]+[counts[k-1]+(n-1)*counts[k]for k in range(1,12)]
    class_size=counts[3];reduction=[]
    for case in eleven:
        c=case['coverage'];baseline=len(c['frames'])*class_size;actual=c['constructed']
        reduction.append(dict(mode=case['mode'],all_three_cycle_permutations=class_size,frames=len(c['frames']),naive_full_frame_population=baseline,
            cyclic_forest_population=actual,exact_population_ratio=str(Fraction(baseline,actual)),
            qualification='candidate-population reduction,NOT a measured time or discovery multiplier'))
    source_names=['src/bin/genus_three_chords.rs','scripts/genus3/run_genus_three_chords.py','scripts/genus3/verify_genus_three_chords.py','scripts/genus3/test_genus_three_chords.py',
        'scripts/genus3/explore_genus_three_saturated.py','scripts/genus3/test_genus_three_saturated.py','scripts/genus3/verify_genus_three_twists.py','scripts/genus3/test_genus_three_twists.py',
        'scripts/genus3/run_genus_three.py','scripts/genus3/test_genus_three.py','scripts/genus3/record_genus_three_twists.py',
        'evidence/genus_three_twists_20261022/MATHEMATICS.md','evidence/genus_three_twists_20261022/README.md']
    pilot=None
    if a.producer:
        cmd=[str(a.producer),'--d','12','--lane','collision','--index-profile','1,9,9,9','--seed','2091207','--attempts','1000000','--seconds','30',
            '--trace-early','1000','--closure-seconds','60','--cap','1000000','--out',str(a.pilot_out)]
        subprocess.run(cmd,check=True,capture_output=True,text=True)
        checked=audit_guided_native(a.pilot_out,a.chords/'independent_auditor',native_components=True)
        write(a.pilot_out/'CHECKS.json',checked)
        pilot=dict(command=cmd,binary_sha256=sha(a.producer),checks_sha256=sha(a.pilot_out/'CHECKS.json'),
            counts=checked['construction'],native=checked['independent_native'],scope='bounded sampled open12 profile,not a census or a promised discovery')
    result=dict(passed=True,goal_found=False,new_nonisotrivial_genus3=0,new_components=0,pilot=pilot,
        degree11_necessary_negative_states=115390,new_degree15_negative_states=44800,fixedpoint_sampled_witnesses=fixed['records'],
        pure_math=['fixed-point twist lemma; whole odd-degree half-rank family','balanced single-cycle partition theorem','general support bounds'],
        exact_scope='whole identity-free source-genus3 degree11 via scalar proof and two exact slices; three pure ordinary3^5 degree15 passports; not a global genus3 exclusion',
        future_degree14=dict(indices=[6,6,7,13],necessary_index7_type='2^7',fixed_point_product_bound=14,
            proposed_involution_thirds=factorial(14)//(2**6*factorial(6)*factorial(2)),
            status='mathematical next-target arithmetic only,not enumerated'),
        raw_roots=dict(chords=str(a.chords.resolve()),saturated=str(a.saturated.resolve()),fixedpoint=str(a.fixedpoint.resolve())),
        raw_receipts={str(q.resolve()):sha(q)for q in (a.chords/'RESULT.json',a.chords/'TESTS.json',a.saturated/'RESULT.json',a.saturated/'TESTS.json',a.fixedpoint/'CHECKS.json')},
        chord_cases=chords['cases'],saturated_cases=saturated['cases'],population_reduction=reduction,tests=tests,
        source_sha256={name:sha(project/name)for name in source_names},human_mathematical_review='not recorded; elementary proofs retained for review; no new geometric classifier rule')
    write(out/'CHECKS.json',result)
    exclusions=dict(schema=1,source_genus=3,identity_free=True,excluded_degrees=[11],
        closed_index_budget=dict(degree=13,indices=[6,6,6,12],reason='entire fixed-point half-rank family,not just involutions'),
        focused_open_budget=dict(degree=12,indices=[1,9,9,9]),
        mathematics=dict(path=str(out/'MATHEMATICS.md'),sha256=sha(out/'MATHEMATICS.md')),
        checks=dict(path=str(out/'CHECKS.json'),sha256=sha(out/'CHECKS.json')),
        finite_degree11_audits=[dict(path=str(a.chords/f'{c["mode"]}_11/EXACT_AUDIT.json'),sha256=c['audit_sha256'])for c in eleven])
    write(out/'EXCLUSIONS.json',exclusions)
    print(json.dumps(dict(passed=True,degree11_negative_states=115390,degree15_negative_states=44800,population_reduction=reduction),indent=2))


if __name__=='__main__':main()
