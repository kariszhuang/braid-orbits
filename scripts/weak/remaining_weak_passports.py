#!/usr/bin/env python3
"""Exact coverage map for specified XXYY/XXXY passport profiles, not all tuples.

Reuse the production scanner's passport generator and validated exclusions.
An unexcluded passport is unclassified by this map; it is NOT evidence that
no previous random candidate had that passport. Sampled negative results never
become complete exclusions here.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys

PROJECT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(PROJECT))
from run_campaign import finite_passport_exclusions, structured_slices


def build(degrees,genera,excluded):
    result=[]
    for group in ('A','S'):
        for row in structured_slices(group,degrees,genera,quotient=True):
            key=(group,row['d'],tuple(sorted(row['passport'].split('/'))))
            source=excluded.get(key)
            result.append(dict(group=group,degree=row['d'],genus=row['genus'],passport=row['passport'],
                mechanism=row['mechanisms'],coverage='resolved_selected_passport' if source else 'not_classified_by_complete_coverage_ledger',
                evidence=source,raw_fixed_first_pairs=str(row['total']),
                estimated_centralizer_rows=str(row['estimated_pairs']),
                estimate_is_not_actual_scan_size=True))
    return sorted(result,key=lambda r:(r['coverage']=='resolved_selected_passport',r['degree'],
                                       int(r['estimated_centralizer_rows']),r['group'],r['passport']))


def selftest():
    rows=build([9],[1],{})
    assert rows and all(r['genus']==1 and r['coverage'].startswith('not_classified') for r in rows)
    row=rows[0];key=(row['group'],9,tuple(sorted(row['passport'].split('/'))))
    with_exclusion=build([9],[1],{key:'verified_fixture'})
    assert sum(r['coverage']=='resolved_selected_passport' for r in with_exclusion)==1
    assert all('component' not in r for r in rows)
    print('remaining passport scope checks passed')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--degrees',default='9,10,11,12,13,14,15,16,17')
    parser.add_argument('--genera',default='0,1,2,3')
    parser.add_argument('--out',type=Path);parser.add_argument('--selftest',action='store_true')
    args=parser.parse_args()
    if args.selftest:selftest();return
    if args.out is None:parser.error('out required')
    degrees=sorted(set(map(int,args.degrees.split(','))));genera=sorted(set(map(int,args.genera.split(','))))
    if not degrees or min(degrees)<4 or max(degrees)>17 or not genera or min(genera)<0 or max(genera)>12:
        parser.error('degrees 4..17 and genera 0..12')
    exclusions=finite_passport_exclusions();rows=build(degrees,genera,exclusions)
    sources={str(Path(__file__)):hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             str(PROJECT/'scripts/weak/run_campaign.py'):hashlib.sha256((PROJECT/'scripts/weak/run_campaign.py').read_bytes()).hexdigest()}
    for source in sorted({r['evidence'] for r in rows if r['evidence']}):
        path=Path(source)
        if not path.is_absolute():path=PROJECT/path
        if path.is_dir():
            for file in sorted(path.glob('attempt_*/SUMMARY.txt')):
                if (file.parent/'VERIFIED').exists():sources[str(file)]=hashlib.sha256(file.read_bytes()).hexdigest()
        else:sources[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
    document=dict(degrees=degrees,genera=genera,passports=rows,
                  statuses=dict(Counter(r['coverage'] for r in rows)),source_sha256=sources,
                  scope='all generated nonidentity XXYY/XXXY passports in requested degrees/genera; not arbitrary four-types, tuple enumeration, or universal no-counterexample coverage',
                  sampling_status='not inferred; construction exclusions are not whole-passport exclusions')
    with args.out.open('x') as stream:json.dump(document,stream,indent=2);stream.write('\n')


if __name__=='__main__':main()
