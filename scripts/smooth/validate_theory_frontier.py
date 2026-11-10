"""Bounded theory/discovery checks and paired verified-bank workload timings.

This is not the delayed unbounded search. Every pilot receives independent
closure/connectivity audit; caps remain pending. One CPU/shared campaign lock.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import csv
import fcntl
import json
from pathlib import Path
import random
import shutil
import subprocess
import time

import sympy
from audit_guided_native import audit_guided_native,known_keys
from benchmark_forest_second import interval
from profile_smooth_theory import old_profiles,profile,structural_summary
from run_smooth_degree9 import sha
from verified_smooth_bank import VerifiedSmoothBank

PROJECT=Path(__file__).resolve().parents[2]
PRIOR=Path('data/discovery_design_final_checks_20261013')


def save(path,data):path.write_text(json.dumps(data,indent=2)+'\n')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',type=Path,required=True);p.add_argument('--auditor',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--repeats',type=int,default=7)
    args=p.parse_args();root=args.out.resolve();root.mkdir()
    if not Path('data').is_dir() or shutil.disk_usage(root).free<30*2**30:raise ValueError('T7 guard')
    lock=(PROJECT/'.run-control/.campaign.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    binary=args.binary.resolve();auditor=args.auditor.resolve()
    env=__import__('os').environ.copy();env.update(DESIGN_BINARY=str(binary),GUIDED_AUDIT_BINARY=str(auditor))
    test_logs=[]
    for name in ['test_theory_frontier','test_discovery_design','test_guided_native','test_forest_fast','test_forest_launchd','test_forest_recovery']:
        r=subprocess.run([__import__('sys').executable,str(_scripts_path.find(f'{name}.py'))],env=env,check=True,capture_output=True,text=True)
        path=root/f'{name}.log';path.write_text(r.stdout+r.stderr);test_logs.append(dict(name=name,sha256=sha(path)))
        print('PASS',name,flush=True)
    text=subprocess.check_output([str(binary),'--check-controls',str(PRIOR/'controls.tsv')],text=True)
    (root/'CONTROLS.log').write_text(text);print(text.strip(),flush=True)
    old=old_profiles();save(root/'THEORY_494.json',dict(components=old,summary=structural_summary(old)))
    print('OLD_THEORY',json.dumps(structural_summary(old),separators=(',',':'))[:1800],flush=True)
    bank=VerifiedSmoothBank(root/'VERIFIED_BANKS',PRIOR/'known_states.tsv',sha(PRIOR/'known_states.tsv'))
    old_keys=known_keys(PRIOR/'known_states.tsv');new={};pilots=[];leads=[]

    def job(folder,n,seed,attempts,bank_path,seconds=8):
        command=[binary,'--d',n,'--lane','quotient','--frame-reuse','1','--seed',seed,
            '--seconds',seconds,'--attempts',attempts,'--cap','3000','--closure-seconds','1',
            '--compact-frames','--export-negative-cache','--known-bank',bank_path,'--out',folder]
        began=time.monotonic();subprocess.run(list(map(str,command)),check=True,stdout=subprocess.DEVNULL)
        producer=time.monotonic()-began;check=audit_guided_native(folder,auditor,collect_records=True)
        total=time.monotonic()-began;records=check.pop('_closed_records')
        save(folder/'CHECKS.json',check);(folder/'VERIFIED').write_text('Independent bounded validation\n')
        return check,records,dict(producer=producer,total=total)

    for n in [12,15,16,18,20,21,24,25]:
        folder=root/f'pilot_{n}_quotient';snapshot,_=bank.snapshot(n)
        check,records,timing=job(folder,n,2026101300+n,100,snapshot,5)
        for item in check['closed_identities']:
            key=tuple(tuple(p) for p in item['key'])
            if (n,key) in old_keys:continue
            identity=json.dumps([n,item['key']],separators=(',',':'))
            if identity in new:continue
            value=profile(item['summary']['seed'],n)
            if value['group_order']!=int(item['summary']['group_order']) or value['genus']!=int(item['summary']['genus']):raise ValueError('independent profile mismatch')
            new[identity]=dict(value,states=int(item['summary']['states']),raw_folder=str(folder/item['folder']),
                summary_sha256=sha(folder/item['folder']/'SUMMARY.tsv'),
                table_sha256=sha(folder/item['folder']/'STATES.tsv'),checks_sha256=sha(folder/'CHECKS.json'))
        if records:bank.publish(folder,check,records)
        with (folder/'RESOLUTIONS.tsv').open() as stream:
            for row in csv.DictReader(stream,delimiter='\t'):
                if row['status'] in ['capped','timeout']:
                    leads.append(dict(degree=n,seed=row['seed'],source=str(folder),candidate=row['candidate'],checks_sha256=sha(folder/'CHECKS.json')))
        pilots.append(dict(degree=n,counts=check['construction'],known_aliases=check['known_aliases'],**timing))
        print('PILOT',json.dumps(pilots[-1]),flush=True)

    # The measured degree12 pilot is resolved before testing repeated work.
    # Same frozen proposal seed for both arms; reuse is the intended treatment.
    rng=random.Random(2026101321);timings=[]
    expanded,_=bank.snapshot(12)
    for repeat in range(args.repeats):
        order=['static','shared'];rng.shuffle(order);pair={};reports={};folders={}
        for variant in order:
            folder=root/f'bank_{repeat}_{variant}';folders[variant]=folder
            check,records,timing=job(folder,12,2026101312,100,
                PRIOR/'known_states.tsv' if variant=='static' else expanded,30)
            reports[variant]=check;pair[variant]=timing
        for field in ['constructed','scope','transitive','passport','short_reject']:
            if reports['static']['construction'][field]!=reports['shared']['construction'][field]:raise ValueError('bank changed proposal population')
        if sha(folders['static']/'FRAMES.tsv')!=sha(folders['shared']/'FRAMES.tsv'):raise ValueError('bank changed quotient frames')
        if int(reports['shared']['construction']['closed']) or int(reports['shared']['construction']['unresolved']):raise ValueError('verified bank did not retire repeated components')
        if reports['shared']['known_aliases']==0:raise ValueError('empty bank benchmark treatment')
        timings.append(pair);print('BANK_TIMING',repeat,json.dumps(pair),flush=True)

    entries=list(new.values());save(root/'CATALOG.json',dict(components=entries,summary=structural_summary(entries),sympy_version=sympy.__version__))
    # Immutable launch input = preceding494 + the independently disjoint new sets.
    known=root/'known_states.tsv';known.write_bytes((PRIOR/'known_states.tsv').read_bytes())
    with known.open('a') as output:
        for e in entries:
            with (Path(e['raw_folder'])/'STATES.tsv').open() as stream:
                for row in csv.DictReader(stream,delimiter='\t'):output.write(str(e['degree'])+'\t'+row['tuple']+'\n')
    expected=58956+sum(e['states'] for e in entries)
    if len(known_keys(known))!=expected:raise ValueError('component/state union overlap')
    (root/'controls.tsv').write_bytes((PRIOR/'controls.tsv').read_bytes())
    prior_leads=json.loads((PRIOR/'saved_leads.json').read_text())
    save(root/'saved_leads.json',prior_leads+leads)
    result=dict(passed=True,design_passed=True,theory_passed=True,controls=66,control_states=25192,
        binary_sha256=sha(binary),audit_binary_sha256=sha(auditor),known_bank_sha256=sha(known),
        expanded_known_components=494+len(entries),expanded_known_states=expected,
        original_catalog_root='data/discovery_design_catalog_20261013',
        theory_494_sha256=sha(root/'THEORY_494.json'),new_catalog_sha256=sha(root/'CATALOG.json'),
        pilots=pilots,tests=test_logs,new_summary=structural_summary(entries),
        saved_partial_occurrences=len(leads),saved_leads_sha256=sha(root/'saved_leads.json'),
        bank_benchmark=dict(repeats=args.repeats,timings=timings,
            pipeline=interval([t['static']['total']/t['shared']['total'] for t in timings]),
            scope='same100-proposal mixed degree12 quotient stream after independent verified publication; not general discovery throughput'),
        limitations='Same generalized imprimitive mechanism; full independent closures only, capped states unresolved. No primitive/nonisotrivial classification or literature novelty.')
    save(root/'RESULTS.json',result);print('RESULTS='+str(root/'RESULTS.json'),flush=True);lock.close()

if __name__=='__main__':main()
