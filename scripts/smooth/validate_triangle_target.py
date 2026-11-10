#!/usr/bin/env python3
"""Focused prelaunch proof checks for the wide, three-shard TRIANGLE search."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import csv
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import time

from audit_guided_native import audit_component_native, audit_guided_native
from odd_order_spin import cycle_lengths, parse_seed
from run_smooth_degree9 import sha
from run_smooth_mission import values
from run_s18_linked_roots import state_key

PROJECT=Path(__file__).resolve().parents[2]
CATALOG=PROJECT/'evidence/smooth_positive_catalog_20261012/CATALOG.json'

def save(path,value):
    path=Path(path);tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(value,indent=2)+'\n');tmp.replace(path)

def read_rows(path):
    with Path(path).open(newline='') as stream:return list(csv.DictReader(stream,delimiter='\t'))

def audit_triangle_leads(folder,degree,lo,hi,min_genus,exact_profile=None,exact_genus=None):
    checked=0
    with (Path(folder)/'LEADS.tsv').open() as stream:
        for line in stream:
            _,raw=line.rstrip('\n').split('\t',1);t=parse_seed(raw,degree)
            v=values(t);ix=sorted(degree-len(cycle_lengths(p)) for p in t)
            if not v['smooth'] or v['genus']<min_genus:raise ValueError('lead fails smoothness/genus filter')
            if exact_genus is not None and v['genus']!=exact_genus:raise ValueError('lead outside exact source-genus target')
            if ix[1]+ix[2]>degree-1 or ix[0]+ix[3]<degree:raise ValueError('lead outside TRIANGLE budget')
            if not lo<=ix[0]<=hi:raise ValueError('lead in wrong minimum-index shard')
            if exact_profile is not None and ix!=list(exact_profile):raise ValueError('lead outside exact target passport')
            high=ix[3]
            if not any(degree-len(cycle_lengths(p))==high and 1 not in cycle_lengths(p) for p in t):
                raise ValueError('lead has no fixed-point-free maximum-index entry')
            if any(vv==tuple(range(degree)) for vv in t):raise ValueError('identity entry')
            checked+=1
    return checked

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--binary',type=Path,required=True);ap.add_argument('--auditor',type=Path,required=True)
    ap.add_argument('--prior-run',type=Path,required=True);ap.add_argument('--compiled-source',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--attempts',type=int,default=6000);ap.add_argument('--pilot-seconds',type=int,default=45)
    args=ap.parse_args();binary=args.binary.resolve();auditor=args.auditor.resolve();prior=args.prior_run.resolve()
    root=args.out.resolve();root.mkdir(parents=True,exist_ok=False)
    lock=(PROJECT/'.run-control/.campaign.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    try:
        controls=prior/'CONTROLS.tsv';full_bank=prior/'KNOWN_STATES.tsv'
        prior_validation=json.loads((prior/'VALIDATION.json').read_text())
        if not prior_validation.get('passed') or sha(full_bank)!=prior_validation['known_bank_sha256']:
            raise ValueError('prior verified-bank provenance failed')
        shutil.copy2(controls,root/'controls.tsv')
        parent_bank_sha=sha(full_bank);bank_rows=0;degree_mass={}

        catalog=json.loads(CATALOG.read_text())['components']
        control_rows=[line.split('\t',2) for line in controls.read_text().splitlines()]
        if len(catalog)!=66 or len(control_rows)!=66 or sum(int(r[1]) for r in control_rows)!=25192:
            raise ValueError('legacy positive-control mass changed')
        if [(int(r[0]),int(r[1]),r[2]) for r in control_rows]!=[
            (int(c['degree']),int(c['states']),c['representative']) for c in catalog]:
            raise ValueError('control rows differ from the bound 66-component catalog')

        exported=root/'controls-export.tsv'
        output=subprocess.check_output([str(binary),'--check-controls',str(controls),'--export-known',str(exported)],text=True)
        m=re.search(r'controls=(\d+) states=(\d+) all_closed=true all_compatible_keys=true',output)
        if not m or (int(m.group(1)),int(m.group(2)))!=(66,25192):raise ValueError('wide-binary control closure/key check failed')
        for item in catalog:
            d=int(item['degree']);ix=sorted(item['representative_values']['indices'])
            if ix[0]+ix[3]>d-1 or ix[1]+ix[2]<=d-1:
                raise ValueError(f'control {item["id"]} is not in the STAR regime')
        if sum(1 for _ in exported.open())!=25192:raise ValueError('checked control-state export mass')
        # This is only a duplicate-skip subset, not a complete historical
        # catalog. Every included state is a verified STAR state, so it cannot
        # equal a TRIANGLE candidate; omitted rows can only cause repeated work.
        shutil.copy2(exported,root/'known_states.tsv')
        with (root/'known_states.tsv').open() as stream:
            for line in stream:
                d=int(line.split('\t',1)[0]);bank_rows+=1;degree_mass[d]=degree_mass.get(d,0)+1
        if bank_rows!=25192:raise ValueError('small verified STAR skip-bank mass')

        union=set();native_seconds=0.0
        for i,item in enumerate(catalog):
            began=time.monotonic();check=audit_component_native(item['raw_folder'],auditor)
            native_seconds+=time.monotonic()-began
            summary=check['summary'];keys=check['state_keys'];n=int(item['degree'])
            if not check['independent_six_move_closure'] or len(keys)!=int(item['states']):raise ValueError(f'native control not closed: {item["id"]}')
            if check['table_sha256']!=item['table_sha256'] or int(summary['states'])!=int(item['states']):
                raise ValueError(f'native control binding mismatch: {item["id"]}')
            current={(n,key) for key in keys}
            if union & current:raise ValueError('the 66 controls are not distinct components')
            union.update(current)
            if i%16==0:print(f'NATIVE_CONTROLS {i+1}/66',flush=True)
        if len(union)!=25192:raise ValueError('native control-state union mass')

        pilots=[];shards=[(16,7,7,8,3,'centralizer',(7,7,7,15)),
            (17,6,12,8,3,'mixed',None),(27,13,31,8,3,'joint',None),
            (16,7,7,8,3,'strong',(7,7,7,15)),
            (16,7,7,8,3,'pairwise_low',(7,7,7,15)),
            (16,7,7,8,3,'middle_repair',(7,7,7,15))]
        env=os.environ.copy();env['DESIGN_BINARY']=str(binary);env['GUIDED_AUDIT_BINARY']=str(auditor)
        for j,(d,lo,hi,strength,exact_genus,method,exact_profile) in enumerate(shards):
            folder=root/f'pilot_{j}_d{d}';command=[str(binary),'--d',str(d),'--lane','triangle',
                '--triangle-min-genus','3',
                '--triangle-low-min',str(lo),
                '--triangle-low-max',str(7 if exact_profile else 63),
                '--triangle-passport-min',str(lo),'--triangle-passport-max',str(hi),
                '--triangle-strength-divisor',str(strength),'--triangle-third-method',method,
                '--frame-reuse','1','--seed',str(2026103100+j),
                '--seconds',str(args.pilot_seconds),'--attempts',str(args.attempts),'--closure-seconds','4',
                '--cap','20000','--full-closure-budget','--compact-frames','--known-bank',str(root/'known_states.tsv'),
                '--out',str(folder)]
            if exact_genus:command+=['--target-genus',str(exact_genus)]
            if exact_profile:command+=['--triangle-exact-passport',','.join(map(str,exact_profile))]
            subprocess.run(command,check=True,capture_output=True,text=True,timeout=args.pilot_seconds+240)
            independent=audit_guided_native(folder,auditor,native_components=True)
            save(folder/'CHECKS.json',independent);(folder/'VERIFIED').write_text('Focused TRIANGLE prelaunch pilot; every promoted closure independently checked.\n')
            checked_leads=audit_triangle_leads(folder,d,lo,hi,3,exact_profile,exact_genus or None)
            if checked_leads!=int(independent['construction']['scope']):
                # Scope is the set that reaches exact transitivity/passport screening;
                # leads are only starts that survive all inexpensive rejection tests.
                if checked_leads>int(independent['construction']['scope']):raise ValueError('pilot lead mass exceeds scope mass')
            pilots.append(dict(degree=d,passport_minimum=[lo,hi],third_method=method,strength_divisor=strength,
                exact_genus=exact_genus or None,construction=independent['construction'],
                passport_profiles=independent['triangle_filters']['passport_profiles'],
                exact_passport=exact_profile,
                accepted_leads_checked=checked_leads,closed_components=len(independent['closed_identities']),
                unresolved=independent['independent_native']['unresolved'],checks_sha256=sha(folder/'CHECKS.json')))
            print('TRIANGLE_PILOT',json.dumps(pilots[-1]),flush=True)

        # Exercise the packed-permutation upper endpoint without admitting a
        # closure candidate: genus 100 is a safe impossible pilot threshold.
        max_degree_job=root/'smoke_d63';subprocess.run([str(binary),'--d','63','--lane','triangle',
            '--triangle-min-genus','100','--triangle-low-min','2','--triangle-low-max','63',
            '--triangle-passport-min','13','--triangle-passport-max','31',
            '--triangle-strength-divisor','8','--frame-reuse','1','--seed','2026103163',
            '--seconds','5','--attempts','1','--closure-seconds','1','--cap','1000',
            '--out',str(max_degree_job)],check=True,capture_output=True,text=True,timeout=30)
        smoke=read_rows(max_degree_job/'COUNTS.tsv')
        if len(smoke)!=1 or int(smoke[0]['degree'])!=63:raise ValueError('wide degree-63 smoke did not run')

        # Exercise the launcher shard schedule over the full supported r1 range.
        env['PYTHONPATH']=os.pathsep.join(str(d) for d in sorted((PROJECT/'scripts').iterdir()) if d.is_dir())
        schedule=subprocess.run([sys.executable,'-m','unittest','test_structural_novelty.StructuralNoveltyTests.test_triangle_schedule_skips_reported_exclusions_and_covers_higher_r1_bands'],
            cwd=PROJECT,capture_output=True,text=True,env=env,timeout=120)
        if schedule.returncode:raise RuntimeError('triangle shard schedule test failed\n'+schedule.stdout+schedule.stderr)
        (root/'SHARD_TEST.log').write_text(schedule.stdout+schedule.stderr)
        compiled_source=args.compiled_source.resolve()
        compiled_hashes={str(p.relative_to(compiled_source)):sha(p) for p in sorted((compiled_source/'src').rglob('*.rs'))}
        launch_scripts=['scripts/smooth/run_forest_fast.py','scripts/smooth/audit_guided_native.py',
            'scripts/smooth/verified_smooth_bank.py','scripts/smooth/validate_triangle_target.py',
            'scripts/smooth/test_structural_novelty.py']
        launch_script_hashes={rel:sha(PROJECT/rel) for rel in launch_scripts}
        result=dict(passed=True,design_passed=True,native_components_passed=True,
            binary_sha256=sha(binary),audit_binary_sha256=sha(auditor),controls=66,control_states=25192,
            compiled_source_hashes=compiled_hashes,launch_script_hashes=launch_script_hashes,
            known_bank_sha256=sha(root/'known_states.tsv'),known_bank_source_sha256=parent_bank_sha,
            known_bank_rows=bank_rows,known_bank_degree_mass=degree_mass,
            known_bank_provenance='A verified subset of 25,192 states from the 66 independently closed STAR control components. STAR/TRIANGLE index regimes are disjoint in positive genus, so it cannot filter a target TRIANGLE state; omitted historical rows cause only redundant work.',
            full_native_control_audit=dict(components=66,states=len(union),seconds=native_seconds,distinct=True),
            wide_binary_controls_output=output.strip(),triangle_pilots=pilots,
            degree_63_smoke=dict(degree=63,constructed=smoke[0]['constructed'],
                genus_reject=read_rows(max_degree_job/'TRIANGLE_FILTERS.tsv')[0]['genus_reject']),
            shard_test_sha256=sha(root/'SHARD_TEST.log'),
            scope='Lane0 targets exact source genus 3 and the degree-16 passport (7,7,7,15); lanes1/2 sample exact source-genus-3 TRIANGLE passports for d>=17 in disjoint r1 bands. Constructor sampling is biased, not a census. Full six-move closure and independent native audit remain mandatory.',
            triangle_theory='For a source-genus-g TRIANGLE passport, r1>=2g and h<=d-1 imply the elementary degree floor; reviewed C12/C13 additionally require h>=floor(d/2)+2g+1. Exact genus-3 TRIANGLE work excludes degrees<=15; the degree-16 table leaves (7,7,7,15) as the search target, with the three-involution subcase already excluded and the other reported computational rows still requiring certificate audit. The reviewed C12 theorem excludes a maximum-index entry of type 15*2^((d-15)/2) in genus 3. C08 says a full all-smooth TRIANGLE orbit in source genus>=2 is non-isotrivial; C09 requires a fixed-point-free maximum-index entry.')
        save(root/'RESULTS.json',result);print('VALIDATION='+str(root/'RESULTS.json'),flush=True)
    finally:
        lock.close()

if __name__=='__main__':main()
