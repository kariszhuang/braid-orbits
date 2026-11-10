"""Validate independent native closures against bound Python-checked catalogs.

Checks every state/edge in all837 already audited components. Seven paired
same-certificate timings compare complete Python and independent Rust audits.
No discovery/census rerun or change to a prepared/active campaign.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
from collections import Counter
import fcntl
import json
from pathlib import Path
import random
import shutil
import subprocess
import time

from audit_guided_native import audit_component_native,audit_guided_native,known_keys
from benchmark_forest_second import interval
from profile_smooth_actions import action_profile
from profile_smooth_theory import structural_summary
from run_smooth_mission import audit_component
from run_smooth_degree9 import sha

PROJECT=Path(__file__).resolve().parents[2]


def save(path,data):path.write_text(json.dumps(data,indent=2)+'\n')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--validation',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--binary',type=Path,required=True);p.add_argument('--auditor',type=Path,required=True)
    args=p.parse_args();root=args.out.resolve();root.mkdir();prior=args.validation.resolve()
    lock=(PROJECT/'.run-control/.campaign.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    binary=args.binary.resolve();auditor=args.auditor.resolve();source=root/'source';source.mkdir();hashes={}
    for path in [*sorted((PROJECT/'src').rglob('*.rs')),*sorted((PROJECT/'scripts').rglob('*.py')),
                 *sorted((PROJECT/'evidence/smooth_theory_frontier_20261013').glob('*'))]:
        rel=path.relative_to(PROJECT);target=source/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(path,target);hashes[str(rel)]=sha(target)
    save(root/'SOURCE_HASHES.json',hashes)
    env=__import__('os').environ.copy();env.update(DESIGN_BINARY=str(binary),GUIDED_AUDIT_BINARY=str(auditor))
    tests=[]
    for name in ['test_native_components','test_theory_frontier','test_discovery_design','test_guided_native','test_forest_fast','test_forest_launchd','test_forest_recovery']:
        result=subprocess.run([__import__('sys').executable,str(_scripts_path.find(f'{name}.py'))],env=env,check=True,capture_output=True,text=True)
        path=root/f'{name}.log';path.write_text(result.stdout+result.stderr);tests.append(dict(name=name,sha256=sha(path)))
        print('PASS',name,flush=True)
    text=subprocess.check_output([str(binary),'--check-controls',str(prior/'controls.tsv')],text=True)
    (root/'CONTROLS.log').write_text(text);print(text.strip(),flush=True)
    old=json.loads((prior/'THEORY_494.json').read_text())['components']
    new=json.loads((prior/'CATALOG.json').read_text())['components'];all_entries=old+new
    expected=known_keys(prior/'known_states.tsv');union=set();entries=[];began=time.monotonic()
    for i,e in enumerate(all_entries):
        folder=Path(e['raw_folder'])
        check=audit_component_native(folder,auditor);keys=check.pop('state_keys');n=e['degree']
        if check['table_sha256']!=e['table_sha256'] or not check['independent_six_move_closure']:raise ValueError('catalog audit binding')
        current={(n,key) for key in keys}
        if len(current)!=e['states'] or current&union:raise ValueError('distinct native component/state mass')
        union.update(current);entries.append(dict(e,action=action_profile(folder,e['table_sha256'])))
        if i%100==0:print('NATIVE_COMPONENTS',i+1,flush=True)
    if union!=expected or len(union)!=166056:raise ValueError('native lex-key union differs from complete Python catalog')
    all_native_seconds=time.monotonic()-began
    # A full independent audit of the old frozen73.7-second pipeline's input.
    job=prior/'pilot_21_quotient';began=time.monotonic();checked=audit_guided_native(job,auditor,native_components=True)
    native_job_seconds=time.monotonic()-began;old_check=json.loads((job/'CHECKS.json').read_text())
    for field in ['construction','components','closed_identities','known_aliases','input_banks','exported_negative_bank','sources']:
        if json.loads(json.dumps(checked[field]))!=old_check[field]:raise ValueError('native/Python full guided comparison: '+field)
    # Seven independent timed repetitions, randomized order, same largest
    # degree21 certificate. Every repeat also compares all canonical keys.
    selected=max([e for e in new if e['degree']==21],key=lambda e:e['states']);folder=Path(selected['raw_folder'])
    rng=random.Random(2026101337);timings=[]
    for repeat in range(7):
        order=['python','native'];rng.shuffle(order);reports={};pair={}
        for variant in order:
            began=time.monotonic();reports[variant]=(audit_component(folder) if variant=='python' else audit_component_native(folder,auditor))
            pair[variant]=time.monotonic()-began
        if reports['python']!=reports['native']:raise ValueError('same-certificate complete audit/key comparison')
        timings.append(pair);print('COMPONENT_TIMING',repeat,json.dumps(pair),flush=True)
    for name in ['known_states.tsv','controls.tsv','saved_leads.json']:shutil.copy2(prior/name,root/name)
    prior_result=json.loads((prior/'RESULTS.json').read_text())
    for name,digest in hashes.items():
        if sha(PROJECT/name)!=digest:raise ValueError('validation source changed while checking: '+name)
    groups=Counter((e['quotient_genus'],e['action']['pure_braid_fiber_size']) for e in entries)
    save(root/'THEORY_ALL.json',dict(components=entries,summary=structural_summary(entries),
        pure_braid_patterns=[dict(quotient_genus=h,pure_fiber=s,components=v) for (h,s),v in sorted(groups.items())]))
    result=dict(prior_result,binary_sha256=sha(binary),audit_binary_sha256=sha(auditor),native_components_passed=True,
        measured_constructor_binary_sha256=prior_result['binary_sha256'],
        measured_constructor_auditor_sha256=prior_result['audit_binary_sha256'],
        all_native_components=len(entries),all_native_states=len(union),all_native_seconds=all_native_seconds,
        native_python_canonical_union_equal=True,full_job_native_seconds=native_job_seconds,
        native_component_benchmark=dict(certificate=str(folder),states=selected['states'],timings=timings,
            audit=interval([t['python']/t['native'] for t in timings]),
            scope='complete audit of same1512-state degree21 certificate, startup/export/key checks included; not general discovery-rate gain'),
        theory_all_sha256=sha(root/'THEORY_ALL.json'),source_hashes_sha256=sha(root/'SOURCE_HASHES.json'),
        native_tests=tests,validation_root=str(prior),validation_sha256=sha(prior/'RESULTS.json'),
        limitations='Full837 independent components/key union checked; no primitive existence, all-degree census or geometry/literature classification.')
    save(root/'RESULTS.json',result);print('RESULTS='+str(root/'RESULTS.json'),flush=True);lock.close()

if __name__=='__main__':main()
