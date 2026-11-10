"""Bound new constructors/runner/keys to finite regressions and preserved proofs.

One bounded validation process, no unbounded activation. Original proof sources
are retained; user-owned source/report changes are not rewritten or rebound.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import fcntl
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

from audit_guided_native import audit_guided_native,audit_component_native,prime_known_cache,known_keys
from run_smooth_degree9 import sha

PROJECT=Path(__file__).resolve().parents[2]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',type=Path,required=True);p.add_argument('--auditor',type=Path,required=True)
    p.add_argument('--inputs',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    args=p.parse_args();root=args.out.resolve();root.mkdir()
    lock=(PROJECT/'.run-control/.campaign.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    binary=args.binary.resolve();auditor=args.auditor.resolve();inputs=args.inputs.resolve()
    source=root/'source';source.mkdir();hashes={}
    for path in [*sorted((PROJECT/'src').rglob('*.rs')),*sorted((PROJECT/'scripts').rglob('*.py'))]:
        rel=path.relative_to(PROJECT);target=source/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(path,target);hashes[str(rel)]=sha(target)
    (root/'SOURCE_HASHES.json').write_text(json.dumps(hashes,indent=2)+'\n')
    shutil.copy2(binary,root/'smooth_forest_fast');shutil.copy2(auditor,root/'smooth_negative_audit')
    env=os.environ.copy();env.update(DESIGN_BINARY=str(binary),GUIDED_AUDIT_BINARY=str(auditor))
    tests=[]
    for name in ['test_mechanism_breaker','test_native_components','test_theory_frontier','test_discovery_design',
                 'test_guided_native','test_forest_fast','test_forest_launchd','test_forest_recovery']:
        result=subprocess.run([sys.executable,str(_scripts_path.find(f'{name}.py'))],env=env,capture_output=True,text=True)
        path=root/f'{name}.log';path.write_text(result.stdout+result.stderr);tests.append(dict(name=name,sha256=sha(path)))
        if result.returncode:print(result.stdout+result.stderr,flush=True);raise RuntimeError('regression failed: '+name)
        print('PASS',name,flush=True)
    text=subprocess.check_output([str(binary),'--check-controls',str(inputs/'controls.tsv')],text=True)
    (root/'CONTROLS.log').write_text(text);print(text.strip(),flush=True)
    theory_path=Path('data/theory_frontier_native_validation_20261013_v2/THEORY_ALL.json')
    entries=json.loads(theory_path.read_text())['components'];mass=0;began=time.monotonic()
    for i,e in enumerate(entries):
        result=audit_component_native(Path(e['raw_folder']),auditor)
        if result['table_sha256']!=e['table_sha256'] or not result['independent_six_move_closure'] or len(result['state_keys'])!=e['states']:raise ValueError('existing catalog component changed')
        mass+=e['states']
        if i%100==0:print('EXISTING_CLOSURES',i+1,flush=True)
    if len(entries)!=837 or mass!=166056:raise ValueError('prior catalog mass')
    catalog_seconds=time.monotonic()-began
    proofs={};proof_records=[]
    for name in ['psl11','m11']:
        folder=Path('data')/f'primitive_{name}_pilot_20261013'
        check=audit_guided_native(folder,auditor,native_components=True)
        (root/f'{name}_GUIDED.json').write_text(json.dumps(check,indent=2)+'\n')
        proof=folder/'EXACT_POPULATION.json';report=json.loads(proof.read_text())
        if not report['passed'] or not report['independent_complete_population'] or report['counts']['qualifying']!=report['counts']['literal_negative']:raise ValueError('primitive finite proof')
        if any(sha(folder/name)!=digest for name,digest in report['sources'].items()):raise ValueError('primitive proof input changed')
        proofs[f'{name}_11']=dict(path=str(proof),sha256=sha(proof))
        proof_records.append(dict(group=name,proof_sha256=sha(proof),counts=report['counts']))
    (root/'primitive_exclusions.json').write_text(json.dumps(proofs,indent=2)+'\n')
    # A new constructor workload, not a repeated old search population.
    pilots=[]
    for n in [12,15,18]:
        folder=root/f'mutation_{n}';began=time.monotonic()
        subprocess.run([str(binary),'--d',str(n),'--lane','mutate','--mutation-seeds',str(inputs/'mutation_seeds.tsv'),
            '--frame-reuse','1','--seed',str(2026101300+n),'--seconds','5','--attempts','10000',
            '--known-bank',str(inputs/'known_states.tsv'),'--cap','3000','--closure-seconds','1','--out',str(folder)],check=True,capture_output=True)
        check=audit_guided_native(folder,auditor,collect_records=True,native_components=True);check.pop('_closed_records')
        (folder/'CHECKS.json').write_text(json.dumps(check,indent=2)+'\n');(folder/'VERIFIED').write_text('Independent bounded pilot\n')
        pilots.append(dict(degree=n,counts=check['construction'],seconds=time.monotonic()-began))
        print('MUTATION_PILOT',json.dumps(pilots[-1]),flush=True)
    # Validate exact bulk keys against the previous independently checked key
    # population; this is a derived cache, not a new smoothness certificate.
    initial=Path('data/theory_frontier_final_checks_20261013/known_states.tsv')
    expected=known_keys(initial);cache_report=prime_known_cache(initial,auditor)
    if known_keys(initial)!=expected or cache_report['records']!=166056:raise ValueError('native bulk keys differ')
    for name in ['known_states.tsv','controls.tsv','mutation_seeds.tsv','profile_queue.json','INPUTS.json']:
        shutil.copy2(inputs/name,root/name)
    projection=Path('data/marked_projection_pilot_20261013/PROFILES.json')
    result=dict(passed=True,design_passed=True,mechanism_passed=True,native_components_passed=True,
        controls=66,control_states=25192,binary_sha256=sha(binary),audit_binary_sha256=sha(auditor),
        known_bank_sha256=sha(root/'known_states.tsv'),mutation_seeds_sha256=sha(root/'mutation_seeds.tsv'),
        profile_queue_sha256=sha(root/'profile_queue.json'),primitive_exclusions_sha256=sha(root/'primitive_exclusions.json'),
        input_ledger_sha256=sha(root/'INPUTS.json'),tests=tests,existing_components=837,existing_states=mass,
        existing_full_native_audit_seconds=catalog_seconds,primitive_exact_proofs=proof_records,
        mutation_pilots=pilots,bulk_known_keys=cache_report,marked_projection_pilot_sha256=sha(projection),
        sources_sha256=sha(root/'SOURCE_HASHES.json'),
        limitations='Finite named group exclusions; sampled mutations; no new primitive certificate or geometry claim. Profiling queue not finished.')
    if any(sha(PROJECT/name)!=digest for name,digest in hashes.items()):raise ValueError('source changed during validation')
    (root/'RESULTS.json').write_text(json.dumps(result,indent=2)+'\n');print('VALIDATION='+str(root/'RESULTS.json'),flush=True)

if __name__=='__main__':main()
