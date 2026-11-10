"""Bind ramified-lift discovery/proofs and final six-slot launch inputs."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import csv
import fcntl
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

from audit_guided_native import audit_guided_native,audit_component_native
from profile_mechanism_breaker import install_bases,ramified_double_profile
from verify_smooth_geometry import parse_tuple,compose,inverse,cycles
from run_s18_linked_roots import state_key
from run_smooth_degree9 import sha,rows

PROJECT=Path(__file__).resolve().parents[2]

def pure(q):return tuple(2*q[x//2]+x%2 for x in range(2*len(q)))

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--prior',type=Path,required=True)
    p.add_argument('--inputs',type=Path,required=True);p.add_argument('--catalog',type=Path,required=True)
    p.add_argument('--binary',type=Path,required=True);p.add_argument('--auditor',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);args=p.parse_args();args.out.mkdir()
    root=args.out.resolve();prior=args.prior.resolve();inputs=args.inputs.resolve()
    if not Path('data').is_dir() or shutil.disk_usage(root).free<30*2**30:raise ValueError('T7 guard')
    lock=(PROJECT/'.run-control/.campaign.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    binary=args.binary.resolve();auditor=args.auditor.resolve();old=json.loads((prior/'RESULTS.json').read_text())
    env=os.environ.copy();env.update(DESIGN_BINARY=str(binary),GUIDED_AUDIT_BINARY=str(auditor))
    logs=[]
    for name in ['test_mechanism_breaker','test_native_components','test_theory_frontier','test_discovery_design','test_guided_native','test_forest_fast','test_forest_launchd','test_forest_recovery']:
        result=subprocess.run([sys.executable,str(_scripts_path.find(f'{name}.py'))],env=env,capture_output=True,text=True)
        path=root/f'{name}.log';path.write_text(result.stdout+result.stderr)
        if result.returncode:print(result.stdout+result.stderr,flush=True);raise ValueError('final regression: '+name)
        logs.append(dict(name=name,sha256=sha(path)));print('PASS',name,flush=True)
    text=subprocess.check_output([str(binary),'--check-controls',str(prior/'controls.tsv')],text=True)
    (root/'CONTROLS.log').write_text(text);print(text.strip(),flush=True)
    existing=json.loads(Path('data/theory_frontier_native_validation_20261013_v2/THEORY_ALL.json').read_text())['components']
    for e in existing:
        report=audit_component_native(Path(e['raw_folder']),auditor)
        if not report['independent_six_move_closure'] or report['table_sha256']!=e['table_sha256'] or len(report['state_keys'])!=e['states']:raise ValueError('final old catalog audit')
    print('FINAL_EXISTING_CLOSURES',len(existing),sum(e['states'] for e in existing),flush=True)
    data=json.loads(args.catalog.read_text());new=data['components']
    if len(new)!=4 or sum(e['states'] for e in new)!=5184:raise ValueError('ramified catalog mass')
    for e in new:
        if sha(Path(e['raw_folder'])/'STATES.tsv')!=e['table_sha256'] or not e['all_rows_and_six_quotient_edges_checked'] or e['kernel_order']!=128 or not e['kernel_code']['equal_parities_on_three_coarse_blocks']:
            raise ValueError('ramified full projection/kernel proof')
    finite=Path('data/ramified_double_finite_20261013')
    check=audit_guided_native(finite,auditor,collect_records=True,native_components=True);check.pop('_closed_records')
    if int(check['construction']['constructed'])!=144 or int(check['construction']['closed'])!=4 or int(check['construction']['unresolved']):raise ValueError('finite lift census')
    if {json.dumps(e['key']) for e in check['closed_identities']}!={json.dumps(e['component_identity_key']) for e in new}:raise ValueError('finite lifts do not match mutation discoveries')
    union=set()
    for e in new:
        with (Path(e['raw_folder'])/'STATES.tsv').open() as stream:
            keys={state_key(parse_tuple(row['tuple'],18)) for row in csv.DictReader(stream,delimiter='\t')}
        if len(keys)!=1296 or union&keys:raise ValueError('ramified component disjointness')
        union.update(keys)
    base=next(parse_tuple(line.split('\t',1)[1],9) for line in (inputs/'double_bases.tsv').read_text().splitlines() if line.startswith('9\t'))
    a,b=pure(base[0]),pure(base[1]);direct=set()
    for mask in range(512):
        c=tuple(2*base[2][x//2]+((x%2)^((mask>>(x//2))&1)) for x in range(18));d=inverse(compose(compose(a,b),c))
        if 18-len(cycles(c))!=13 or 18-len(cycles(d))!=13:continue
        t=(a,b,c,d);key=state_key(t)
        if key not in union:raise ValueError('unresolved direct double voltage')
        direct.add(key)
    if len(direct)!=144:raise ValueError('independent direct voltage cardinality')
    finite_report=dict(passed=True,all512_direct_binary_voltages_checked=True,qualifying144=len(direct),
        components4=len(new),states5184=len(union),every_qualifying_voltage_in_complete_component=True,
        independent_six_move_audit=check,source_catalog_sha256=sha(args.catalog))
    (root/'DOUBLE_EXACT.json').write_text(json.dumps(finite_report,indent=2)+'\n')
    # A small prefix of the genuinely new degree24 lift construction; caps are
    # retained for the long run, never relabeled as positive counts.
    pilot=root/'double24_prefix'
    subprocess.run([str(binary),'--d','24','--lane','double_lift','--double-bases',str(inputs/'double_bases.tsv'),
        '--frame-reuse','1','--seconds','5','--attempts','4','--cap','3000','--closure-seconds','1','--out',str(pilot)],check=True,capture_output=True)
    pilot_check=audit_guided_native(pilot,auditor,native_components=True)
    (pilot/'CHECKS.json').write_text(json.dumps(pilot_check,indent=2)+'\n');(pilot/'VERIFIED').write_text('Independent bounded prefix\n')
    for name in ['controls.tsv','primitive_exclusions.json']:shutil.copy2(prior/name,root/name)
    for name in ['mutation_seeds.tsv','double_bases.tsv']:shutil.copy2(inputs/name,root/name)
    shutil.copy2(inputs/'known_states.tsv',root/'known_states.tsv')
    with (root/'known_states.tsv').open('a') as stream:
        for e in new:
            for row in rows(Path(e['raw_folder'])/'STATES.tsv'):stream.write('18\t'+row['tuple']+'\n')
    queue=json.loads((inputs/'profile_queue.json').read_text());queue['components']+=new
    (root/'profile_queue.json').write_text(json.dumps(queue,indent=2)+'\n')
    source_hashes={str(p.relative_to(PROJECT)):sha(p) for p in [*sorted((PROJECT/'src').rglob('*.rs')),*sorted((PROJECT/'scripts').rglob('*.py'))]}
    (root/'SOURCE_HASHES.json').write_text(json.dumps(source_hashes,indent=2)+'\n')
    result=dict(old,passed=True,mechanism_passed=True,design_passed=True,native_components_passed=True,
        binary_sha256=sha(binary),audit_binary_sha256=sha(auditor),
        known_bank_sha256=sha(root/'known_states.tsv'),profile_queue_sha256=sha(root/'profile_queue.json'),
        mutation_seeds_sha256=sha(root/'mutation_seeds.tsv'),double_bases_sha256=sha(root/'double_bases.tsv'),
        primitive_exclusions_sha256=sha(root/'primitive_exclusions.json'),
        prior_validation_sha256=sha(prior/'RESULTS.json'),final_sources_sha256=sha(root/'SOURCE_HASHES.json'),
        tests=logs,ramified_catalog_sha256=sha(args.catalog),new_ramified_components=4,new_ramified_states=5184,
        double_exact_sha256=sha(root/'DOUBLE_EXACT.json'),double24_prefix=pilot_check['construction'],
        profile_components=len(queue['components']),known_states=960524,
        limitations='New ramified-lift mechanism in this project, still imprimitive. Read actual degree24 closed/unresolved counters; its structural/group-order profile remains pending. Geometry unclassified. No weak/primitive/literature novelty claim.')
    (root/'RESULTS.json').write_text(json.dumps(result,indent=2)+'\n');print('FINAL='+str(root/'RESULTS.json'),flush=True)

if __name__=='__main__':main()
