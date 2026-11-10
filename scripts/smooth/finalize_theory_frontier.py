"""Bind final publication guard and mixed-degree partial audits to launch inputs.

Preserves all earlier benchmarks and catalogs. No unbounded activation.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import json
from pathlib import Path
import shutil
import subprocess

from audit_guided_native import audit_guided_native
from run_smooth_degree9 import sha

PROJECT=Path(__file__).resolve().parents[2]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--prior',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--binary',type=Path,required=True);p.add_argument('--auditor',type=Path,required=True)
    args=p.parse_args();root=args.out.resolve();root.mkdir();prior=args.prior.resolve()
    binary=args.binary.resolve();auditor=args.auditor.resolve()
    env=__import__('os').environ.copy();env.update(DESIGN_BINARY=str(binary),GUIDED_AUDIT_BINARY=str(auditor))
    logs=[]
    for name in ['test_native_components','test_theory_frontier','test_discovery_design','test_guided_native','test_forest_fast','test_forest_launchd','test_forest_recovery']:
        result=subprocess.run([__import__('sys').executable,str(_scripts_path.find(f'{name}.py'))],env=env,check=True,capture_output=True,text=True)
        path=root/f'{name}.log';path.write_text(result.stdout+result.stderr);logs.append(dict(name=name,sha256=sha(path)))
    old=json.loads((prior/'RESULTS.json').read_text())
    if sha(binary)!=old['binary_sha256'] or sha(auditor)!=old['audit_binary_sha256']:raise ValueError('native executable changed')
    pilot=Path(old['validation_root'])/'pilot_24_quotient'
    check=audit_guided_native(pilot,auditor,collect_records=True,native_components=True)
    records=check.pop('_closed_records');baseline=json.loads((pilot/'CHECKS.json').read_text())
    for name in ['construction','components','closed_identities','known_aliases','input_banks','exported_negative_bank','sources']:
        if json.loads(json.dumps(check[name]))!=baseline[name]:raise ValueError('full positive/partial Python comparison: '+name)
    (root/'DEGREE24_NATIVE_CHECKS.json').write_text(json.dumps(check,indent=2)+'\n')
    for name in ['known_states.tsv','controls.tsv','saved_leads.json']:shutil.copy2(prior/name,root/name)
    source_hashes={str(p.relative_to(PROJECT)):sha(p) for p in [*sorted((PROJECT/'src').rglob('*.rs')),*sorted((PROJECT/'scripts').rglob('*.py'))]}
    (root/'SOURCE_HASHES.json').write_text(json.dumps(source_hashes,indent=2)+'\n')
    result=dict(old,final_publication_binding_passed=True,final_tests=logs,
        final_source_hashes_sha256=sha(root/'SOURCE_HASHES.json'),
        degree24_native_positive_partial_mass=check['construction'],
        final_checks_sha256=sha(root/'DEGREE24_NATIVE_CHECKS.json'),
        final_script_change='bind verified publication records to independent report; permutation/audit kernels and measured binaries unchanged',
        prior_native_validation_root=str(prior),prior_native_validation_sha256=sha(prior/'RESULTS.json'))
    (root/'RESULTS.json').write_text(json.dumps(result,indent=2)+'\n');print('FINAL_CHECKS='+str(root/'RESULTS.json'),flush=True)

if __name__=='__main__':main()
