"""Check the final launcher-only adjustments and bind preserved pilot leads."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

from run_smooth_degree9 import sha

PROJECT=Path(__file__).resolve().parents[2]

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--prior',type=Path,required=True)
    p.add_argument('--binary',type=Path,required=True);p.add_argument('--auditor',type=Path,required=True)
    args=p.parse_args();root=args.prior.resolve();old=json.loads((root/'RESULTS.json').read_text())
    if sha(args.binary)!=old['binary_sha256'] or sha(args.auditor)!=old['audit_binary_sha256']:raise ValueError('kernel changed after final checks')
    env=os.environ.copy();env.update(DESIGN_BINARY=str(args.binary.resolve()),GUIDED_AUDIT_BINARY=str(args.auditor.resolve()))
    logs=[]
    for name in ['test_mechanism_breaker','test_forest_fast','test_forest_launchd','test_forest_recovery']:
        result=subprocess.run([sys.executable,str(_scripts_path.find(f'{name}.py'))],env=env,capture_output=True,text=True)
        path=root/f'launch_{name}.log';path.write_text(result.stdout+result.stderr)
        if result.returncode:print(result.stdout+result.stderr);raise ValueError('launch regression')
        logs.append(dict(name=name,sha256=sha(path)));print('PASS',name,flush=True)
    pilot=root/'double24_prefix'
    if not (pilot/'VERIFIED').exists():raise ValueError('unverified double24 saved starts')
    check=json.loads((pilot/'CHECKS.json').read_text())
    if any(sha(pilot/name)!=digest for name,digest in check['sources'].items()):raise ValueError('saved double24 journal changed')
    with (pilot/'RESOLUTIONS.tsv').open() as stream:
        leads=[dict(degree=24,seed=r['seed'],source=str(pilot),candidate=r['candidate'],checks_sha256=sha(pilot/'CHECKS.json'))
               for r in csv.DictReader(stream,delimiter='\t') if r['status'] in ['capped','timeout']]
    (root/'saved_leads.json').write_text(json.dumps(leads,indent=2)+'\n')
    bindings={str(p.relative_to(PROJECT)):sha(p) for p in [*sorted((PROJECT/'src').rglob('*.rs')),*sorted((PROJECT/'scripts').rglob('*.py'))]}
    (root/'LAUNCH_SOURCES.json').write_text(json.dumps(bindings,indent=2)+'\n')
    result=dict(old,launcher_only_checks=logs,launcher_changes='sanitize underscore lane names for interrupted degree recovery; keep two discovery slots alongside four replay/audit slots',
        parent_final_results_sha256=sha(root/'RESULTS.json'),final_launch_sources_sha256=sha(root/'LAUNCH_SOURCES.json'),
        saved_leads_sha256=sha(root/'saved_leads.json'),saved_leads=len(leads))
    (root/'LAUNCH_VALIDATION.json').write_text(json.dumps(result,indent=2)+'\n');print('LAUNCH_CHECKS='+str(root/'LAUNCH_VALIDATION.json'),flush=True)

if __name__=='__main__':main()
