"""Final input-guard/control rebinding; no repeat of the discovery population.

Replaying degree18's measured seed is a known-membership regression: it must
create NO new component after importing the enlarged certified-state bank.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import json
from pathlib import Path
import subprocess
from audit_guided_native import audit_guided_native
from run_smooth_degree9 import sha

PROJECT=Path(__file__).resolve().parents[2]

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--catalog',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--binary',type=Path,required=True);p.add_argument('--auditor',type=Path,required=True)
    args=p.parse_args();root=args.out.resolve();root.mkdir();catalog=args.catalog.resolve()
    binary=args.binary.resolve();auditor=args.auditor.resolve();logs=[]
    commands=[['python3','scripts/smooth/test_discovery_design.py'],['python3','scripts/smooth/test_guided_native.py'],
              ['python3','scripts/smooth/test_forest_fast.py'],['python3','scripts/smooth/test_forest_launchd.py'],
              ['python3','scripts/smooth/test_forest_recovery.py'],
              [str(binary),'--check-controls',str(catalog/'controls.tsv')]]
    import os
    env=dict(os.environ,DESIGN_BINARY=str(binary),GUIDED_AUDIT_BINARY=str(auditor))
    for i,command in enumerate(commands):
        result=subprocess.run(command,cwd=PROJECT,env=env,capture_output=True,text=True,check=True)
        path=root/f'check_{i}.log';path.write_text(result.stdout+result.stderr)
        logs.append(dict(command=command,log=str(path),sha256=sha(path)))
        print('PASS',command,flush=True)
    folder=root/'expanded_known_regression'
    subprocess.run(list(map(str,[binary,'--d','18','--lane','solved','--frame-reuse','1','--seed','202610148',
                   '--seconds','20','--attempts','300','--known-bank',catalog/'known_states.tsv',
                   '--compact-frames','--export-negative-cache','--out',folder])),check=True,stdout=subprocess.DEVNULL)
    check=audit_guided_native(folder,auditor)
    if check['closed_identities'] or int(check['construction']['unresolved']) or check['known_aliases']==0:
        raise ValueError('expanded bank failed to retire known components')
    (folder/'CHECKS.json').write_text(json.dumps(check,indent=2)+'\n')
    for name in ['known_states.tsv','controls.tsv']:(root/name).write_bytes((catalog/name).read_bytes())
    prior=json.loads((catalog/'RESULTS.json').read_text())
    result=dict(prior,binary_sha256=sha(binary),audit_binary_sha256=sha(auditor),
                final_controls=logs,expanded_bank_regression=check['construction'],
                final_known_aliases=check['known_aliases'],
                measured_binary_sha256=prior['binary_sha256'],
                final_guard_change='reject non-product-one known-bank tuples; valid proposal kernels unchanged',
                catalog_sha256=sha(catalog/'CATALOG.json'))
    (root/'RESULTS.json').write_text(json.dumps(result,indent=2)+'\n')
    print('FINAL_CHECKS='+str(root/'RESULTS.json'),flush=True)

if __name__=='__main__':main()
