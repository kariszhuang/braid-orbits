"""Bounded exact regression receipts, never starts or stops production."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess
from audit_guided_native import audit_component_native
from genus_three_high_controls import controls
from run_genus_three import sha,write


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['binary','auditor','tests','out','controls','source']:
        p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--compact-binary',type=Path)
    a=p.parse_args();a.out.mkdir();tests={}
    for group in ['g3_partition_tests','return_solver_tests','genus_three_tests','forest_fast_tests',
                  'structural_novelty_tests','discovery_design_tests','smooth_tests']:
        run=subprocess.run([str(a.tests),group,'--nocapture','--test-threads','1'],capture_output=True,text=True)
        (a.out/(group+'.txt')).write_text(run.stdout+run.stderr)
        if run.returncode:raise ValueError(f'Rust test group failed: {group}')
        tests[group]=run.stdout.split('test result:')[1].strip()
    subprocess.run([str(a.binary),'--controls',str(a.controls),'--verifier','macro','--compare','--check-cyclic',
                    '--out',str(a.out/'old_controls')],check=True,capture_output=True,text=True)
    old=[]
    for folder in sorted((a.out/'old_controls').glob('control_*')):
        v=audit_component_native(folder,a.auditor);v.pop('state_keys');old.append(v)
    if a.compact_binary:
        subprocess.run([str(a.compact_binary),'--controls',str(a.controls),'--verifier','macro','--compare','--check-cyclic',
                        '--out',str(a.out/'compact_controls')],check=True,capture_output=True,text=True)
        for folder in sorted((a.out/'compact_controls').glob('control_*')):
            previous=audit_component_native(a.out/'old_controls'/folder.name,a.auditor)
            compact=audit_component_native(folder,a.auditor)
            if previous['state_keys']!=compact['state_keys']:raise ValueError('compact/wide state-set mismatch')
    high=[]
    for case in controls():
        n=case['degree'];folder=a.out/f'high_{n}'
        subprocess.run([str(a.binary),'--d',str(n),'--tuple',case['text'],'--out',str(folder),
                        '--closure-seconds','60','--cap','10000'],check=True,capture_output=True,text=True)
        v=audit_component_native(folder,a.auditor);v.pop('state_keys')
        if not v['independent_six_move_closure']:raise ValueError('positive control did not close')
        # The independent auditor uses literal arrays and all six moves; compare
        # the ordinary producer's entire canonical state set as another check.
        summary=v['summary'];control=a.out/f'high_input_{n}.tsv'
        control.write_text(f"{n}\t{summary['states']}\t{case['text']}\n")
        subprocess.run([str(a.binary),'--controls',str(control),'--verifier','macro','--compare','--check-cyclic'],
                        check=True,capture_output=True,text=True)
        if n>34 and summary['group_order']!='unknown':raise ValueError('large order overflow guard')
        high.append(v)
    corrupt=[]
    for kind in ['out_of_range_label','invalid_edge']:
        folder=a.out/('corrupt_'+kind);shutil.copytree(a.out/'high_63',folder)
        rows=(folder/'STATES.tsv').read_text().splitlines();fields=rows[1].split('\t')
        if kind=='out_of_range_label':
            fields[-1]=re.sub(r'(?<!\d)63(?!\d)','64',fields[-1],count=1)
        else:
            edges=fields[-2].split(',');edges[0]='999999';fields[-2]=','.join(edges)
        rows[1]='\t'.join(fields);(folder/'STATES.tsv').write_text('\n'.join(rows)+'\n')
        failed=subprocess.run([str(a.auditor),'--component',str(folder),'63',str(a.out/(kind+'.tsv'))],capture_output=True,text=True)
        if failed.returncode==0:raise ValueError('corrupt independent certificate was accepted')
        corrupt.append(dict(kind=kind,returncode=failed.returncode,error=failed.stderr.strip()))
    result=dict(passed=True,binary_sha256=sha(a.binary),audit_binary_sha256=sha(a.auditor),
        rust_tests=tests,old_controls=old,high_controls=high,adversarial_certificate_checks=corrupt,
        compiled_source_files={str(p.relative_to(a.source)):sha(p)for p in sorted(a.source.rglob('*.rs'))},
        compiler_variants=dict(wide='rustc --edition=2021 --cfg smooth_wide_keys -O',compact='rustc --edition=2021 -O',auditor='rustc --edition=2021 -O; standalone MAX=63'),
        qualification='positive controls are established isotrivial examples, not discoveries; ordinary vs macro state sets and independent six-edge closure')
    if a.compact_binary:result['compact_binary_sha256']=sha(a.compact_binary)
    write(a.out/'VALIDATION.json',result)
    print(json.dumps(dict(passed=True,old_controls=len(old),old_states=sum(int(v['summary']['states'])for v in old),
        high_controls=len(high),high_states=sum(int(v['summary']['states'])for v in high)),indent=2))


if __name__=='__main__':main()
