#!/usr/bin/env python3
"""Seal C213 validation and independently replay specific cyclic witnesses.

Uses one reserved CPU slot. Does not launch/stop a campaign or publish data.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
import sys

from audit_guided_native import audit_guided_native
from run_genus_three import sha, write

PROJECT=Path(__file__).resolve().parents[2]
SOURCES=('src/smooth_cyclic_screen.rs','src/smooth_genus_three.rs',
    'src/bin/smooth_negative_audit.rs','scripts/genus3/run_genus_three.py',
    'scripts/genus3/test_genus_three.py','scripts/genus3/verify_genus_three_boundary.py',
    'scripts/genus3/test_genus_three_boundary.py','scripts/genus3/watch_genus_three.py',
    'scripts/genus3/test_watch_genus_three.py','scripts/genus3/validate_genus_three_math.py',
    'scripts/genus3/validate_genus_three_cyclic.py','scripts/genus3/benchmark_genus_three_cyclic.py',
    'scripts/genus3/record_genus_three_cyclic.py')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('build',type=Path)
    p.add_argument('early',type=Path)
    p.add_argument('late',type=Path)
    p.add_argument('out',type=Path)
    a=p.parse_args();a.out.mkdir()
    validation=json.loads((a.build/'VALIDATION.json').read_text())
    if not validation['passed']:raise ValueError('incomplete source validation')
    for name in ('wide','compact','auditor'):
        key={'wide':'binary_sha256','compact':'compact_binary_sha256','auditor':'audit_binary_sha256'}[name]
        if sha(a.build/'bin'/name)!=validation[key]:raise ValueError('binary changed')
    for rel in SOURCES[:3]:
        if sha(PROJECT/rel)!=validation['compiled_source_files'][rel]:raise ValueError('scientific source changed after validation')
    benchmarks={}
    for name,root in [('early',a.early),('late',a.late)]:
        result=json.loads((root/'CHECKS.json').read_text())
        if not result['passed'] or result['binary_sha256']!=validation['binary_sha256']:
            raise ValueError('benchmark does not bind final binary')
        benchmarks[name]=dict(path=str(root.resolve()),sha256=sha(root/'CHECKS.json'),summaries=result['summaries'])
    env=dict(os.environ,G3_CYCLIC_VALIDATION=str(a.build.resolve()))
    tests=[]
    for filename in ('test_genus_three.py','test_genus_three_boundary.py','test_watch_genus_three.py'):
        cmd=[sys.executable,'-m','unittest','discover','-s',str(PROJECT/'scripts/genus3'),'-p',filename,'-v']
        run=subprocess.run(cmd,check=True,capture_output=True,text=True,env=env)
        if 'skipped' in run.stderr:raise ValueError('required native test skipped')
        tests.append(dict(command=cmd,stdout=run.stdout,stderr=run.stderr,passed=True))
    journals=[]
    for n in (12,16):
        folder=a.build/f'cyclic_witnesses_{n}'
        cmd=[str(a.build/'bin/compact'),'--d',str(n),'--lane','collision','--seed','21309',
            '--seconds','60','--attempts','1000000','--cyclic-before-deep','--trace-cyclic','1000',
            '--closure-seconds','15','--cap','1000000','--out',str(folder)]
        subprocess.run(cmd,check=True,capture_output=True,text=True)
        checked=audit_guided_native(folder,a.build/'bin/auditor',native_components=True)
        with (folder/'CYCLIC.tsv').open() as stream:stats=next(csv.DictReader(stream,delimiter='\t'))
        if int(stats['persisted'])<100:raise ValueError('insufficient specific cyclic regression witnesses')
        write(folder/'CHECKS.json',checked)
        journals.append(dict(command=cmd,path=str(folder.resolve()),cyclic=stats,native=checked['independent_native'],
            checks_sha256=sha(folder/'CHECKS.json')))
    proof_bindings={rel:sha(PROJECT/rel) for rel in (
        'evidence/genus_three_pure_math_20261022/MATHEMATICS.md',
        'evidence/genus_three_star_obstruction_20261022/MATHEMATICS.md',
        'evidence/genus_three_dihedral_obstruction_20261022/MATHEMATICS.md')}
    receipt=dict(passed=True,claim='C213',source_genus=3,goal_found=False,
        validation=dict(path=str((a.build/'VALIDATION.json').resolve()),sha256=sha(a.build/'VALIDATION.json'),
            rust_tests=validation['rust_tests'],cyclic_rust_tests=4,
            controls=len(validation['old_controls'])+len(validation['high_controls']),
            control_states=sum(int(x['summary']['states']) for x in validation['old_controls']+validation['high_controls'])),
        benchmarks=benchmarks,tests=tests,cyclic_journals=journals,
        code_sha256={rel:sha(PROJECT/rel) for rel in SOURCES},mathematical_inputs=proof_bindings,
        finite_population=dict(degree=16,indices=[7,7,7,15],frames=['8+8','7+9'],
            matchings=[64*5040,36*5040],total=504000,
            convention='labelled matching candidates in normalized frames, NOT504000 distinct quotient states',
            completeness='exact coverage per persisted cursor interval; unresolved closures retained'),
        policy=dict(cyclic='late before closure; early mode diagnostic, not an assumed speed gain',
            geometry='no new geometric exclusion or public relabeling',gpu=False,workers_total_max=6),
        development_failures=['v1 test incorrectly assumed all genus0 tuples are all-smooth; corrected to genus0 with a full-cycle entry. Failed snapshot retained.',
            'Python interval test incorrectly required a survivor in a tiny interval; expanded intervals, constructor unchanged.'],
        limitations=['No global genus3 nonexistence proof','No new nonisotrivial example established by a faster filter',
            'Kernel measurements do not estimate time-to-discovery','Human mathematical proof review not recorded'])
    write(a.out/'CHECKS.json',receipt)
    release=dict(validation,cyclic_release_receipt=dict(path=str((a.out/'CHECKS.json').resolve()),sha256=sha(a.out/'CHECKS.json')))
    write(a.out/'VALIDATION.json',release)
    print(json.dumps(dict(passed=True,controls=receipt['validation']['controls'],states=receipt['validation']['control_states'],
        specific_cyclic_witnesses=sum(int(x['cyclic']['persisted']) for x in journals))))


if __name__=='__main__':main()
