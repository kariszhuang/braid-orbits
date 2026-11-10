#!/usr/bin/env python3
"""One-core frozen exact slices plus independent complete audits. Not a pool."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

from run_genus_three import sha,write


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('out',type=Path);parser.add_argument('auditor',type=Path)
    args=parser.parse_args();root=args.out.resolve();root.mkdir();project=Path(__file__).resolve().parents[2]
    names=['src/bin/genus_three_chords.rs','scripts/genus3/run_genus_three_chords.py','scripts/genus3/verify_genus_three_chords.py','scripts/genus3/test_genus_three_chords.py']
    # Freeze the verifier's Python dependencies too; the compiled native auditor
    # remains the previously admitted independent binary, not a new copy of ours.
    names+=[str(q.relative_to(project)) for q in sorted((project/'scripts').rglob('*.py')) if str(q.relative_to(project)) not in names]
    for name in names:
        dest=root/'source'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(project/name,dest)
    auditor=root/'independent_auditor';shutil.copyfile(args.auditor,auditor);auditor.chmod(0o700)
    command=['rustc','--edition=2021','-O','-C','codegen-units=1',str(root/'source/src/bin/genus_three_chords.rs'),'-o',str(root/'chords')]
    start=time.perf_counter();subprocess.run(command,check=True,capture_output=True,text=True)
    build_seconds=time.perf_counter()-start;results=[]
    for n in (7,9,11):
        for mode in ('transposition','double'):
            folder=root/f'{mode}_{n}';command=[str(root/'chords'),str(n),str(folder),mode]
            subprocess.run(command,check=True,capture_output=True,text=True)
            command=[sys.executable,str(root/'source/scripts/genus3/verify_genus_three_chords.py'),str(folder),str(auditor)]
            subprocess.run(command,check=True,capture_output=True,text=True)
            checked=json.loads((folder/'EXACT_AUDIT.json').read_text());results.append(dict(degree=n,mode=mode,coverage=checked['coverage'],maximum_word=checked['maximum_word'],independent_seconds=checked['independent_seconds'],audit_sha256=sha(folder/'EXACT_AUDIT.json')))
            print(n,mode,checked['status'],checked['coverage']['qualified'],flush=True)
    env=dict(os.environ,G3_CHORD_BINARY=str(root/'chords'),G3_CHORD_AUDITOR=str(auditor))
    tests=subprocess.run([sys.executable,'-m','unittest','discover','-s',str(root/'source/scripts/genus3'),'-p','test_genus_three_chords.py','-v'],check=True,capture_output=True,text=True,env=env)
    write(root/'TESTS.json',dict(passed=True,stdout=tests.stdout,stderr=tests.stderr))
    write(root/'RESULT.json',dict(passed=True,build_seconds=build_seconds,cases=results,binary_sha256=sha(root/'chords'),audit_binary_sha256=sha(auditor),tests_sha256=sha(root/'TESTS.json'),
        sources={str(q.relative_to(root/'source')):sha(q)for q in (root/'source').rglob('*')if q.is_file()and q.suffix in ('.rs','.py')},
        scope='degree7/9 controls,degree11 two finite full-product slices; combine with the separate scalar and single-cycle proofs before claiming a genus3 degree11 exclusion'))


if __name__=='__main__':main()
