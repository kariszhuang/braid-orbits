#!/usr/bin/env python3
"""Bounded one-core exact passport audit. No production pool or new theorem filter."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import time
from run_genus_three import sha,write
from verify_genus_three_triangle import audit

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('out',type=Path);p.add_argument('auditor',type=Path)
    a=p.parse_args();project=Path(__file__).resolve().parents[2];a.out.mkdir();src=a.out/'source';src.mkdir();events=[]
    for name in ['src/bin/genus_three_triangle.rs','scripts/genus3/verify_genus_three_triangle.py','scripts/genus3/run_genus_three_triangle.py']:
        dest=src/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(project/name,dest)
    cmd=['rustc','--edition=2021','-O','-C','codegen-units=1',str(src/'src/bin/genus_three_triangle.rs'),'-o',str(a.out/'triangle')]
    began=time.perf_counter();subprocess.run(cmd,check=True,capture_output=True,text=True);events.append(dict(command=cmd,seconds=time.perf_counter()-began))
    cases=[]
    for n in [5,9,13]:
        folder=a.out/f'degree_{n}';cmd=[str(a.out/'triangle'),str(n),str(folder)]
        began=time.perf_counter();subprocess.run(cmd,check=True,capture_output=True,text=True)
        produced=time.perf_counter()-began;began=time.perf_counter();checked=audit(folder,a.auditor)
        audited=time.perf_counter()-began;write(folder/'EXACT_AUDIT.json',checked)
        cases.append(dict(degree=n,producer_seconds=produced,audit_seconds=audited,
            coverage=checked['coverage'],audit_sha256=sha(folder/'EXACT_AUDIT.json')))
        print(n,checked['coverage'],flush=True)
    write(a.out/'RESULT.json',dict(passed=True,events=events,cases=cases,auditor_sha256=sha(a.auditor),binary_sha256=sha(a.out/'triangle'),
        source_sha256={str(q.relative_to(src)):sha(q)for q in src.rglob('*')if q.is_file()},
        scope='degrees5/9 independent regressions; degree13 one complete genus3 passport,not a whole degree; no geometry classifier change'))

if __name__=='__main__':main()
