"""Single-CPU standalone intrinsic-cover build from an already clean snapshot."""
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
import time

from run_genus_three import sha,write

PROJECT=Path(__file__).resolve().parents[2]


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('base',type=Path);p.add_argument('out',type=Path);a=p.parse_args()
    a.out.mkdir();source=a.out/'source';(a.out/'bin').mkdir();(a.out/'logs').mkdir()
    old=json.loads((a.base/'VALIDATION.json').read_text())
    for rel,digest in old['compiled_source_files'].items():
        src=a.base/'source'/rel
        if sha(src)!=digest:raise ValueError('base source changed '+rel)
        dest=source/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dest)
    overlays=['src/intrinsic_s3.rs','src/intrinsic_search.rs','src/bin/genus_three_intrinsic.rs','src/bin/smooth_negative_audit.rs']
    for rel in overlays:
        dest=source/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(PROJECT/rel,dest)
    # Compile-time include_str! certificate fixtures from the frozen base.
    for rel in('evidence/certificates/CERTIFIED_A25_weak.txt','evidence/certificates/CERTIFIED_A25_order3_weak.txt','evidence/spectral_census/s18_position_20260926/NEGATIVE_VECTOR.tsv'):
        src=a.base/'source'/rel;dest=source/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dest)
    events=[]
    for name,entry,flags in [('intrinsic','src/bin/genus_three_intrinsic.rs',['--cfg','smooth_wide_keys']),('auditor','src/bin/smooth_negative_audit.rs',[]),('coordinate-tests','src/intrinsic_s3.rs',['--test'])]:
        command=['rustc','--edition=2021','-O','-C','codegen-units=1',*flags,str(source/entry),'-o',str(a.out/'bin'/name)];start=time.monotonic()
        with(a.out/'logs'/f'{name}.out').open('w')as out,(a.out/'logs'/f'{name}.err').open('w')as err:
            subprocess.run(command,check=True,stdout=out,stderr=err,env=dict(os.environ,CARGO_MANIFEST_DIR=str(source)),timeout=180)
        events.append(dict(command=command,seconds=time.monotonic()-start))
    tests=subprocess.run([str(a.out/'bin/coordinate-tests'),'--test-threads','1'],check=True,capture_output=True,text=True,timeout=30)
    write(a.out/'BUILD.json',dict(passed=True,base=str(a.base.resolve()),overlays=overlays,events=events,coordinate_tests=tests.stdout,
        sources={str(p.relative_to(source)):sha(p)for p in source.rglob('*.rs')},binaries={name:sha(a.out/'bin'/name)for name in('intrinsic','auditor','coordinate-tests')},
        preserved_clean_main_sha256=sha(source/'src/main.rs'),qualification='constructor and auditor compiled independently; no scientific closure conclusion from compilation'))
    print(json.dumps(dict(passed=True,out=str(a.out))))


if __name__=='__main__':main()
