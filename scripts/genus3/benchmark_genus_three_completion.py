#!/usr/bin/env python3
"""One bounded CPU slot. No launch, geometry classification, or exclusions.

Freeze the admitted clean source, overlay only the two constructor sources,
then compare late labelled leads/sec under the same literal word bank. A
different proposal law is never described as a whole-search speed multiplier.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import csv
import json
from pathlib import Path
import random
import shutil
import subprocess
import time
from audit_guided_native import audit_guided_native
from run_genus_three import sha,write

PROJECT=Path(__file__).resolve().parents[2]

def run(argv,logs,name):
    start=time.perf_counter()
    with (logs/(name+'.out')).open('w') as out,(logs/(name+'.err')).open('w') as err:
        subprocess.run(list(map(str,argv)),stdout=out,stderr=err,check=True)
    return time.perf_counter()-start

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('base',type=Path);p.add_argument('out',type=Path)
    p.add_argument('--seconds',type=int,default=2);p.add_argument('--repeats',type=int,default=3)
    p.add_argument('--degrees',nargs='+',type=int,default=[12,15,18,21,24])
    a=p.parse_args()
    if a.seconds<1 or a.repeats<1 or any(not 11<=d<=25 for d in a.degrees):p.error('invalid bounded pilot')
    a.out.mkdir();source=a.out/'source';logs=a.out/'logs';logs.mkdir();(a.out/'bin').mkdir()
    old=json.loads((a.base/'VALIDATION.json').read_text())
    for rel,digest in old['compiled_source_files'].items():
        src=a.base/'source'/rel
        if sha(src)!=digest:raise ValueError('base source changed: '+rel)
        dst=source/rel;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst)
    for rel in ['src/smooth_genus_three.rs','src/smooth_joint_completion.rs']:
        shutil.copyfile(PROJECT/rel,source/rel)
    bank=Path(json.loads((a.base/'CHECKPOINT.json').read_text())['word_bank'])
    shutil.copyfile(bank,a.out/'WORDS.txt')
    auditor=a.base/'smooth_negative_audit';assert sha(auditor)==old['audit_binary_sha256']
    commands=[]
    def build(name,args):
        argv=['rustc','--edition=2021','-O',*args];seconds=run(argv,logs,name)
        commands.append(dict(name=name,command=list(map(str,argv)),seconds=seconds))
    build('joint_tests',['--test',source/'src/smooth_joint_completion.rs','-o',a.out/'bin/joint_tests'])
    commands.append(dict(name='joint_tests',seconds=run([a.out/'bin/joint_tests','--test-threads','1'],logs,'joint_tests_run')))
    build('compact',[source/'src/bin/smooth_genus_three.rs','-o',a.out/'bin/compact'])
    binary=a.out/'bin/compact';rows=[];order=random.Random(2026102209)
    for degree in a.degrees:
        for rep in range(a.repeats):
            modes=['collision','greedy','csp'];order.shuffle(modes);block={}
            for mode in modes:
                folder=a.out/f'{degree}_{rep}_{mode}'
                argv=[str(binary),'--d',str(degree),'--lane','collision' if mode=='collision' else 'joint-genus3',
                      '--seed',str(20271022+rep),'--seconds',str(a.seconds),'--cap','1000000',
                      '--closure-seconds','10','--words',str(a.out/'WORDS.txt'),'--out',str(folder)]
                if mode=='csp':argv+=['--joint-csp','--joint-nodes','5000']
                wall=run(argv,logs,folder.name)
                checked=audit_guided_native(folder,auditor,native_components=True);write(folder/'CHECKS.json',checked)
                with (folder/'COUNTS.tsv').open() as f:counts=next(csv.DictReader(f,delimiter='\t'))
                with (folder/'STAGES.tsv').open() as f:stages=next(csv.DictReader(f,delimiter='\t'))
                leads={line.split('\t',1)[1] for line in (folder/'LEADS.tsv').read_text().splitlines()}
                joint=None
                if (folder/'JOINT.tsv').exists():
                    with (folder/'JOINT.tsv').open() as f:joint=next(csv.DictReader(f,delimiter='\t'))
                block[mode]=dict(command=argv,wall_seconds=wall,counts=counts,stages=stages,
                    distinct_labelled_leads=len(leads),late_leads_per_second=len(leads)/wall,joint=joint,
                    independent=checked['independent_native'],checks_sha256=sha(folder/'CHECKS.json'))
            rows.append(dict(degree=degree,replicate=rep,results=block))
            print(degree,rep,{k:round(v['late_leads_per_second'],3) for k,v in block.items()},flush=True)
            write(a.out/'PROGRESS.json',dict(blocks=rows))
    write(a.out/'RESULT.json',dict(passed=True,commands=commands,blocks=rows,
        compiled_source_files={str(q.relative_to(source)):sha(q)for q in sorted(source.rglob('*.rs'))},
        binary_sha256=sha(binary),audit_binary_sha256=sha(auditor),words_sha256=sha(a.out/'WORDS.txt'),
        scope='different biased proposals,one auxiliary CPU and five continuing search workers; labelled leads are not orbits; no speedup asserted for zero-v-zero'))

if __name__=='__main__':main()
