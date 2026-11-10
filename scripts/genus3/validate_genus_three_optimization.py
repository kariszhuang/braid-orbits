#!/usr/bin/env python3
"""Freeze/build/test on one CPU while five search workers continue. No launch."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import csv
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import time
from audit_guided_native import audit_guided_native
from run_genus_three import sha,write
from verify_genus_three_exclusions import audit_exclusions

PROJECT=Path(__file__).resolve().parents[2]
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('base',type=Path);p.add_argument('review',type=Path);p.add_argument('out',type=Path)
    a=p.parse_args();a.out.mkdir();source=a.out/'source';logs=a.out/'logs';logs.mkdir();(a.out/'bin').mkdir();events=[]
    old=json.loads((a.base/'VALIDATION.json').read_text())
    for rel,h in old['compiled_source_files'].items():
        src=a.base/'source'/rel;assert sha(src)==h
        dst=source/rel;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst)
    for rel in ['src/smooth_genus_three.rs','src/smooth_joint_completion.rs']:
        shutil.copyfile(PROJECT/rel,source/rel)
    for rel in ['evidence/certificates/CERTIFIED_A25_weak.txt','evidence/certificates/CERTIFIED_A25_order3_weak.txt',
                'evidence/spectral_census/s18_position_20260926/NEGATIVE_VECTOR.tsv']:
        dest=source/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(PROJECT/rel,dest)
    def run(name,cmd,env=None):
        start=time.perf_counter()
        with (logs/(name+'.out')).open('w')as out,(logs/(name+'.err')).open('w')as err:
            subprocess.run(list(map(str,cmd)),check=True,stdout=out,stderr=err,env=env or dict(os.environ,CARGO_MANIFEST_DIR=str(source)))
        elapsed=time.perf_counter()-start;events.append(dict(name=name,command=list(map(str,cmd)),seconds=elapsed));return elapsed
    for name,flags in [('compact',[]),('wide',['--cfg','smooth_wide_keys']),('tests',['--cfg','smooth_wide_keys','--test'])]:
        run('build_'+name,['rustc','--edition=2021','-O','-C','codegen-units=1',*flags,source/'src/bin/smooth_genus_three.rs','-o',a.out/'bin'/name])
    auditor=a.base/'smooth_negative_audit';assert sha(auditor)==old['audit_binary_sha256']
    controls=a.out/'CONTROLS.tsv';controls.write_text(''.join(f"{v['summary']['degree']}\t{v['summary']['states']}\t{v['summary']['seed']}\n"for v in old['old_controls']))
    run('full_validation',[sys.executable,PROJECT/'scripts/genus3/validate_genus_three_math.py','--binary',a.out/'bin/wide','--compact-binary',a.out/'bin/compact',
        '--auditor',auditor,'--tests',a.out/'bin/tests','--out',a.out/'checks','--controls',controls,'--source',source])
    run('joint_regressions',[a.out/'bin/tests','joint_completion::tests','--test-threads','1'])
    env=dict(os.environ,G3_TRIE_FIXTURE=str(a.review/'TRIE_FIXTURE.tsv'),G3_TRIE_WORDS=str(Path(json.loads((a.base/'CHECKPOINT.json').read_text())['word_bank'])))
    run('trie_benchmark',[a.out/'bin/tests','--ignored','--nocapture','benchmark_cached_trie_indices','--test-threads','1'],env)
    # Isolate the trie change inside the same compiled binary and proposal law.
    order=random.Random(208);blocks=[]
    for n,lane,attempts in [(18,'partial-fiber',300000),(24,'full-fiber',300000),(33,'general-fiber',100000),(21,'collision',1000000),(21,'mutation',1000000)]:
        for rep in range(5):
            modes=['literal','cached'];order.shuffle(modes);pair={}
            for mode in modes:
                folder=a.out/f'{n}_{lane}_{rep}_{mode}';binary=a.out/'bin'/('compact'if n<=25 else 'wide')
                cmd=[binary,'--d',str(n),'--lane',lane,'--seed',str(42821+rep),'--attempts',str(attempts),'--seconds','600',
                     '--closure-seconds','15','--cap','1000000','--words',env['G3_TRIE_WORDS'],'--out',folder]
                if mode=='literal':cmd+=['--literal-trie']
                if lane in ['partial-fiber','full-fiber']:cmd+=['--skip-covered-fiber']
                if lane=='mutation':cmd+=['--mutation-seeds',a.base/'MUTATION_SEEDS.tsv']
                wall=run(folder.name,cmd);checks=audit_guided_native(folder,auditor,native_components=True)
                with (folder/'COUNTS.tsv').open()as f:c=next(csv.DictReader(f,delimiter='\t'))
                with (folder/'STAGES.tsv').open()as f:s=next(csv.DictReader(f,delimiter='\t'))
                checks['isotrivial_exclusions']=audit_exclusions(folder/'ISOTRIVIAL.tsv',int(s['isotrivial_outside_target']))
                write(folder/'CHECKS.json',checks)
                pair[mode]=dict(wall=wall,counts=c,stages=s,leads_sha256=sha(folder/'LEADS.tsv'),resolutions_sha256=sha(folder/'RESOLUTIONS.tsv'),checks_sha256=sha(folder/'CHECKS.json'))
            assert {k:v for k,v in pair['literal']['counts'].items()if k!='seconds'}=={k:v for k,v in pair['cached']['counts'].items()if k!='seconds'}
            for key in ['stages','leads_sha256','resolutions_sha256']:assert pair['literal'][key]==pair['cached'][key],(n,lane,key)
            blocks.append(dict(degree=n,lane=lane,replicate=rep,results=pair));print(n,lane,rep,round(pair['literal']['wall']/pair['cached']['wall'],3),flush=True)
    focused=[]
    for rep in range(3):
        folder=a.out/f'focused13_{rep}';cmd=[a.out/'bin/compact','--d','13','--lane','collision',
            '--index-profile','6,6,6,12','--seed',str(913081+rep),'--seconds','2',
            '--trace-early','1000','--words',env['G3_TRIE_WORDS'],'--closure-seconds','15','--cap','1000000','--out',folder]
        wall=run(folder.name,cmd);checked=audit_guided_native(folder,auditor,native_components=True);write(folder/'CHECKS.json',checked)
        with (folder/'COUNTS.tsv').open()as f:counts=next(csv.DictReader(f,delimiter='\t'))
        focused.append(dict(command=list(map(str,cmd)),seconds=wall,counts=counts,native=checked['independent_native'],checks_sha256=sha(folder/'CHECKS.json')))
    validation=json.loads((a.out/'checks/VALIDATION.json').read_text());validation['joint_regression_log_sha256']=sha(logs/'joint_regressions.out');validation['full_command_receipt']=str(a.out/'RESULT.json')
    write(a.out/'VALIDATION.json',validation)
    write(a.out/'RESULT.json',dict(passed=True,events=events,paired=blocks,focused=focused,binary_sha256=sha(a.out/'bin/wide'),compact_binary_sha256=sha(a.out/'bin/compact'),
        auditor_sha256=sha(auditor),log_sha256={str(q.relative_to(a.out)):sha(q)for q in logs.iterdir()},
        qualification='same-law trie-only comparison; all first witnesses,counters,RNG and journals equal; CSP experimental only,not admitted'))

if __name__=='__main__':main()
