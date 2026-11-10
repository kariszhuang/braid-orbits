#!/usr/bin/env python3
"""Sequential finite checks for the small-passport mission, with frozen evidence."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

PROJECT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(PROJECT))
from run_campaign import atomic_json, discovery_command, exact_coverage, require_results_volume, small_passport_plan


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--binary',type=Path,default=PROJECT/'target/smart/release/hurwitz_breaker')
    args=parser.parse_args();require_results_volume(args.out);args.out.mkdir(parents=True,exist_ok=False)
    binary=args.out/'bin/hurwitz_breaker';binary.parent.mkdir();shutil.copy2(args.binary,binary)
    known=args.out/'known';known.mkdir()
    for path in sorted((PROJECT/'evidence/certificates').glob('CERTIFIED_*weak.txt')):shutil.copy2(path,known/path.name)
    inputs=[*sorted((PROJECT/'src').rglob('*.rs')),PROJECT/'Cargo.toml',PROJECT/'Cargo.lock',PROJECT/'scripts/weak/run_campaign.py',
            PROJECT/'scripts/weak/test_campaign_runner.py',PROJECT/'scripts/weak/analyze_escape_words.py',PROJECT/'scripts/weak/odd_order_spin.py',Path(__file__)]
    document=dict(binary_sha256=sha(binary),source_sha256={},jobs=[],
        scope='bounded sequential technical validation and selected passport timing; no degree census or discovery ETA')
    for path in inputs:
        relative=path.relative_to(PROJECT);dest=args.out/'source'/relative;dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(path,dest);document['source_sha256'][str(relative)]=sha(dest)
    plan=small_passport_plan();atomic_json(args.out/'SMALL_PASSPORT_PLAN.json',plan)
    document['population']=dict(passports=len(plan['rows']),excluded=len(plan['excluded']),
        by_group_genus={f'{g}{genus}':n for (g,genus),n in Counter((r['group'],r['genus']) for r in plan['rows']).items()},
        plan_sha256=sha(args.out/'SMALL_PASSPORT_PLAN.json'))
    def run(name,cmd,expected=0):
        cmd=list(map(str,cmd));start=time.monotonic();log=args.out/f'{name}.log'
        with log.open('x') as stream:result=subprocess.run(cmd,cwd=PROJECT,stdout=stream,stderr=subprocess.STDOUT,timeout=180)
        row=dict(name=name,command=cmd,returncode=result.returncode,expected_code=expected,
                 seconds=time.monotonic()-start,log_sha256=sha(log))
        document['jobs'].append(row);atomic_json(args.out/'CHECKS.json',document)
        if result.returncode!=expected:raise RuntimeError(f'{name} failed: {log}')
        print(name,result.returncode,flush=True)
        return row
    run('runner_checks',[sys.executable,PROJECT/'scripts/weak/test_campaign_runner.py'])
    run('escape_selftest',[sys.executable,PROJECT/'scripts/weak/analyze_escape_words.py','--selftest'])
    for cert in sorted(known.glob('CERTIFIED_*weak.txt')):run('verify_'+cert.stem,[binary,'verify','--certificate',cert])
    targets=plan['rows'][:8]
    # A higher-genus passport also checks a finite count followed by a resumed tail.
    dense=next(r for r in plan['rows'] if r['group']=='A' and r['genus']==4 and
               sorted(r['passport'].split('/'))==sorted(['3.3.3']*3+['2.2.5']))
    targets+=[dense]
    # Pilot more than the cheapest tail before giving a throughput estimate.
    for group,genus,fraction in [('A',1,0.5),('A',1,0.9),('A',1,1.0),('S',0,0.0),('S',0,0.5),('S',1,0.5)]:
        population=sorted((r for r in plan['rows'] if r['group']==group and r['genus']==genus),
                          key=lambda r:(r['estimated_pairs'],r['name']))
        target=population[min(len(population)-1,int(fraction*(len(population)-1)))]
        if target not in targets:targets.append(target)
    for i,row in enumerate(targets):
        folder=args.out/f'pilot_{i:02d}';folder.mkdir();attempt=folder/'attempt_000';attempt.mkdir()
        cmd=discovery_command(row,binary,attempt,15,known,known/'CERTIFIED_A25_weak.txt',dict(next=0),seed=20261009)
        if i==8:cmd[cmd.index('--count')+1]='128'
        atomic_json(attempt/'SCAN_INPUT.json',dict(binary_sha256=sha(binary),passport=row['passport'],enumeration='centralizer-v1'))
        entry=run(f'pilot_{i:02d}',cmd);(attempt/'VERIFIED').touch()
        entry['coverage']=exact_coverage(folder,sha(binary),row['passport'])
        entry['summary']=dict(line.split(' = ',1) for line in (attempt/'SUMMARY.txt').read_text().splitlines() if ' = ' in line)
        if i==8 and not entry['coverage']['enumeration_complete']:
            second=folder/'attempt_001';second.mkdir();cursor=entry['coverage']['next']
            cmd=discovery_command(row,binary,second,15,known,known/'CERTIFIED_A25_weak.txt',dict(next=cursor),seed=20261009)
            cmd[cmd.index('--count')+1]='128'
            atomic_json(second/'SCAN_INPUT.json',dict(binary_sha256=sha(binary),passport=row['passport'],enumeration='centralizer-v1'))
            resumed=run('resume_dense',cmd);(second/'VERIFIED').touch()
            resumed['coverage']=exact_coverage(folder,sha(binary),row['passport'])
            assert resumed['coverage']['next']>=cursor
        proofs=attempt/'POSITIVE_WITNESSES.tsv'
        run(f'words_{i:02d}',[sys.executable,PROJECT/'scripts/weak/analyze_escape_words.py','--table',proofs,'--out',folder/'WORD_AUDIT.json'])
        atomic_json(args.out/'CHECKS.json',document)
    row=dict(kind='passport',group='S',d=4,passport='4/4/4/4',cache_positive_witnesses=True)
    folder=args.out/'zero_rows';folder.mkdir();attempt=folder/'attempt_000';attempt.mkdir()
    cmd=discovery_command(row,binary,attempt,1,known,known/'CERTIFIED_A25_weak.txt',dict(next=0))
    atomic_json(attempt/'SCAN_INPUT.json',dict(binary_sha256=sha(binary),passport=row['passport'],enumeration='centralizer-v1'))
    entry=run('zero_rows',cmd);(attempt/'VERIFIED').touch()
    entry['coverage']=exact_coverage(folder,sha(binary),row['passport'])
    assert entry['coverage']['total']==0 and entry['coverage']['search_complete']
    guard=args.out/'bad_limit'
    cmd=discovery_command(row,binary,guard,1,known,known/'CERTIFIED_A25_weak.txt',dict(next=0))
    cmd[cmd.index('--probe-witness-limit')+1]='1025';run('bad_limit_guard',cmd,2)
    document['all_expected_exits']=True;atomic_json(args.out/'CHECKS.json',document)
    print(json.dumps(document['population']))


if __name__=='__main__':main()
