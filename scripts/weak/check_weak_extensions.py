#!/usr/bin/env python3
"""Bounded sequential technical validation beside the unchanged main campaign.

Never acquires/rewrites campaign state, launches another pool, or inherits a
live candidate file. One child at a time, fixed small budgets, fresh outputs,
immutable binary/source snapshots. No unbounded launch is hidden in this tool.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

PROJECT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(PROJECT))
from run_campaign import require_results_volume


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary',type=Path,default=PROJECT/'target/smart/release/hurwitz_breaker')
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--graphs',action='store_true',help='also perform exact saved-graph analysis in pinned SymPy')
    parser.add_argument('--stability',action='store_true',help='also validate policy-ten coupled actions, broader libraries, retirement and cross-job positive proofs')
    args=parser.parse_args();require_results_volume(args.out)
    args.out.mkdir(parents=True,exist_ok=False)
    binary=args.out/'bin/hurwitz_breaker';binary.parent.mkdir();shutil.copy2(args.binary,binary)
    sources=[*sorted((PROJECT/'src').rglob('*.rs')),PROJECT/'Cargo.toml',PROJECT/'Cargo.lock',
             PROJECT/'scripts/weak/run_campaign.py',Path(__file__),PROJECT/'scripts/weak/analyze_weak_mechanisms.py',
             PROJECT/'scripts/weak/weak_graph_extensions.py',PROJECT/'scripts/weak/odd_order_spin.py',
             PROJECT/'scripts/weak/remaining_weak_passports.py',PROJECT/'scripts/weak/verify_weak_graph_extensions.py']
    inputs={str(p.relative_to(PROJECT)):sha(p) for p in sources}
    for p in sources:
        dest=args.out/'source'/p.relative_to(PROJECT);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dest)
    certs=sorted((PROJECT/'evidence/certificates').glob('*weak.txt'))
    known=args.out/'known';known.mkdir()
    for cert in certs:shutil.copy2(cert,known/cert.name)
    plan=dict(started=dt.datetime.now(dt.timezone.utc).isoformat(),binary_sha256=sha(binary),source_sha256=inputs,
              known_sha256={p.name:sha(p) for p in known.iterdir()},jobs=[],
              scope='bounded single-worker construction checks and saved-graph analysis; no degree census or throughput/discovery-time claim')
    def save():
        # This generated runtime manifest is replaced atomically, not a source edit.
        temporary=args.out/'CHECKS.json.tmp'
        with temporary.open('w') as stream:json.dump(plan,stream,indent=2);stream.write('\n')
        os.replace(temporary,args.out/'CHECKS.json')
    save()
    def run(name,command,timeout=180,expected_code=0):
        command=list(map(str,command));log=args.out/f'{name}.log'
        start=time.monotonic()
        with log.open('x') as stream:
            result=subprocess.run(command,cwd=PROJECT,stdout=stream,stderr=subprocess.STDOUT,timeout=timeout)
        entry=dict(name=name,command=command,returncode=result.returncode,expected_code=expected_code,seconds=time.monotonic()-start,
                   log=str(log),log_sha256=sha(log))
        summary=args.out/name/'SUMMARY.txt'
        if summary.exists():entry['summary']=dict(line.split(' = ',1) for line in summary.read_text().splitlines() if ' = ' in line)
        plan['jobs'].append(entry);save();print(name,result.returncode,flush=True)
        if result.returncode!=expected_code:raise RuntimeError(f'{name} unexpected return code; preserved {log}')
    if args.stability:run('project_checks',['make','check'])
    common=[binary,'smart-hunt','--weak-only','--threads',1,'--seconds',2,'--seeds',128,
            '--candidate-seconds',1,'--steps',32,'--cap',10000,'--seed',20261004,
            '--screen-steps',32,'--known-certificates',known,'--skip-known-starts',
            '--local-block-max',6,'--local-combination-cap',4096,'--local-frame-cap',128]
    for family,group,d,genus in [('independent-symmetry','S',9,1),('independent-symmetry','A',10,1),
                               ('independent-glue','S',10,1),('double-coset-glue','A',11,2),
                               ('coordinated-glue','S',9,1),('coordinated-glue','A',12,2)]:
        name=f'{family}_{group}{d}'
        run(name,common+['--family',family,'--group',group,'--d',d,'--glue-genus',genus,
                         '--pair-genus-max',2,'--boundary-centralizer-cap',512,
                         '--coordinated-rounds',2,'--repair-beam',4,'--out',args.out/name])
    for group,relaxed in [('A',False),('S',True)]:
        name=f'coordinated_fiber_{group}18_relaxed{int(relaxed)}'
        run(name,common+['--family','coordinated-fiber','--group',group,'--d',18,
            '--certificate',known/'CERTIFIED_S18_genus1_weak.txt','--coordinated-rounds',2,'--repair-beam',4,
            *(['--allow-partition-change'] if relaxed else []),'--out',args.out/name])
    name='braid_actions'
    run(name,[binary,'braid-library','--independent-actions','--d',9,'--modulus',3,'--draws',4,
              '--seed',20261004,'--out',args.out/name])
    action=next((args.out/name).glob('action_*.txt'),None)
    if action is None:raise RuntimeError('control modular representation was not constructed')
    run('explicit_rigid',[binary,'rigid-glue','--braid-action',action,'--entry-order',2,
        '--max-nodes',2000,'--seconds',2,'--cap',10000,'--out',args.out/'explicit_rigid'])
    run('small_state_actions',[binary,'braid-library','--state-actions','--d',4,'--seconds',2,'--max-tests',10000,
                               '--out',args.out/'small_state_actions'])
    if args.stability:
        for i,spec in enumerate(('9:3:2','2:1:0,9:3:2','12:4:2','9:3:2,8:2:1')):
            name=f'coupled_{i}'
            run(name,[binary,'rigid-glue','--representation-blocks',spec,'--representation-seed',71,
                '--entry-order',0,'--seconds',2,'--max-nodes',2000000,'--cap',200000,
                '--known-certificates',known,'--diagnose','--out',args.out/name])
        representation=next(j['summary']['representation key'] for j in plan['jobs'] if j['name']=='coupled_0')
        ledger=args.out/'COMPLETED_REPRESENTATIONS.tsv'
        ledger.write_text('representation_key\tentry_order\n'+representation+'\t2\n')
        rigid=[binary,'rigid-glue','--representation-blocks','9:3:2','--representation-seed',71,
            '--seconds',2,'--max-nodes',2000000,'--cap',200000,'--exclude-completed-representations',ledger]
        run('representation_restricted_skip',rigid+['--entry-order',2,'--out',args.out/'representation_restricted_skip'])
        if plan['jobs'][-1]['summary'].get('completed representation skipped')!='true':
            raise RuntimeError('completed restricted realization was not recognized')
        run('representation_unrestricted_not_skipped',rigid+['--entry-order',0,'--out',args.out/'representation_unrestricted_not_skipped'])
        if plan['jobs'][-1]['summary'].get('completed representation skipped')=='true':
            raise RuntimeError('restricted completion incorrectly excluded unrestricted entries')
        for group,d in (('S',9),('A',10),('A',11)):
            name=f'broader_{group}{d}'
            command=common.copy()
            command[command.index('--local-block-max')+1]=d
            command[command.index('--local-frame-cap')+1]=2048
            command[command.index('--local-combination-cap')+1]=100000
            run(name,command+['--family','independent-symmetry','--group',group,'--d',d,
                '--pure-braid-start','--pure-power',2,'--point-orbit-max',d,
                '--cache-positive-witnesses','--out',args.out/name])
        for cap in (1,512):
            name=f'retire_cap{cap}'
            run(name,common+['--family','double-coset-glue','--group','S','--d',10,
                '--glue-genus',1,'--pair-genus-max',2,'--boundary-centralizer-cap',cap,
                '--retire-boundary-plans','--cache-positive-witnesses','--out',args.out/name])
        seed='[(2 3)(5 7)(6 8)(9 10), (3 4)(5 6)(7 9)(8 10), (1 6 10)(2 9)(3 8)(4 5), (1 10 7)(2 5)(3 9)(4 8)]'
        replay=args.out/'inputs/positive_control.txt';replay.parent.mkdir();replay.write_text(seed+'\n')
        base=[binary,'smart-hunt','--weak-only','--group','S','--d',10,'--threads',1,
            '--seconds',5,'--candidate-seconds',5,'--steps',0,'--cap',100000,
            '--replay',replay,'--replay-only','--cache-positive-witnesses']
        run('positive_export',base+['--out',args.out/'positive_export'])
        words=args.out/'positive_export/POSITIVE_WITNESSES.tsv'
        if len(words.read_text().splitlines())!=2:raise RuntimeError('literal positive proof was not exported')
        run('positive_import',base+['--positive-witnesses',words,'--out',args.out/'positive_import'])
        values=plan['jobs'][-1]['summary']
        if values.get('positive witness cache hits')!='1' or values.get('hard unsolved candidates')!='0':
            raise RuntimeError('cross-job positive cache did not resolve its control')
        bad_seed=next(l[5:] for l in (known/'CERTIFIED_S18_genus1_weak.txt').read_text().splitlines() if l.startswith('seed='))
        false_word=args.out/'inputs/false_positive.tsv'
        false_word.write_text('degree\tmoves\tseed\n18\t\t'+bad_seed+'\n')
        run('false_positive_guard',[binary,'smart-hunt','--weak-only','--group','S','--d',18,
            '--seconds',1,'--threads',1,'--cache-positive-witnesses','--positive-witnesses',false_word,
            '--out',args.out/'false_positive_guard'],expected_code=2)
        if 'false positive witness' not in (args.out/'false_positive_guard.log').read_text():
            raise RuntimeError('invalid-word control failed for an unrelated reason')
        run('mutation_retirement_guard',common+['--family','coordinated-glue','--group','S','--d',9,
            '--glue-genus',1,'--retire-boundary-plans','--out',args.out/'mutation_retirement_guard'],expected_code=2)
        if 'only valid for unchanged, unrepaired gluing' not in (args.out/'mutation_retirement_guard.log').read_text():
            raise RuntimeError('mutation-retirement guard failed for an unrelated reason')
    for cert in certs:run('verify_'+cert.stem,[binary,'verify','--certificate',known/cert.name])
    if args.graphs:
        run('mechanism',[binary,'mechanism-audit','--certificates',','.join(str(known/p.name) for p in certs),
                         '--out',args.out/'mechanism'])
        run('extended_graphs',['uv','run','--isolated','--no-project','--python','3.13','--with','sympy==1.14.0',
            'python',args.out/'source/scripts/weak/analyze_weak_mechanisms.py','--audit',args.out/'mechanism','--extended',
            '--out',args.out/'EXTENDED_GRAPHS.json'],timeout=180)
        run('verify_extended_graphs',['uv','run','--isolated','--no-project','--python','3.13','--with','sympy==1.14.0',
            'python',args.out/'source/scripts/weak/verify_weak_graph_extensions.py','--audit',args.out/'mechanism',
            '--analysis',args.out/'EXTENDED_GRAPHS.json','--out',args.out/'STRUCTURAL_VERIFICATION.json'],timeout=180)
    run('spin_tests',[sys.executable,args.out/'source/scripts/weak/odd_order_spin.py','--selftest'])
    run('spin_even_order_guard',[sys.executable,args.out/'source/scripts/weak/odd_order_spin.py','--certificate',known/'CERTIFIED_A25_weak.txt',
                                 '--out',args.out/'SPIN_GUARD.json'])
    run('passport_tests',[sys.executable,args.out/'source/scripts/weak/remaining_weak_passports.py','--selftest'])
    plan['completed']=dt.datetime.now(dt.timezone.utc).isoformat();save()
    print('EXTENSIONS VALIDATED',args.out,flush=True)


if __name__=='__main__':main()
