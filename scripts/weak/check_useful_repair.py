#!/usr/bin/env python3
"""Locked, sequential bounded checks/paired comparisons; no competing hunt."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import csv
import fcntl
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
import run_campaign as runner
from analyze_escape_words import analyze


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--predecessor',type=Path,required=True)
    parser.add_argument('--binary',type=Path,default=ROOT/'target/smart/release/hurwitz_breaker')
    parser.add_argument('--scope-source',type=Path,help='reuse a completed frozen predecessor scope audit, with checkpoint guards and logical revalidation')
    args=parser.parse_args()
    with (runner.CONTROL/'.campaign.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        args.out.mkdir(parents=True,exist_ok=False)
        records=[];scope_reuse=None
        sources=[str(p.relative_to(ROOT)) for p in sorted((ROOT/'src').rglob('*.rs'))]+[
            str(p.relative_to(ROOT)) for p in sorted((ROOT/'scripts').rglob('*.py'))]+[
            'Cargo.toml','Cargo.lock','Makefile','scripts/weak/run_campaign.py','scripts/weak/test_campaign_runner.py','scripts/weak/test_campaign_cache.py']
        hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources}
        binary_hash=hashlib.sha256(args.binary.read_bytes()).hexdigest()
        (args.out/'source').mkdir();(args.out/'bin').mkdir()
        for p in sources:shutil.copy2(ROOT/p,args.out/'source'/f'{hashes[p]}_{Path(p).name}')
        frozen=args.out/'bin'/binary_hash;shutil.copy2(args.binary,frozen);args.binary=frozen
        predecessor_hashes={p:hashlib.sha256((args.predecessor/p).read_bytes()).hexdigest()
            for p in ('CONFIG.json','STATE.json','CHECKPOINT.json')}
        def save():
            runner.atomic_json(args.out/'CHECKS.json',dict(checks=records,source_sha256=hashes,
                binary_sha256=binary_hash,predecessor=str(args.predecessor),workers=1,
                predecessor_sha256=predecessor_hashes,
                scope_reuse=scope_reuse,
                scope='bounded controls/pilots and dependent technical comparisons; no discovery forecast'))
        def run(name,command,timeout=180):
            begin=time.monotonic()
            with (args.out/f'{name}.log').open('w') as log:
                code=subprocess.run(list(map(str,command)),cwd=ROOT,stdout=log,stderr=log,timeout=timeout).returncode
            records.append(dict(name=name,command=list(map(str,command)),returncode=code,seconds=time.monotonic()-begin));save()
            if code:raise RuntimeError(f'{name} failed; preserve {args.out}')
            print(f'{name}: passed',flush=True)
        certs=sorted((ROOT/'evidence/certificates').glob('CERTIFIED_*_weak.txt'))
        run('project_checks',['make','check'],timeout=600)
        run('spectral_regressions',['make','check-spectra'],timeout=180)
        assert len(certs)==9
        for i,cert in enumerate(certs):
            run(f'control_{i}',[args.binary,'certificate-profile','--certificate',cert,
                '--pure-stability','--out',args.out/f'control_{i}.tsv'])
            pure=(args.out/f'control_{i}.tsv.pure.txt').read_text()
            d=int(runner.certificate_header(cert)['d'])
            assert f'factors through S4 = {str(d==25).lower()}' in pure
            assert f'maximum passport fiber size = {1 if d==25 else 10}' in pure
        for group,d in [('S',10),('A',12),('S',14),('A',17)]:
            for hidden in (False,True):
                name=f'pilot_{group}{d}_{"hidden" if hidden else "admit"}'
                run(name,[args.binary,'smart-hunt','--group',group,'--d',d,'--family','joint-pure',
                    '--joint-admission-aware',*(['--joint-hidden-fiber'] if hidden else []),
                    '--joint-repair-rounds',4,'--weak-only','--all-pairs-intransitive',
                    '--threads',1,'--seconds',3,'--candidate-seconds',1,'--cap',5000,
                    '--local-block-max',min(d,11),'--local-frame-cap',128,'--point-orbit-max',d,
                    '--local-combination-cap',10000,'--screen-steps',64,'--seed',20261010,
                    '--skip-known-starts','--known-certificates',ROOT/'evidence/certificates',
                    '--cache-positive-witnesses','--construction-witness-limit',32,'--out',args.out/name])
                summary=(args.out/name/'SUMMARY.txt').read_text()
                assert 'joint admission aware = true' in summary
                if hidden:assert 'joint complete pure basis required = true' in summary
                for cert in (args.out/name).glob('CERTIFIED_*.txt'):
                    run(name+'_verify',[args.binary,'verify','--certificate',cert])
                (args.out/name/'VERIFIED').touch()
                proof=analyze([args.out/name/'POSITIVE_WITNESSES.tsv'])
                runner.atomic_json(args.out/f'{name}_ESCAPES.json',proof)
        for group,d,name in [('A',25,'CERTIFIED_A25_XXYY_genus4_weak.txt'),('S',18,'CERTIFIED_S18_genus1_weak.txt')]:
            label=f'paired_{group}{d}'
            run(label,[args.binary,'smart-hunt','--group',group,'--d',d,'--joint-benchmark',
                '--certificate',ROOT/'evidence/certificates'/name,'--benchmark-starts',32,
                '--benchmark-repeats',3,'--seed',20261010,'--out',args.out/label],timeout=300)
            with (args.out/label/'BENCHMARK.tsv').open() as stream:rows=list(csv.DictReader(stream,delimiter='\t'))
            assert len(rows)==6 and {r['arm'] for r in rows}=={'legacy','admission-aware'}
        if args.scope_source:
            source=args.scope_source;proof_path=source.parent/'CHECKS.json'
            proof=json.loads(proof_path.read_text());audit=json.loads(source.read_text())
            expected={f'control_{i}' for i in range(9)}|{f'pilot_{g}{d}_{arm}'
                for g,d in [('S',10),('A',12),('S',14),('A',17)] for arm in ('admit','hidden')}|{'paired_A25','paired_S18'}
            if (proof.get('predecessor_sha256')!=predecessor_hashes or
                Path(proof['predecessor']).resolve()!=args.predecessor.resolve() or
                not expected<={r['name'] for r in proof['checks']} or any(r['returncode'] for r in proof['checks'])):
                raise RuntimeError('Scope source lacks completed matching predecessor checks')
            exclusions={(r['group'],r['d'],tuple(sorted(r['passport'].split('/')))):r['evidence'] for r in audit['exclusions']}
            if runner.small_degree_scope_audit(exclusions)!=audit:raise RuntimeError('Scope source failed logical revalidation')
            scope_reuse=dict(source=str(source),sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                checks=str(proof_path),checks_sha256=hashlib.sha256(proof_path.read_bytes()).hexdigest(),
                scope='reuse of raw interval validation on the same unchanged stopped predecessor; not progress TSV reuse')
        else:
            exclusions=runner.inherited_passport_exclusions(args.predecessor)
            audit=runner.small_degree_scope_audit(exclusions)
        runner.atomic_json(args.out/'SMALL_DEGREE_SCOPE_AUDIT.json',audit)
        if any(hashlib.sha256((args.predecessor/p).read_bytes()).hexdigest()!=sha
               for p,sha in predecessor_hashes.items()):raise RuntimeError('Predecessor changed during scope audit')
        save();print(f'All {len(records)} bounded checks passed.',flush=True)


if __name__=='__main__':main()
