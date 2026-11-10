#!/usr/bin/env python3
"""Bounded paired CPU/Metal controls using an already frozen campaign runtime.

Requires the campaign stopped; does not launch a second discovery pool. Full
CPU shadow auditing, saved flag equality and input SHA256 equality are required.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import fcntl
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT))
import run_campaign as runner


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--campaign', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    root = args.campaign.resolve()
    with (runner.CONTROL / '.campaign.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        runtime = root / 'hybrid_runtime'
        for name, expected in json.loads((runtime / 'NATIVE_HASHES.json').read_text()).items():
            assert hashlib.sha256((runtime / name).read_bytes()).hexdigest() == expected
        plan = json.loads((runtime / 'PLAN.json').read_text())
        for name, expected in plan['source_hashes'].items():
            assert hashlib.sha256((runtime / 'source' / name).read_bytes()).hexdigest() == expected
        assert plan['screen_mode'] == 2
        args.out.mkdir(parents=True, exist_ok=False)
        results = []
        binary = root / 'bin' / json.loads((root / 'STATE.json').read_text())['binary_sha256']
        for d, family in ((9, 0), (18, 3)):
            arms = []
            for backend in ('cpu', 'gpu'):
                target = args.out / f'd{d}_f{family}_{backend}'
                target.mkdir()
                command = list(map(str, [runtime / 'hybrid_pilot', runtime / 'libscreen.dylib',
                               runtime / 'source/scripts/gpu/gpu_screen.metal', target, root / 'known',
                               backend, d, family, 4096, 2, 2, 20261010 + d, 5000, 25, 2, 1]))
                runner.atomic_json(target / 'COMMAND.json', command)
                begin = time.monotonic()
                with (target / 'RUN.log').open('w') as log:
                    subprocess.run(command, cwd=PROJECT, stdout=log, stderr=log, check=True, timeout=180)
                summary = json.loads((target / 'SUMMARY.json').read_text())
                assert summary['shadow_mismatches'] == 0 and len(summary['rows']) == 2
                for cert in target.glob('CERTIFIED_*.txt'):
                    subprocess.run([str(binary), 'verify', '--certificate', str(cert)], check=True, stdout=subprocess.DEVNULL)
                arms.append(dict(backend=backend, wall_seconds=time.monotonic() - begin,
                                 pipeline_seconds=summary['pipeline_seconds'],
                                 shadow_audit_seconds=summary['shadow_audit_seconds']))
            cpu = args.out / f'd{d}_f{family}_cpu'
            gpu = args.out / f'd{d}_f{family}_gpu'
            for batch in range(2):
                prefix = f'batch_{batch:04d}'
                assert json.loads((cpu / f'{prefix}.input.json').read_text()) == json.loads((gpu / f'{prefix}.input.json').read_text())
                assert (cpu / f'{prefix}.flags.bits').read_bytes() == (gpu / f'{prefix}.flags.bits').read_bytes()
            results.append(dict(degree=d, family=family, matched_inputs=True, matched_all_flags=True, arms=arms))
        runner.atomic_json(args.out / 'CHECKS.json', dict(cases=results, shadow_mismatches=0,
            scope='Two paired technical controls, 4096 proposals x two batches per arm; not a throughput benchmark or a degree census.'))
        print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
