#!/usr/bin/env python3
"""One bounded finite frontier; never changes the live campaign.

Freeze/build these sources, share literally validated completion libraries,
and independently audit each completed passport before recording progress.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

PROJECT = Path(__file__).resolve().parents[2]
INPUTS = ('src/genus_three_finite.rs', 'src/bin/genus_three_finite.rs',
          'src/genus_three_factors.rs',
          'src/bin/genus_three_finite_audit.rs',
          'scripts/genus3/verify_genus_three_finite.py', 'scripts/genus3/test_genus_three_finite.py',
          'scripts/genus3/run_genus_three_finite.py')
QUICK = [(code, 'A') for code in 'ABCDEF']
FRONTIER = QUICK + [(code, 'B') for code in 'BCDEF'] + [('C', 'C'), ('D', 'C'), ('D', 'D')]
REMAINING = [('B','B','B'),('C','B','B'),('D','B','B'),('C','B','C')]


def sha(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def save(path, value):
    tmp = path.with_suffix('.tmp')
    tmp.write_text(json.dumps(value, indent=2) + '\n')
    tmp.replace(path)


def run(root, name, command, env=None):
    with (root/'logs'/f'{name}.out').open('w') as out, (root/'logs'/f'{name}.err').open('w') as err:
        subprocess.run(command, stdout=out, stderr=err, check=True, env=env)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('out', type=Path)
    p.add_argument('--suite', choices=['quick', 'frontier', 'remaining'], default='remaining')
    p.add_argument('--workers', type=int, choices=range(1,5), default=4,
                   help='up to4 additional slots alongside the six-slot broad search')
    p.add_argument('--graphs', action='store_true', help='also repeat the pure/mixed complete graph controls and paired baseline')
    args = p.parse_args()
    root = args.out.resolve()
    if shutil.disk_usage(root.parent).free < 10 * 1024**3:
        p.error('requires at least 10 GiB free; no work launched')
    root.mkdir()
    for name in ('logs', 'jobs', 'source', 'bin'):
        (root/name).mkdir()
    for relative in INPUTS:
        dest = root/'source'/relative
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(PROJECT/relative, dest)
    bindings = {relative: sha(root/'source'/relative) for relative in INPUTS}
    cases = REMAINING if args.suite=='remaining' else [(c,'A',last) for c,last in (QUICK if args.suite=='quick' else FRONTIER)]
    plan = dict(scope='only degree15 special4x2 and named ordinary profiles; not an entire degree',
                quotient='simultaneous S15 conjugation', cases=cases, graphs=args.graphs,
                cpu_workers=args.workers, reserved_broad_slots=6, total_cpu_budget=10, gpu=False, disk_guard_GiB=10,
                public_deployment=False, live_campaign_untouched=True, source_sha256=bindings,
                compiler=subprocess.check_output(['rustc', '--version'], text=True).strip(),
                flags=['--edition=2021', '-C', 'opt-level=3', '-C', 'target-cpu=native'])
    save(root/'PLAN.json', plan)
    binary = root/'bin/genus_three_finite'
    auditor = root/'bin/genus_three_finite_audit'
    run(root, 'build', ['rustc', *plan['flags'], str(root/'source/src/bin/genus_three_finite.rs'), '-o', str(binary)])
    run(root, 'audit_build', ['rustc', *plan['flags'], str(root/'source/src/bin/genus_three_finite_audit.rs'), '-o', str(auditor)])
    tests = root/'bin/rust_tests'
    run(root, 'test_build', ['rustc', '--edition=2021', '-O', '--test', str(root/'source/src/genus_three_finite.rs'), '-o', str(tests)])
    run(root, 'rust_tests', [str(tests),'--test-threads','1'])
    env=dict(os.environ,GENUS_THREE_FINITE_BIN=str(binary),GENUS_THREE_FINITE_AUDIT_BIN=str(auditor))
    run(root, 'python_tests', [sys.executable, str(root/'source/scripts/genus3/test_genus_three_finite.py')],env=env)
    receipt = dict(status='running', binary_sha256=sha(binary), audit_binary_sha256=sha(auditor),source_sha256=bindings, jobs=[])
    save(root/'RECEIPT.json', receipt)
    verifier = root/'source/scripts/genus3/verify_genus_three_finite.py'
    libraries = {}
    jobs = [('graph_pure', ['--baseline', '--repeat', '5']), ('graph_mixed', ['--mixed'])] if args.graphs else []
    jobs += [(f'slice_{code}_{first}{last}', ['--slice', code, '--first-factor',first,'--remaining',last]) for code, first,last in cases]
    for name, options in jobs:
        assert all(sha(root/'source'/s) == h for s, h in bindings.items()), 'frozen source changed'
        folder = root/'jobs'/name
        if name.startswith('slice_'):
            key=(options[3],options[5])
            if key in libraries:
                options = options + ['--library', str(libraries[key])]
            options += ['--threads',str(args.workers)]
        if shutil.disk_usage(root).free<10*1024**3:
            raise RuntimeError('10GiB guard; completed proofs retained')
        start = time.monotonic()
        run(root, name, [str(binary), str(folder), *options])
        run(root, name+'_audit', [sys.executable, str(verifier), str(folder), '--out', str(folder/'AUDIT.json'),
                                '--native',str(auditor),'--workers',str(args.workers)])
        summary = json.loads((folder/'SUMMARY.json').read_text())
        audit = json.loads((folder/'AUDIT.json').read_text())
        if name.startswith('slice_'):
            libraries[key] = folder/'COMPLETIONS.tsv'
        entry = dict(name=name, summary=summary, audit=audit, wall_seconds=time.monotonic()-start,
                     files_sha256={f.name: sha(f) for f in folder.iterdir() if f.is_file()})
        receipt['jobs'].append(entry)
        save(root/'RECEIPT.json', receipt)
        print(f'{name}: {summary["candidate_states"]:,} states; {summary["unresolved"]} unresolved; audited {entry["wall_seconds"]:.2f}s', flush=True)
        if summary['unresolved']:
            receipt['status']='needs_full_closure';save(root/'RECEIPT.json',receipt)
            print('Short words left unresolved states; retained tokens require ordinary full braid closure',flush=True)
            return
    receipt['status'] = 'complete'
    receipt['selected_layout_states'] = sum(j['summary']['candidate_states'] for j in receipt['jobs'] if j['name'].startswith('slice_'))
    receipt['resolved_passports'] = sum(j['summary']['unresolved'] == 0 for j in receipt['jobs'] if j['name'].startswith('slice_'))
    save(root/'RECEIPT.json', receipt)
    print(f'Complete: {receipt["resolved_passports"]} passports, {receipt["selected_layout_states"]:,} selected-layout states', flush=True)


if __name__ == '__main__':
    main()
