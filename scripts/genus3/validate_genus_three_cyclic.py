#!/usr/bin/env python3
"""One-slot clean build/regression. Does not stop, launch, or publish research."""
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

from run_genus_three import sha, write

PROJECT = Path(__file__).resolve().parents[2]
OVERLAYS = ('src/smooth_genus_three.rs', 'src/smooth_cyclic_screen.rs',
            'src/bin/smooth_negative_audit.rs')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('base', type=Path)
    p.add_argument('out', type=Path)
    a = p.parse_args()
    a.out = a.out.resolve()
    a.out.mkdir()
    source = a.out/'source'
    logs = a.out/'logs'
    logs.mkdir()
    (a.out/'bin').mkdir()
    old = json.loads((a.base/'VALIDATION.json').read_text())
    for rel, digest in old['compiled_source_files'].items():
        src = a.base/'source'/rel
        if sha(src) != digest:
            raise ValueError('base source binding changed: '+rel)
        dest = source/rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, dest)
    # Never compile the unrelated user-dirty main.rs.
    for rel in OVERLAYS:
        dest = source/rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(PROJECT/rel, dest)
    for rel in ('evidence/certificates/CERTIFIED_A25_weak.txt',
                'evidence/certificates/CERTIFIED_A25_order3_weak.txt',
                'evidence/spectral_census/s18_position_20260926/NEGATIVE_VECTOR.tsv'):
        dest = source/rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(PROJECT/rel, dest)
    events = []

    def run(name, cmd):
        started = time.perf_counter()
        with (logs/(name+'.out')).open('w') as out, (logs/(name+'.err')).open('w') as err:
            subprocess.run(list(map(str, cmd)), check=True, stdout=out, stderr=err,
                           env=dict(os.environ, CARGO_MANIFEST_DIR=str(source)))
        events.append(dict(name=name, command=list(map(str, cmd)), seconds=time.perf_counter()-started))

    for name, flags, entry in (
        ('compact', [], 'smooth_genus_three.rs'),
        ('wide', ['--cfg', 'smooth_wide_keys'], 'smooth_genus_three.rs'),
        ('tests', ['--cfg', 'smooth_wide_keys', '--test'], 'smooth_genus_three.rs'),
        ('auditor', [], 'smooth_negative_audit.rs'),
    ):
        run('build_'+name, ['rustc', '--edition=2021', '-O', '-C', 'codegen-units=1',
                           *flags, source/'src/bin'/entry, '-o', a.out/'bin'/name])
    run('cyclic_regressions', [a.out/'bin/tests', 'cyclic_screen_tests', '--nocapture', '--test-threads', '1'])
    controls = a.out/'CONTROLS.tsv'
    controls.write_text(''.join(f"{v['summary']['degree']}\t{v['summary']['states']}\t{v['summary']['seed']}\n" for v in old['old_controls']))
    run('full_validation', [sys.executable, PROJECT/'scripts/genus3/validate_genus_three_math.py',
        '--binary', a.out/'bin/wide', '--compact-binary', a.out/'bin/compact',
        '--auditor', a.out/'bin/auditor', '--tests', a.out/'bin/tests',
        '--out', a.out/'checks', '--controls', controls, '--source', source])
    validation = json.loads((a.out/'checks/VALIDATION.json').read_text())
    validation['cyclic_regression_sha256'] = sha(logs/'cyclic_regressions.out')
    write(a.out/'VALIDATION.json', validation)
    write(a.out/'BUILD.json', dict(passed=True, events=events, overlays=list(OVERLAYS),
        baseline=str(a.base.resolve()), preserved_main_sha256=sha(source/'src/main.rs'),
        source_sha256=validation['compiled_source_files'],
        qualification='one-slot clean build and controls, not a discovery or runtime speed claim'))
    print(json.dumps(dict(passed=True, out=str(a.out))))


if __name__ == '__main__':
    main()
