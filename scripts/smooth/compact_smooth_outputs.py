"""Losslessly pack repeated smooth certificates and finished profile outputs.

Unique certificates, unverified jobs, capped starts, ledgers and frozen inputs
stay loose. Delete only individually verified archive members. Restore never
overwrites a different file. No search is started or resumed by this script.
"""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import tarfile

PROJECT = Path(__file__).resolve().parents[2]


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def identity(item):
    return (item['degree'], tuple(tuple(p) for p in item.get('key', item.get('component_identity_key', []))))


def archive_plan(root, protected=()):
    root = root.resolve()
    protected = {Path(p).resolve() for p in protected}
    keepers = {}
    index = json.loads((root / 'NEW_COMPONENTS.json').read_text())
    for item in json.loads((root / 'PROFILE_INPUT.json').read_text())['components']:
        if 'component_identity_key' in item:
            keepers[identity(item)] = Path(item['raw_folder']).resolve()
    for item in index:
        keepers[identity(item)] = (root / item['job'] / item['folder']).resolve()
    files, directories, repeated, profiles = set(), set(), 0, 0

    def add_folder(folder):
        children = list(folder.iterdir())
        if any(not p.is_file() or p.is_symlink() or p.resolve() in protected for p in children):
            raise ValueError('unsafe or protected archive folder: ' + str(folder))
        files.update(children)
        directories.add(folder)
        sidecar = folder.with_name('._' + folder.name)
        if sidecar.is_file() and sidecar.resolve() not in protected:
            files.add(sidecar)

    for job in sorted(root.iterdir()):
        if not job.is_dir():
            continue
        if job.name.startswith('profile_'):
            report = job / 'PROFILES.json'
            if not report.exists():
                continue
            value = json.loads(report.read_text())
            if value['next'] != value['start'] + 1 or len(value['components']) != 1:
                continue
            inner = job / f"profile_{value['start']:06}.json"
            if not inner.exists() or json.loads(inner.read_text()) != value['components'][0]:
                raise ValueError('profile copies disagree: ' + str(job))
            if not all(p.name in [report.name, inner.name, '._' + report.name, '._' + inner.name]
                       for p in job.iterdir()):
                continue
            add_folder(job)
            profiles += 1
            for log in [job.with_suffix('.log'), job.with_name('._' + job.name + '.log')]:
                if log.is_file() and log.resolve() not in protected:
                    files.add(log)
            continue
        if not job.name.startswith(('guided_', 'replay_')) or not (job / 'VERIFIED').exists():
            continue
        checks = json.loads((job / 'CHECKS.json').read_text())
        for item in checks.get('closed_identities', []):
            if item['summary']['status'] != 'closed_smooth':
                continue
            folder = (job / item['folder']).resolve()
            keeper = keepers.get(identity(item))
            if keeper is None or keeper == folder or not (keeper / 'STATES.tsv').is_file():
                continue
            if not folder.is_dir() or not (folder / 'STATES.tsv').is_file():
                continue
            if not all(p.name in ['STATES.tsv', 'SUMMARY.tsv', '._STATES.tsv', '._SUMMARY.tsv']
                       for p in folder.iterdir()):
                continue
            if any(p.resolve() in protected for p in folder.iterdir()):
                continue
            add_folder(folder)
            repeated += 1
    resolved = sorted(p.resolve() for p in files)
    if len(resolved) != len(set(resolved)) or any(not p.is_relative_to(root) for p in resolved):
        raise ValueError('archive targets escape source root')
    return resolved, sorted(directories, reverse=True), dict(
        unique_new_certificates_preserved=len(index), repeated_certificate_directories=repeated,
        finished_profile_directories=profiles, archived_files=len(resolved))


def compact(root, output, files, directories):
    root = root.resolve()
    files = [p.resolve() for p in files]
    directories = [p.resolve() for p in directories]
    output.mkdir()
    archive = output / 'REPEATED_OUTPUTS.tar.gz'
    records = []
    free_before = shutil.disk_usage(root).free
    with tarfile.open(archive, 'w:gz', compresslevel=1) as tar:
        for i, path in enumerate(files):
            before = path.stat()
            digest = sha(path)
            info = tar.gettarinfo(str(path), arcname=str(path.relative_to(root)))
            if not info.isfile():
                raise ValueError('nonregular archive member')
            with path.open('rb') as stream:
                tar.addfile(info, stream)
            if path.stat().st_size != before.st_size or path.stat().st_mtime_ns != before.st_mtime_ns:
                raise ValueError('source changed during archive')
            records.append(dict(path=info.name, sha256=digest, size=before.st_size,
                                mtime_ns=before.st_mtime_ns, allocated_bytes=before.st_blocks * 512))
            if i % 1000 == 0:
                print('ARCHIVED', i, '/', len(files), flush=True)
    expected = {e['path']: e for e in records}
    seen = set()
    with tarfile.open(archive, 'r:gz') as tar:
        for member in tar:
            if not member.isfile() or member.name not in expected or member.name in seen:
                raise ValueError('unexpected archive member')
            h = hashlib.sha256()
            with tar.extractfile(member) as stream:
                for block in iter(lambda: stream.read(1024 * 1024), b''):
                    h.update(block)
            if h.hexdigest() != expected[member.name]['sha256'] or member.size != expected[member.name]['size']:
                raise ValueError('archive content mismatch; originals retained')
            seen.add(member.name)
    if seen != set(expected):
        raise ValueError('incomplete archive; originals retained')
    report = dict(root=str(root), archive=str(archive), archive_sha256=sha(archive),
                  all_members_verified_before_removal=True, members=records,
                  original_allocated_bytes=sum(e['allocated_bytes'] for e in records),
                  archive_bytes=archive.stat().st_size, free_bytes_before=free_before)
    (output / 'RECOVERY.json').write_text(json.dumps(report, indent=2) + '\n')
    removed = []
    for path, record in zip(files, records):
        current = path.stat()
        if current.st_size != record['size'] or current.st_mtime_ns != record['mtime_ns']:
            raise ValueError('source changed; remaining originals retained')
        path.unlink()
        removed.append(record['path'])
    for folder in directories:
        folder.rmdir()
    report.update(removed_files=len(removed), removed_directories=len(directories),
                  free_bytes_after=shutil.disk_usage(root).free)
    (output / 'RESULTS.json').write_text(json.dumps(report, indent=2) + '\n')
    return report


def restore(archive, root):
    root = root.resolve()
    report = json.loads((archive.parent / 'RECOVERY.json').read_text())
    if sha(archive) != report['archive_sha256']:
        raise ValueError('archive checksum changed')
    records = {r['path']: r for r in report['members']}
    with tarfile.open(archive, 'r:gz') as tar:
        members = tar.getmembers()
        seen = set()
        existed = set()
        # Check ALL existing targets before writing anything. On exFAT macOS
        # creates AppleDouble sidecars while normal files are restored. Restore
        # archived sidecars last, replacing only ones generated by this call.
        for member in members:
            target = root / member.name
            if (not member.isfile() or member.name not in records or member.name in seen
                    or not target.resolve().is_relative_to(root) or target.is_symlink()):
                raise ValueError('unsafe recovery member')
            seen.add(member.name)
            if target.exists():
                if sha(target) != records[member.name]['sha256']:
                    raise ValueError('recovery would overwrite different data')
                existed.add(member.name)
        if seen != set(records):
            raise ValueError('incomplete recovery archive')
        members.sort(key=lambda m: (Path(m.name).name.startswith('._'), -len(Path(m.name).parts), m.name))
        for member in members:
            if member.name in existed:
                continue
            target = root / member.name
            target.parent.mkdir(parents=True, exist_ok=True)
            generated_sidecar = Path(member.name).name.startswith('._') and target.exists()
            with tar.extractfile(member) as source, target.open('wb' if generated_sidecar else 'xb') as sink:
                shutil.copyfileobj(source, sink)
            if sha(target) != records[member.name]['sha256']:
                raise ValueError('recovery content mismatch')
            os.utime(target, ns=(records[member.name]['mtime_ns'], records[member.name]['mtime_ns']))
        for member in members:
            if sha(root / member.name) != records[member.name]['sha256']:
                raise ValueError('final recovery checksum mismatch')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--out', type=Path)
    parser.add_argument('--execute', action='store_true')
    parser.add_argument('--restore', type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    if root.parent != Path('data') or not root.name.startswith('smooth_fast_'):
        raise ValueError('only an exact managed smooth run is admitted')
    lock = (PROJECT / '.run-control/.campaign.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    state = json.loads((root / 'STATE.json').read_text())
    if state['status'] not in ['failed', 'stopped', 'complete']:
        raise ValueError('run is not inactive')
    launch = json.loads((root / 'LAUNCH.json').read_text())
    try:
        os.kill(launch['pid'], 0)
    except ProcessLookupError:
        pass
    else:
        raise ValueError('controller is still alive or PID has been reused')
    if args.restore:
        restore(args.restore.resolve(), root)
        print('RESTORED', flush=True)
        return
    manifest = json.loads((PROJECT / 'evidence/manifest.json').read_text())
    protected = [(PROJECT / e['path']).resolve() for e in manifest['sources']]
    files, directories, summary = archive_plan(root, protected)
    summary.update(allocated_bytes=sum(p.stat().st_blocks * 512 for p in files),
                   logical_bytes=sum(p.stat().st_size for p in files))
    print(json.dumps(summary, indent=2), flush=True)
    if args.execute:
        if args.out is None or args.out.resolve().parent != root.parent or args.out.exists():
            raise ValueError('new explicit sibling archive directory required')
        report = compact(root, args.out.resolve(), files, directories)
        report['selection'] = summary
        (args.out / 'RESULTS.json').write_text(json.dumps(report, indent=2) + '\n')
        print('COMPACTED', json.dumps({k: report[k] for k in ['archive', 'removed_files', 'removed_directories', 'free_bytes_before', 'free_bytes_after']}), flush=True)


if __name__ == '__main__':
    main()
