"""Merge independently certified new closures into an existing public catalog.

Only NEW_COMPONENTS is exported; existing published representatives and their
annotations are preserved. No search is started and no files are deployed.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
from collections import Counter
import json
from pathlib import Path

from clean_smooth_catalog import clean_component
from export_research_catalog import identifier, read, sha, smooth_record


def merge(previous, additions):
    old = {e['id']: e for e in previous['components']}
    if len(old) != previous['component_count']:
        raise ValueError('duplicate or inconsistent previous components')
    seen = set()
    added = []
    for entry in additions:
        identity = entry['id']
        if identity in seen:
            raise ValueError('duplicate discovery identity')
        seen.add(identity)
        if identity in old:
            if any(entry[field] != old[identity][field]
                   for field in ['degree', 'genus', 'states', 'group_order']):
                raise ValueError('existing identity has inconsistent mathematics')
            continue
        old[identity] = entry
        added.append(identity)
    return dict(component_count=len(old), components=list(old.values())), added


def verify_bindings(bindings):
    """Recheck every immutable certificate input before publishing the export."""
    for binding in bindings:
        for name in ['state_table', 'summary', 'proof']:
            reference = binding[name]
            if sha(Path(reference['path'])) != reference['sha256']:
                raise ValueError('certificate input changed during export')


def extend(root, previous, out, discovery_snapshot=None):
    if out.exists() or out.with_suffix('.checks.json').exists():
        raise ValueError('refusing to overwrite export/evidence')
    status = read(root / 'STATE.json')
    live_index = root / 'NEW_COMPONENTS.json'
    if discovery_snapshot is not None and discovery_snapshot.resolve() == live_index.resolve():
        raise ValueError('snapshot must be separate from the live discovery index')
    if discovery_snapshot is None and status['status'] not in ['failed', 'stopped', 'complete']:
        raise ValueError('requires an inactive immutable campaign')
    # A frozen index may select already VERIFIED jobs while other jobs continue.
    # It never freezes, stops or relabels the live campaign itself.
    index_path = discovery_snapshot if discovery_snapshot is not None else live_index
    raw = index_path.read_bytes()
    discoveries = json.loads(raw)
    additions = []
    bindings = []
    for index, item in enumerate(discoveries, 1):
        folder = root / item['job'] / item['folder']
        if folder.resolve().parent.parent != root.resolve():
            raise ValueError('unsafe certificate path')
        proof = folder.parent / 'CHECKS.json'
        report = read(proof)
        matches = [c for c in report['components']
                   if c['summary'] == item['summary']]
        if len(matches) != 1:
            raise ValueError('ambiguous or missing certificate summary')
        value = dict(degree=item['degree'], states=int(item['summary']['states']),
                     representative=item['summary']['seed'], raw_folder=str(folder),
                     table_sha256=matches[0]['table_sha256'],
                     summary_sha256=sha(folder / 'SUMMARY.tsv'),
                     component_identity_key=item['key'], checks_sha256=sha(proof))
        record = smooth_record(value, root / 'smooth_negative_audit')
        entry = clean_component(record)
        if entry['id'] != identifier(item['degree'], item['key']):
            raise ValueError('identity changed during export')
        additions.append(entry)
        bindings.append(dict(id=entry['id'], **record['evidence']))
        if index % 10 == 0:
            print('VERIFIED_NEW', index, '/', len(discoveries), flush=True)
    result, added = merge(previous, additions)
    added_set = set(added)
    if raw != index_path.read_bytes():
        raise ValueError('discovery index changed during export')
    verify_bindings(bindings)
    checks = dict(source_root=str(root.resolve()), state_at_export=status,
                  discovery_index_path=str(index_path.resolve()),
                  discovery_index_sha256=sha(index_path),
                  frozen_discovery_snapshot=discovery_snapshot is not None,
                  examined=len(additions), added=len(added),
                  overlaps=len(additions) - len(added),
                  component_count=result['component_count'],
                  total_state_count=sum(e['states'] for e in result['components']),
                  added_ids=added,
                  added_by_degree=dict(Counter(e['degree'] for e in additions
                                               if e['id'] in added_set)),
                  isotrivial_counts=dict(Counter(str(e['isotrivial'])
                                                 for e in result['components'])),
                  bindings=bindings)
    out.parent.mkdir(parents=True, exist_ok=True)
    for path, value in [(out, result), (out.with_suffix('.checks.json'), checks)]:
        with path.open('x') as stream:
            json.dump(value, stream, indent=2, allow_nan=False)
            stream.write('\n')
    return checks


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--previous', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--discovery-snapshot', type=Path,
                        help='separate immutable index of VERIFIED jobs; permits a running campaign')
    args = parser.parse_args()
    checks = extend(args.root.resolve(), read(args.previous), args.out,
                    args.discovery_snapshot)
    print(json.dumps({k: v for k, v in checks.items()
                      if k not in ['bindings', 'added_ids']}, indent=2))
