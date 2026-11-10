"""Export certified stopped genus-three discoveries without starting a search."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path

from clean_smooth_catalog import clean_component
from export_research_catalog import read, sha, smooth_record
from extend_public_smooth_catalog import merge, verify_bindings


def validate_folder(root, item):
    folder = Path(item['raw_folder']).resolve()
    if folder.parent.parent != root.resolve() / 'jobs':
        raise ValueError('certificate is not in the source campaign jobs')
    return folder


def export(root, previous_path, out, workers=2):
    if out.exists() or out.with_suffix('.checks.json').exists():
        raise ValueError('refusing to overwrite an export')
    state = read(root / 'STATE.json')
    if state['status'] not in ('stopped', 'complete', 'failed') or state.get('active_jobs'):
        raise ValueError('requires an inactive source campaign')
    index = root / 'NEW_COMPONENTS.json'
    index_hash = sha(index)
    previous_hash = sha(previous_path)
    previous = read(previous_path)
    items = read(index)
    if len({item['id'] for item in items}) != len(items):
        raise ValueError('duplicate source IDs')

    def convert(item):
        validate_folder(root, item)
        record = smooth_record(item, root / 'smooth_negative_audit')
        entry = clean_component(record)
        if entry['id'] != item['id'] or entry['genus'] != 3:
            raise ValueError('discovery identity/genus changed')
        if 'original_isotrivial' in item and entry['isotrivial'] != item['original_isotrivial']:
            raise ValueError('classification changed; requires review')
        return entry, dict(id=entry['id'], **record['evidence'])

    entries, bindings = [], []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for entry, binding in pool.map(convert, items):
            entries.append(entry)
            bindings.append(binding)
            if len(entries) % 25 == 0:
                print('VERIFIED', len(entries), '/', len(items), flush=True)
    result, added = merge(previous, entries)
    verify_bindings(bindings)
    if sha(index) != index_hash or sha(previous_path) != previous_hash:
        raise ValueError('publication input changed during export')
    added_set = set(added)
    checks = dict(passed=True, source_root=str(root.resolve()), state=state,
        discovery_index=dict(path=str(index.resolve()), sha256=index_hash),
        previous=dict(path=str(previous_path.resolve()), sha256=previous_hash),
        examined=len(entries), added=len(added), overlaps=len(entries)-len(added),
        added_ids=added, component_count=result['component_count'],
        total_state_count=sum(e['states'] for e in result['components']),
        additions_by_degree=dict(Counter(e['degree'] for e in entries if e['id'] in added_set)),
        added_states=sum(e['states'] for e in entries if e['id'] in added_set),
        added_isotrivial=dict(Counter(str(e['isotrivial']) for e in entries if e['id'] in added_set)),
        bindings=bindings,
        scope='Certified full six-move components; existing closure receipts rechecked, displayed arithmetic/block lattice recomputed; not new independent source-group orders or a census.')
    out.parent.mkdir(parents=True, exist_ok=True)
    for path, data in [(out, result), (out.with_suffix('.checks.json'), checks)]:
        with path.open('x') as stream:
            json.dump(data, stream, indent=2, allow_nan=False)
            stream.write('\n')
    return checks


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', required=True, type=Path)
    parser.add_argument('--previous', required=True, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    parser.add_argument('--workers', type=int, choices=(1, 2), default=2)
    args = parser.parse_args()
    checks = export(args.root.resolve(), args.previous.resolve(), args.out, args.workers)
    print(json.dumps({k:v for k,v in checks.items() if k not in ('bindings','added_ids')}, indent=2))
