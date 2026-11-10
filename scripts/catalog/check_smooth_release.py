"""Verify a versioned smooth publication locally, optionally against public URLs."""
import argparse
from collections import Counter
import hashlib
import gzip
import json
from pathlib import Path
from urllib.request import Request, urlopen


def sha(data):
    return hashlib.sha256(data).hexdigest()


def check(release, previous, base_url=None):
    raw = (release / 'all-smooth.json').read_bytes()
    data = json.loads(raw)
    manifest = json.loads((release / 'manifest.json').read_bytes())
    version = data['dataset_version']
    metadata = manifest['versions'][-1]
    assert manifest['latest_version'] == metadata['dataset_version'] == version
    assert sha(raw) == metadata['snapshot']['sha256']
    assert data['component_count'] == manifest['component_count'] == len(data['components'])
    assert data['total_state_count'] == sum(e['states'] for e in data['components'])
    before = json.loads(previous.read_bytes())
    assert before['dataset_version'] == version - 1
    delta = json.loads((release / f'changes/{version}.json').read_bytes())
    assert delta['previous_version'] == version - 1
    assert delta['dataset_version'] == version
    assert delta['component_count'] == data['component_count']
    assert delta['total_state_count'] == data['total_state_count']
    state = {e['id']: e for e in before['components']}
    state.update((e['id'], e) for e in delta['upserted'])
    for entry in delta['deleted']:
        state.pop(entry['id'])
    assert state == {e['id']: e for e in data['components']}
    assert len(state) == len(data['components'])
    for entry in data['components']:
        assert entry['genus'] > 0 and entry['all_states_smooth'] is True
        assert entry['defect_distribution'] == [dict(defect=0, states=entry['states'])]
        assert isinstance(entry['primitive'], bool)
        assert entry['isotrivial'] is None or isinstance(entry['isotrivial'], bool)
    endpoints = {}
    if base_url:
        files = {'/all-smooth.json': release / 'all-smooth.json',
                 '/manifest.json': release / 'manifest.json',
                 metadata['snapshot']['url']: release / metadata['snapshot']['url'].lstrip('/'),
                 metadata['changes']['url']: release / metadata['changes']['url'].lstrip('/'),
                 '/all-smooth-assets/app.mjs': release / 'ui/app.mjs',
                 '/all-smooth-assets/catalog.mjs': release / 'ui/catalog.mjs',
                 '/sources/isotrivial-proof.pdf': release / 'sources/isotrivial-proof.pdf'}
        for endpoint, local in files.items():
            request = Request(base_url.rstrip('/') + endpoint + '?release=' + str(version),
                              headers={'User-Agent': 'python-requests/catalog-verifier',
                                       'Accept-Encoding': 'gzip',
                                       'Cache-Control': 'no-cache'})
            with urlopen(request, timeout=60) as response:
                wire = response.read()
                body = gzip.decompress(wire) if response.headers.get('Content-Encoding')=='gzip' else wire
                assert response.status == 200
            assert body == local.read_bytes(), endpoint
            endpoints[endpoint] = dict(status=200, sha256=sha(body), bytes=len(body),wire_bytes=len(wire))
    return dict(passed=True, dataset_version=version,
                component_count=data['component_count'], total_state_count=data['total_state_count'],
                upserted=len(delta['upserted']), deleted=len(delta['deleted']),
                classification_counts=dict(Counter(e['isotriviality']['status'] for e in data['components'])),
                snapshot_sha256=sha(raw),
                change_sha256=sha((release / f'changes/{version}.json').read_bytes()),
                delta_reconstructs_snapshot=True, endpoints=endpoints)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--release', type=Path, required=True)
    parser.add_argument('--previous', type=Path, required=True)
    parser.add_argument('--base-url')
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    result = check(args.release, args.previous, args.base_url)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        with args.out.open('x') as stream:
            json.dump(result, stream, indent=2)
            stream.write('\n')
    print(json.dumps(result, indent=2))
