"""Package one current snapshot plus selected deltas, with no historical snapshots.

The verified internal release remains intact. The manifest explicitly marks
unavailable historical files and tells older clients to bootstrap the current
snapshot. The current numbered URL is a hard link, not another stored copy.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
from build_public_smooth_catalog import encode, write_json
from package_math_ui import package


def digest(data):return hashlib.sha256(data).hexdigest()


def compact(source,out,changes_from=None):
    if out.exists():raise ValueError('refusing to overwrite a packaged release')
    raw=(source/'all-smooth.json').read_bytes();snapshot=json.loads(raw)
    manifest=json.loads((source/'manifest.json').read_bytes())
    latest=manifest['latest_version']
    if latest!=snapshot['dataset_version']:raise ValueError('snapshot/version mismatch')
    if [v['dataset_version'] for v in manifest['versions']]!=list(range(1,latest+1)):
        raise ValueError('noncontiguous version metadata')
    if digest(raw)!=manifest['versions'][-1]['snapshot']['sha256']:
        raise ValueError('snapshot/hash mismatch')
    first=latest if changes_from is None else changes_from
    if not 1<=first<=latest+1:raise ValueError('invalid delta retention boundary')
    retained=[]
    for v in manifest['versions']:
        n=v['dataset_version']
        if v['snapshot']['url']!=f'/versions/{n}.json' or v['changes']['url']!=f'/changes/{n}.json':
            raise ValueError('unsafe or unexpected version URL')
        if n>=first:
            name=f'changes/{n}.json';body=(source/name).read_bytes()
            if digest(body)!=v['changes']['sha256']:raise ValueError('delta/hash mismatch')
            change=json.loads(body)
            if change['dataset_version']!=n or change['previous_version']!=n-1:
                raise ValueError('delta/version mismatch')
            retained.append((name,body))
    public=copy.deepcopy(manifest)
    public.update(snapshot_retention='current-only',
        minimum_incremental_base_version=first-1,
        older_clients='Download snapshot_url if your version is below minimum_incremental_base_version.')
    for v in public['versions']:
        v['snapshot']['retained']=v['dataset_version']==latest
        v['changes']['retained']=v['dataset_version']>=first
    out.mkdir(parents=True)
    write_json(out/'all-smooth.json',raw)
    (out/'versions').mkdir()
    for suffix in ['', '.gz']:
        os.link(out/('all-smooth.json'+suffix),out/f'versions/{latest}.json{suffix}')
    for name,body in retained:write_json(out/name,body)
    write_json(out/'manifest.json',encode(public))
    with (out/'all-smooth.json.sha256').open('x') as stream:
        stream.write(digest(raw)+'  all-smooth.json\n')
    package(out)
    return public


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--out-dir',type=Path,required=True)
    p.add_argument('--changes-from',type=int,help='oldest retained change version; default: latest only')
    args=p.parse_args()
    result=compact(args.source,args.out_dir,args.changes_from)
    print(json.dumps(dict(latest_version=result['latest_version'],
        minimum_incremental_base_version=result['minimum_incremental_base_version'],
        retained_snapshots=[v['dataset_version'] for v in result['versions'] if v['snapshot']['retained']],
        retained_changes=[v['dataset_version'] for v in result['versions'] if v['changes']['retained']]),indent=2))
