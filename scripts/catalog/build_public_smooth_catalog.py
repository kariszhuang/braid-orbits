"""Publish the approved math-only catalog from a verified export; never starts searches."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
from clean_smooth_catalog import public_component


ENTRY_METADATA={'id','revision','last_updated'}
SCHEMA_VERSION=3 # nullable explicitly deferred orders, scoped quotient lists, Hodge criterion


def content(entry):return {k:v for k,v in entry.items() if k not in ENTRY_METADATA}


def timestamp(value):
    parsed=datetime.fromisoformat(value.replace('Z','+00:00'))
    if parsed.tzinfo is None:raise ValueError('timestamp must include timezone')
    return parsed


def records_by_id(records):
    result={}
    for e in records:
        identity=e.get('id','')
        if not isinstance(identity,str) or not re.fullmatch('[0-9a-f]{64}',identity) or identity in result:
            raise ValueError('missing, invalid or duplicate stable ID')
        result[identity]=e
    return result


def build(source,previous=None,now=None):
    components=[public_component(e) for e in source['components']]
    if len(components)!=source['component_count']:raise ValueError('component count mismatch')
    keys={(e['degree'],tuple(e['tuple'])) for e in components}
    if len(keys)!=len(components):raise ValueError('duplicate representative')
    current=records_by_id(components);old={};removed={}
    now=now or datetime.now(timezone.utc).isoformat(timespec='seconds');timestamp(now)
    if previous:
        if previous['schema_version'] not in (2,3) or type(previous['dataset_version']) is not int or previous['dataset_version']<1:
            raise ValueError('requires versioned previous snapshot')
        old=records_by_id(previous['components']);removed=records_by_id(previous['deleted_components'])
        if old.keys()&removed.keys():raise ValueError('active and deleted ID overlap')
        for e in list(old.values())+list(removed.values()):
            if type(e['revision']) is not int or e['revision']<1:raise ValueError('invalid previous revision')
            if timestamp(e['last_updated'])>timestamp(previous['last_updated']):raise ValueError('entry timestamp beyond snapshot')
        if current.keys()==old.keys() and all(content(e)==content(old[k]) for k,e in current.items()):
            return previous # Publishing unchanged data must not manufacture a new version.
        if timestamp(now)<=timestamp(previous['last_updated']):raise ValueError('changed snapshot must have a later timestamp')
    for identity,e in current.items():
        prior=old.get(identity)
        if prior and content(prior)==content(e):e.update(revision=prior['revision'],last_updated=prior['last_updated'])
        else:e.update(revision=(prior or removed.get(identity) or {}).get('revision',0)+1,last_updated=now)
        removed.pop(identity,None)
    for identity in old.keys()-current.keys():
        removed[identity]=dict(id=identity,revision=old[identity]['revision']+1,last_updated=now)
    components.sort(key=lambda e:(e['degree'],e['genus'],e['states'],e['id']))
    return dict(schema_version=SCHEMA_VERSION,dataset_version=previous['dataset_version']+1 if previous else 1,last_updated=now,
        component_count=len(components),total_state_count=sum(e['states'] for e in components),
        conventions=dict(point_labels='1-based',permutation_product='rightmost permutation acts first',
            states='ordered quadruples modulo simultaneous conjugation by S_d',
            scope='identity-free transitive product-one quadruples; positive genus; smooth at every state',
            contained_in='A_d when all four entries are even, otherwise S_d; not the generated group identity',
            representatives='one tuple per certified complete braid component',fixed_points_count_as_cycles=True,
            braid_moves='(a,b) -> (a b a^-1,a) and its inverse, at all three adjacent positions',
            id='SHA256 of [degree, minimum S_d canonical first-three-permutation key across the complete braid component]; rooted-traversal v1',
            revisions='increase only when published entry data changes, is removed or restored',
            last_updated='entry data-change time; snapshot time is latest dataset change',
            isotrivial='true/false follow the criterion recorded per entry over C: supplied Proposition 1/Claims 6-7 or exact full-orbit Hodge degree for complete all-smooth families; null is unknown, not false',
            group_order='exact JSON integer when computed; null with group_order_status=not_computed otherwise',
            block_quotients='recorded invariant systems; block_quotients_complete=false explicitly marks a non-exhaustive witness list'),
        components=components,deleted_components=sorted(removed.values(),key=lambda e:e['id']))


def delta(snapshot,previous=None):
    old=records_by_id(previous['components']) if previous else {}
    upserted=[e for e in snapshot['components'] if e['id'] not in old or e['revision']!=old[e['id']]['revision']]
    current=records_by_id(snapshot['components'])
    deleted=[e for e in snapshot['deleted_components'] if e['id'] in old and e['id'] not in current]
    return dict(schema_version=snapshot['schema_version'],dataset_version=snapshot['dataset_version'],
        previous_version=previous['dataset_version'] if previous else 0,last_updated=snapshot['last_updated'],
        component_count=snapshot['component_count'],total_state_count=snapshot['total_state_count'],
        upserted=upserted,deleted=deleted)


def encode(value):return (json.dumps(value,indent=2,allow_nan=False)+'\n').encode()
def digest(data):return hashlib.sha256(data).hexdigest()


def write_json(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('xb') as stream:stream.write(data)
    with path.with_name(path.name+'.gz').open('xb') as stream:stream.write(gzip.compress(data,mtime=0))


def release(source,out,previous_dir=None,now=None):
    """Build a complete immutable release, retaining every numbered snapshot/delta."""
    previous=history=None;copies=[]
    if previous_dir:
        previous=json.loads((previous_dir/'all-smooth.json').read_text())
        history=json.loads((previous_dir/'manifest.json').read_text())
        if history['latest_version']!=previous['dataset_version']:raise ValueError('history/snapshot version mismatch')
        for i,v in enumerate(history['versions'],1):
            if v['dataset_version']!=i:raise ValueError('noncontiguous history')
            for kind,folder in [('snapshot','versions'),('changes','changes')]:
                expected=f'/{folder}/{i}.json'
                if v[kind]['url']!=expected:raise ValueError('unsafe history path')
                data=(previous_dir/expected.lstrip('/')).read_bytes()
                if digest(data)!=v[kind]['sha256']:raise ValueError('changed immutable history')
                copies.append((expected.lstrip('/'),previous_dir/expected.lstrip('/')))
        if len(history['versions'])!=previous['dataset_version']:raise ValueError('missing history versions')
        if digest((previous_dir/'all-smooth.json').read_bytes())!=history['versions'][-1]['snapshot']['sha256']:
            raise ValueError('current snapshot/history hash mismatch')
    snapshot=build(source,previous,now)
    if previous and snapshot is previous:return None
    if out.exists():raise ValueError('refusing to overwrite an existing release directory')
    change=delta(snapshot,previous);version=snapshot['dataset_version'];data=encode(snapshot);changes=encode(change)
    versions=list(history['versions']) if history else []
    versions.append(dict(dataset_version=version,last_updated=snapshot['last_updated'],
        snapshot=dict(url=f'/versions/{version}.json',sha256=digest(data)),
        changes=dict(url=f'/changes/{version}.json',sha256=digest(changes)),
        upserted_count=len(change['upserted']),deleted_count=len(change['deleted'])))
    manifest=dict(schema_version=snapshot['schema_version'],latest_version=version,last_updated=snapshot['last_updated'],
        component_count=snapshot['component_count'],snapshot_url='/all-smooth.json',versions=versions)
    out.mkdir(parents=True)
    # Immutable historical bytes can be reused instead of recompressing and
    # storing every old snapshot again for each publication.
    for name,original in copies:
        target=out/name;target.parent.mkdir(parents=True,exist_ok=True)
        for suffix in ['', '.gz']:
            src=Path(str(original)+suffix);dst=Path(str(target)+suffix)
            if suffix and not src.exists():
                dst.write_bytes(gzip.compress(original.read_bytes(),mtime=0));continue
            try:os.link(src,dst)
            except OSError:shutil.copy2(src,dst)
    write_json(out/'all-smooth.json',data)
    (out/'versions').mkdir(exist_ok=True)
    for suffix in ['', '.gz']:
        os.link(out/('all-smooth.json'+suffix),out/f'versions/{version}.json{suffix}')
    write_json(out/f'changes/{version}.json',changes)
    write_json(out/'manifest.json',encode(manifest))
    with (out/'all-smooth.json.sha256').open('x') as stream:stream.write(digest(data)+'  all-smooth.json\n')
    from package_math_ui import package
    package(out)
    return manifest


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True)
    p.add_argument('--out-dir',type=Path,required=True);p.add_argument('--previous-dir',type=Path)
    args=p.parse_args()
    result=release(json.loads(args.source.read_text()),args.out_dir,args.previous_dir)
    print(json.dumps(result or dict(unchanged=True),indent=2))
