#!/usr/bin/env python3
"""Publish a frozen selection of independently certified catalog-growth records.

Preserves old public records. Deferred orders and partial block-system lists
are explicit; no new unverified search results or structural claims enter.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
from collections import Counter
from fractions import Fraction
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

from build_public_smooth_catalog import encode
from clean_smooth_catalog import public_component
from verify_smooth_geometry import parse_tuple, compose, cycles, orbits

REPO=Path(__file__).resolve().parents[2]

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb')as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def uncompressed_sha(path):
    h=hashlib.sha256()
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rb')as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
    return h.hexdigest()

def load_audit(run):
    path=run/'source/evidence/catalog_expansion_20261030/audit.py'
    spec=importlib.util.spec_from_file_location('frozen_growth_audit',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module,path

def convert(r,run,module,audit_path,verify_files=True,artifact_hashes=None):
    cert=r['certificate'];folder=Path(cert['folder']).resolve()
    if not folder.is_relative_to(run.resolve()):raise ValueError('certificate outside owned run')
    meta=folder/'VERIFIED.json'
    if json.loads(meta.read_text())!=r:raise ValueError('immutable verified metadata changed')
    check=cert['independent_native']
    if not check['closed'] or not check['all_six_edges_checked'] or check['states']!=r['states']:
        raise ValueError('not independently closed')
    audit_hash,auditor_hash=artifact_hashes or (sha(audit_path),sha(run/'bin/smooth_negative_audit'))
    if cert['audit_script_sha256']!=audit_hash:raise ValueError('audit source mismatch')
    if cert['auditor_sha256']!=auditor_hash:raise ValueError('auditor binary mismatch')
    paths={'state_table':folder/cert.get('states_file','STATES.tsv'),
           'independent_keys':folder/cert.get('independent_keys_file','INDEPENDENT_KEYS.tsv')}
    if verify_files:
        if uncompressed_sha(paths['state_table'])!=cert['states_sha256']:raise ValueError('state certificate hash mismatch')
        if uncompressed_sha(paths['independent_keys'])!=cert['independent_keys_sha256']:raise ValueError('canonical key certificate hash mismatch')
    n=r['degree'];t=parse_tuple('['+', '.join(r['tuple'])+']',n)
    product=tuple(range(n))
    for p in t:product=compose(product,p)
    if product!=tuple(range(n))or len(orbits(t))!=1:raise ValueError('representative not product-one/transitive')
    systems=[]
    if w:=r.get('block_witness'):
        q=parse_tuple('['+', '.join(w['tuple'])+']',w['block_count'])
        qi=[w['block_count']-len(cycles(p))for p in q]
        excess=[a-w['block_size']*b for a,b in zip(r['indices'],qi)]
        primitive,_=module.block_witness(q)
        systems.append(dict(block_count=w['block_count'],block_size=w['block_size'],
            blocks=w['blocks'],tuple=w['tuple'],genus=w['genus'],indices=qi,
            ramification_excess=excess,primitive=primitive,group_order=None,kernel_order=None,
            orders_status='not_computed'))
    h=r['isotriviality'];l12=Fraction(h['twelve_L']);L=Fraction(h['L'])
    if L!=l12/12 or (L==0)!=r['isotrivial']:raise ValueError('Hodge value/sign mismatch')
    status='isotrivial'if L==0 else'non-isotrivial'
    explanation=(f'The certified full braid orbit has exact normalized Hodge degree L = {L}. '
                 'For complete smooth source families of genus at least two over C, '
                 'L = 0 implies isotriviality and L > 0 implies non-isotriviality.')
    entry=dict(id=r['id'],degree=n,genus=r['genus'],states=r['states'],primitive=r['primitive'],
        group_order=None,group_order_status='not_computed',tuple=r['tuple'],cycle_types=r['cycle_types'],
        indices=r['indices'],defect_distribution=[dict(defect=0,states=r['states'])],
        representative_pair_orbit_sizes={f'g{i+1}_g{i+2}':sorted(map(len,orbits(t[i:i+2])),reverse=True)for i in range(3)},
        quotient_genera=sorted({b['genus']for b in systems}),block_sizes=sorted({b['block_size']for b in systems}),
        block_quotients=systems,block_quotients_complete=bool(r['primitive']),all_states_smooth=True,
        isotrivial=r['isotrivial'],
        isotriviality=dict(status=status,method='full-orbit-hodge-degree-v1',source=None,source_sha256=None,
            explanation=explanation,witnesses=[]),
        hodge_degree=dict(normalization='12L=d-sum(rho(entries))+3*mean(rho(g1*g2))',
            twelve_L=str(l12),L=str(L),pair_product_cycle_distribution=r['product_cycle_distribution']))
    result=public_component(entry)
    binding=dict(id=r['id'],metadata=dict(path=str(meta),sha256=sha(meta)),
        state_table=dict(path=str(paths['state_table']),sha256=sha(paths['state_table']),uncompressed_sha256=cert['states_sha256']),
        independent_keys=dict(path=str(paths['independent_keys']),sha256=sha(paths['independent_keys']),uncompressed_sha256=cert['independent_keys_sha256']))
    return result,binding

def export(snapshot,previous,out):
    selected=json.loads(snapshot.read_text());old=json.loads(previous.read_text())
    declared_runs=selected.get('source_runs')
    if declared_runs is None:declared_runs=[selected['source_run']]
    roots=[Path(p).resolve() for p in declared_runs]
    if not roots or len(set(roots))!=len(roots):raise ValueError('invalid or duplicate source-run roots')
    audit_cache={}
    for run in roots:
        module,path=load_audit(run)
        audit_cache[run]=(module,path,(sha(path),sha(run/'bin/smooth_negative_audit')))
    before={r['id']:r for r in old['components']};merged=dict(before);bindings=[]
    if len(selected['components'])!=selected['component_count']:raise ValueError('snapshot count mismatch')
    if len({r['id']for r in selected['components']})!=selected['component_count']:raise ValueError('duplicate discovery ID')
    for i,r in enumerate(selected['components'],1):
        folder=Path(r['certificate']['folder']).resolve()
        matches=[run for run in roots if folder.is_relative_to(run)]
        if len(matches)!=1:raise ValueError(f'certificate has no unique declared source run: {folder}')
        run=matches[0];module,path,hashes=audit_cache[run]
        entry,binding=convert(r,run,module,path,artifact_hashes=hashes)
        if entry['id']in merged:raise ValueError('unexpected overlap with public baseline')
        merged[entry['id']]=entry;bindings.append(binding)
        if i%200==0:print('VERIFIED_EXPORT',i,'/',selected['component_count'],flush=True)
    for k,r in before.items():assert merged[k]==r
    result=dict(component_count=len(merged),components=list(merged.values()))
    checks=dict(passed=True,added=selected['component_count'],old_records_preserved=len(before),
        component_count=len(merged),new_states=sum(r['states']for r in selected['components']),
        additions_by_classification=dict(Counter(r['isotriviality']['status']for r in selected['components'])),
        snapshot=dict(path=str(snapshot.resolve()),sha256=sha(snapshot)),
        previous=dict(path=str(previous.resolve()),sha256=sha(previous)),bindings=bindings,
        source_runs=[str(run) for run in roots],
        scope='Existing independent full-orbit certificates and immutable source hashes rechecked; no new closure search or full group-order/block-lattice computation.')
    with out.open('xb')as f:f.write(encode(result))
    with out.with_suffix('.checks.json').open('xb')as f:f.write(encode(checks))
    return checks

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--snapshot',type=Path,required=True);p.add_argument('--previous',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();result=export(a.snapshot,a.previous,a.out)
    print(json.dumps({k:v for k,v in result.items()if k!='bindings'},indent=2))
