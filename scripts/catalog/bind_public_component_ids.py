"""Bootstrap stable IDs for the unversioned catalog from its complete certificates."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
from export_research_catalog import read, sha, component_identity, identifier
from run_s18_linked_roots import state_key
from verify_smooth_geometry import parse_tuple


def representative_key(entry):
    return identifier(entry['degree'],state_key(parse_tuple('['+', '.join(entry['tuple'])+']',entry['degree'])))


def bind(source,items,auditor,workers=2):
    entries={representative_key(e):dict(e) for e in source['components']}
    if len(entries)!=source['component_count']:raise ValueError('duplicate representative')
    result=[];ids=set()
    def verify(item):
        key,fresh,summary=component_identity(item,auditor);n=int(summary['degree'])
        t=parse_tuple(item.get('representative',summary['seed']),n)
        return identifier(n,state_key(t)),identifier(n,key),summary,fresh
    fresh_count=0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for lookup,identity,summary,fresh in pool.map(verify,items):
            if lookup not in entries or identity in ids:raise ValueError('missing/duplicate component binding')
            entry=entries.pop(lookup)
            if entry['states']!=int(summary['states']) or entry['genus']!=int(summary['genus']):raise ValueError('component metadata mismatch')
            entry['id']=identity;result.append(entry);ids.add(identity);fresh_count+=fresh
            if len(result)%500==0:print('IDENTIFIED',len(result),flush=True)
    if entries:raise ValueError('unbound exported records')
    result.sort(key=lambda e:(e['degree'],e['genus'],e['states'],e['id']))
    return dict(component_count=len(result),components=result),fresh_count


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['root','source','out']:p.add_argument('--'+name,type=Path,required=True)
    args=p.parse_args();root=args.root.resolve()
    if args.out.exists():raise ValueError('refusing to overwrite identified source')
    if read(root/'STATE.json')['status'] not in ['failed','complete','stopped']:raise ValueError('inactive source required')
    items=list(read(root/'PROFILE_INPUT.json')['components'])
    for e in read(root/'NEW_COMPONENTS.json'):
        folder=root/e['job']/e['folder'];report=read(folder.parent/'CHECKS.json')
        c=next(c for c in report['components'] if c['summary']==e['summary'])
        items.append(dict(raw_folder=str(folder),states=int(e['summary']['states']),representative=e['summary']['seed'],
            component_identity_key=e['key'],table_sha256=c['table_sha256'],checks_sha256=sha(folder.parent/'CHECKS.json')))
    auditor=root/'smooth_negative_audit'
    if not auditor.is_file():raise ValueError('independent auditor unavailable')
    result,fresh=bind(read(args.source),items,auditor,workers=2)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    with args.out.open('x') as stream:json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
    print('BOUND',result['component_count'],'FRESH_MINIMUM_AUDITS',fresh,'SHA256',sha(args.out),flush=True)
