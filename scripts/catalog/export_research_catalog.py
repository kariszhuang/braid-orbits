"""Export positive-genus everywhere-smooth components as filterable math data.

Full smooth graphs remain in original certificates. Every exported component
is bound to a complete independent closure audit; no discovery is launched.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import json
from math import factorial
from pathlib import Path

from audit_guided_native import audit_component_native
from compact_smooth_outputs import sha
from investigate_smooth_quotients import quotient
from run_s18_linked_roots import state_key
from smooth_block_lattice import block_lattice
from verify_smooth_geometry import parse_tuple, cycles, compose, orbits
from clean_smooth_catalog import clean_catalog



def ref(path):return dict(path=str(Path(path).resolve()),sha256=sha(Path(path)))
def read(path):return json.loads(Path(path).read_text())
def first_row(path):
    with Path(path).open() as stream:return next(csv.DictReader(stream,delimiter='\t'))
def notation(p):return ''.join('('+' '.join(str(x+1) for x in c)+')' for c in cycles(p) if len(c)>1) or '()'
def identifier(degree,key):
    return hashlib.sha256(json.dumps([degree,key],separators=(',',':')).encode()).hexdigest()


def tuple_details(t):
    n=len(t[0]);identity=tuple(range(n));product=identity
    for p in t:product=compose(product,p)
    if product!=identity or len(orbits(t))!=1:raise ValueError('not transitive product-one tuple')
    types=[sorted(map(len,cycles(p)),reverse=True) for p in t];indices=[n-len(v) for v in types]
    if sum(indices)%2:raise ValueError('nonintegral genus')
    pairs=[]
    for i in range(3):
        a,b=t[i:i+2];p=compose(a,b);parts=orbits([a,b]);chi=len(cycles(a))+len(cycles(b))+len(cycles(p))-n
        if chi%2:raise ValueError('nonintegral pair genus')
        pairs.append(dict(positions=[i+1,i+2],orbits=[[x+1 for x in c] for c in parts],
            orbit_sizes=sorted(map(len,parts),reverse=True),orbit_count=len(parts),transitive=len(parts)==1,
            product_cycle_notation=notation(p),product_cycles_including_fixed_points=[[x+1 for x in c] for c in cycles(p)],
            product_cycle_count=len(cycles(p)),product_index=n-len(cycles(p)),pair_component_genus_sum=len(parts)-chi//2))
    defect=pairs[0]['product_cycle_count']-pairs[0]['orbit_count']-pairs[2]['orbit_count']+1
    genus=1-n+sum(indices)//2;smooth=pairs[0]['product_index'] in [sum(indices[:2]),sum(indices[2:])]
    if not 0<=defect<=genus or (smooth and defect):raise ValueError('defect bounds/smooth implication')
    if genus!=defect+pairs[0]['pair_component_genus_sum']+pairs[2]['pair_component_genus_sum']:raise ValueError('genus decomposition')
    return dict(cycle_notation='['+', '.join(notation(p) for p in t)+']',
        permutations=[dict(position=i+1,cycle_notation=notation(p),images_1_based=[x+1 for x in p],
            nontrivial_cycles_1_based=[[x+1 for x in c] for c in cycles(p) if len(c)>1],
            fixed_points_1_based=[x+1 for x in range(n) if p[x]==x],cycle_lengths=types[i],
            cycle_multiplicities={str(k):v for k,v in sorted(Counter(types[i]).items(),reverse=True)},index=indices[i]) for i,p in enumerate(t)],
        genus=genus,indices=indices,passport_ordered=types,passport_unordered=sorted(types),
        product_one=True,generated_action_transitive=True,adjacent_pairs=pairs,
        defect=defect,defect_calculation=dict(product_cycles=pairs[0]['product_cycle_count'],
            left_pair_orbits=pairs[0]['orbit_count'],right_pair_orbits=pairs[2]['orbit_count'],constant=1),
        smooth_at_representative=smooth)


def block_details(t):
    lattice=block_lattice(t)
    for system in lattice['proper_systems']:
        blocks=[[x-1 for x in b] for b in system['partition']];labels={x:j for j,b in enumerate(blocks) for x in b}
        q=tuple(tuple(labels[p[b[0]]] for b in blocks) for p in t)
        result=quotient(t,blocks)
        if result['quotient_genus']!=system['genus'] or result['ramification_excess']!=system['ramification_excess']:raise ValueError('block genus reconciliation')
        system['quotient_permutations']=[dict(position=i+1,cycle_notation=notation(p),images_1_based=[x+1 for x in p],
            cycle_lengths=sorted(map(len,cycles(p)),reverse=True)) for i,p in enumerate(q)]
        system['active_block_counts']=[sum(any(p[x]!=x for x in b) for b in blocks) for p in t]
        system['relative_genus_checked']=True
    return lattice


def mechanism(lattice):
    for s in lattice['proper_systems']:
        e=s['ramification_excess']
        if s['block_size']==2 and s['genus']==2 and e==[2,0,0,0] and s['active_block_counts'][0]==2:
            return 'two-active-block ramified double lift'
        for i in s['identity_entries']:
            if s['active_block_counts'][i-1]==1 and all(not v for j,v in enumerate(e) if j!=i-1):return 'single-local/unramified'
        if s['block_size']==2 and sorted(e)==[0,0,1,1]:return 'two-simple-ramification double cover'
    return 'unclassified' # Failure of these sufficient criteria is not primitivity.


def component_identity(item,auditor):
    """Bind the stable minimum identity to a complete closure certificate."""
    folder=Path(item['raw_folder']);summary=first_row(folder/'SUMMARY.tsv');n=int(summary['degree']);count=int(summary['states'])
    if summary['status']!='closed_smooth' or count!=item['states'] or sha(folder/'STATES.tsv')!=item['table_sha256']:raise ValueError('smooth certificate binding')
    if 'summary_sha256' in item and sha(folder/'SUMMARY.tsv')!=item['summary_sha256']:raise ValueError('summary changed')
    proof_path=folder.parent/'CHECKS.json';fresh=False
    if 'component_identity_key' in item:
        key=tuple(tuple(p) for p in item['component_identity_key']);proof=read(proof_path)
        if 'checks_sha256' in item and sha(proof_path)!=item['checks_sha256']:raise ValueError('job proof changed')
        if not (folder.parent/'VERIFIED').exists():raise ValueError('unverified job')
        matches=[e for e in proof['closed_identities'] if e['folder']==folder.name and tuple(tuple(p) for p in e['key'])==key and e['summary']==summary]
        if len(matches)!=1 or not any(c['summary']==summary and c['table_sha256']==item['table_sha256'] and c['independent_six_move_closure'] for c in proof['components']):raise ValueError('independent closed identity missing')
    else:
        audit=audit_component_native(folder,auditor)
        if not audit['independent_six_move_closure'] or len(audit['state_keys'])!=count:raise ValueError('unkeyed certificate not closed')
        key=min(audit['state_keys']);fresh=True
    return key,fresh,summary


def smooth_record(item,auditor):
    key,fresh,summary=component_identity(item,auditor)
    folder=Path(item['raw_folder']);n=int(summary['degree']);count=int(summary['states'])
    proof_path=folder.parent/'CHECKS.json'
    t=parse_tuple(item.get('representative',summary['seed']),n);detail=tuple_details(t)
    if not detail['smooth_at_representative'] or detail['genus']!=int(summary['genus']) or any(p==tuple(range(n)) for p in t):raise ValueError('smooth scope/representative')
    order=int(summary['group_order']);ambient='A' if all(v%2==0 for v in detail['indices']) else 'S';full_order=factorial(n)//(2 if ambient=='A' else 1)
    if ambient!=summary['ambient'] or (order==full_order)!=(summary['full_group']=='true'):raise ValueError('group metadata mismatch')
    blocks=block_details(t)
    return dict(component_id=identifier(n,key),kind='everywhere_smooth',degree=n,genus=detail['genus'],quotient_state_count=count,
        representative=detail,group=dict(order_decimal=str(order),ambient_notation=f'{ambient}_{n}',generates_full_ambient_group=order==full_order,
            order_provenance=item.get('group_order_provenance','producer exact Schreier order; no new independent order computation during export')),
        all_states_smooth=True,defect_distribution=[dict(defect=0,states=count)],
        defect_distribution_basis='Certified smooth at every state; smooth implies zero defect under product one/transitivity.',
        block_lattice=blocks,mechanism_at_representative=mechanism(blocks),
        geometry=dict(status='not finally classified in this export',human_mathematical_review='pending'),
        prior_analysis={k:v for k,v in item.items() if k not in ['component_identity_key','raw_folder','representative','table_sha256','summary_sha256']},
        available_component_profiles=[],evidence=dict(state_table=ref(folder/'STATES.tsv'),summary=ref(folder/'SUMMARY.tsv'),
            proof=ref(proof_path) if proof_path.exists() else None,independent_six_move_closure=True,
            minimum_component_identity='fresh native export audit' if fresh else 'bound previous independent closure identity',
            identity_convention='minimum S_d canonical key across the complete component; computed from ordered first three permutations'))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--workers',type=int,default=2,choices=[1,2]);args=p.parse_args();root=args.root.resolve()
    if args.out.exists():raise ValueError('refusing to overwrite an existing export')
    if read(root/'STATE.json')['status'] not in ['failed','stopped','complete']:raise ValueError('export requires immutable inactive campaign')
    items=list(read(root/'PROFILE_INPUT.json')['components'])
    for e in read(root/'NEW_COMPONENTS.json'):
        folder=root/e['job']/e['folder'];report=read(folder.parent/'CHECKS.json')
        c=next(c for c in report['components'] if c['summary']==e['summary'])
        items.append(dict(degree=e['degree'],states=int(e['summary']['states']),representative=e['summary']['seed'],
            group_order=int(e['summary']['group_order']),raw_folder=str(folder),table_sha256=c['table_sha256'],summary_sha256=sha(folder/'SUMMARY.tsv'),
            component_identity_key=e['key'],checks_sha256=sha(folder.parent/'CHECKS.json'),group_order_provenance='producer exact Schreier order'))
    records=[];seen=set()
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for record in pool.map(lambda item:smooth_record(item,root/'smooth_negative_audit'),items):
            if record['genus']<=0:continue
            if record['component_id'] in seen:raise ValueError('duplicate component in catalog inputs')
            seen.add(record['component_id']);records.append(record)
            if len(records)%250==0:print('VERIFIED',len(records),'/',len(items),flush=True)
    records.sort(key=lambda e:(e['degree'],e['genus'],e['quotient_state_count'],e['component_id']))
    result=clean_catalog(dict(positive_genus_smooth_components=records))
    args.out.parent.mkdir(parents=True,exist_ok=True)
    temp=args.out.with_name(args.out.name+'.tmp')
    with temp.open('x') as stream:json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
    if json.loads(temp.read_text())!=result:raise ValueError('JSON round trip failed')
    temp.replace(args.out)
    print('EXPORT',str(args.out.resolve()),'COMPONENTS',result['component_count'],'BYTES',args.out.stat().st_size,'SHA256',sha(args.out),flush=True)


if __name__=='__main__':main()
