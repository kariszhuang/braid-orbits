"""Produce a compact mathematical catalog for positive-genus smooth components."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
from fractions import Fraction
from functools import lru_cache
import json
from pathlib import Path

from run_s18_linked_roots import state_key
from verify_smooth_geometry import parse_tuple, compose, cycles, orbits


def public_component(entry):
    """Check the displayed representative and normalize the approved public fields."""
    result=dict(entry)
    if not isinstance(result['primitive'],bool):raise ValueError('primitive must be a boolean')
    if result['group_order'] is None:
        if result.get('group_order_status')!='not_computed':
            raise ValueError('missing group order requires explicit not_computed status')
    elif not isinstance(result['group_order'],int) or not 0<result['group_order']<=2**53-1:
        raise ValueError('group order exceeds exact JSON numeric range')
    n=result['degree']; t=parse_tuple('['+', '.join(result['tuple'])+']',n)
    identity=tuple(range(n));product=identity
    for p in t:product=compose(product,p)
    if product!=identity or any(p==identity for p in t) or len(orbits(t))!=1:
        raise ValueError('requires identity-free transitive product-one representative')
    types=[sorted(map(len,cycles(p)),reverse=True) for p in t]
    indices=[n-len(c) for c in types]
    if result['cycle_types']!=types or result['indices']!=indices or sum(indices)!=2*(result['genus']+n-1):
        raise ValueError('representative passport/genus mismatch')
    p=compose(t[0],t[1]);cs=cycles(p);index=n-len(cs)
    if result['genus']<=0 or index not in [sum(indices[:2]),sum(indices[2:])]:
        raise ValueError('requires positive-genus smooth representative')
    if result['defect_distribution']!=[dict(defect=0,states=result['states'])]:
        raise ValueError('requires certified all-zero component distribution')
    if result.get('all_states_smooth',True) is not True:raise ValueError('not everywhere smooth')
    pair_sizes={f'g{i+1}_g{i+2}':sorted(map(len,orbits(t[i:i+2])),reverse=True) for i in range(3)}
    if pair_sizes!=result['representative_pair_orbit_sizes']:raise ValueError('pair orbit mismatch')
    if len(cs)-len(pair_sizes['g1_g2'])-len(pair_sizes['g3_g4'])+1!=0:raise ValueError('defect mismatch')
    result.pop('structure_type',None);result.pop('representative_pair_product_indices',None)
    result['contained_in']=('A' if all(i%2==0 for i in indices) else 'S')+str(n)
    result['all_states_smooth']=True
    result['representative_g1_g2']=dict(
        permutation=''.join('('+' '.join(str(x+1) for x in c)+')' for c in cs if len(c)>1) or '()',
        index=index,cycle_count=len(cs))
    for b in result['block_quotients']:
        partition=b['blocks'];m=b['block_size'];r=b['block_count']
        if len(partition)!=r or any(len(c)!=m for c in partition) or sorted(x for c in partition for x in c)!=list(range(1,n+1)):
            raise ValueError('invalid block partition')
        labels={x-1:i for i,c in enumerate(partition) for x in c}
        q=parse_tuple('['+', '.join(b['tuple'])+']',r)
        for p,v in zip(t,q):
            if any(labels[p[x]]!=v[labels[x]] for x in range(n)):raise ValueError('quotient action mismatch')
        qi=[r-len(cycles(v)) for v in q]
        excess=[a-m*c for a,c in zip(indices,qi)]
        if qi!=b['indices'] or sum(qi)!=2*(b['genus']+r-1) or excess!=b['ramification_excess'] or min(excess)<0:
            raise ValueError('quotient genus/ramification mismatch')
        orders=[result['group_order'],b['group_order'],b['kernel_order']]
        if all(v is not None for v in orders):
            if orders[0]!=orders[1]*orders[2]:raise ValueError('kernel order mismatch')
        elif b.get('orders_status')!='not_computed':
            raise ValueError('missing quotient/kernel order requires explicit status')
    if result['quotient_genera']!=sorted({b['genus'] for b in result['block_quotients']}):raise ValueError('quotient genera mismatch')
    method=result.get('isotriviality',{}).get('method')
    if method=='full-orbit-hodge-degree-v1':
        h=result.get('hodge_degree',{})
        if result['genus']<2 or h.get('normalization')!='12L=d-sum(rho(entries))+3*mean(rho(g1*g2))':
            raise ValueError('invalid full-orbit Hodge scope')
        distribution=h.get('pair_product_cycle_distribution',[])
        total=0;weighted=Fraction(0)
        for row in distribution:
            lens=row['cycle_lengths'];count=row['states']
            if sum(lens)!=n or not lens or any(type(x)is not int or x<1 for x in lens) or type(count)is not int or count<=0:
                raise ValueError('invalid pair-product cycle distribution')
            total+=count;weighted+=count*sum((Fraction(1,x)for x in lens),Fraction(0))
        if total!=result['states']:raise ValueError('Hodge distribution state mass')
        rho=sum((Fraction(1,l)for c in types for l in c),Fraction(0))
        l12=n-rho+3*weighted/total
        if l12<0 or Fraction(h['twelve_L'])!=l12 or Fraction(h['L'])!=l12/12:
            raise ValueError('inconsistent exact Hodge degree')
        expected=l12==0
        if result['isotrivial'] is not expected or result['isotriviality']['status']!=('isotrivial'if expected else'non-isotrivial'):
            raise ValueError('Hodge sign classification mismatch')
    else:
        from classify_smooth_isotrivial import classify
        result.update(classify(result))
    return result


@lru_cache(maxsize=None)
def quotient_order(key):
    from sympy.combinatorics import Permutation,PermutationGroup
    return int(PermutationGroup([Permutation(list(p)) for p in key]).order())


def clean_component(entry):
    if entry['genus']<=0 or not entry['all_states_smooth']:raise ValueError('only positive-genus everywhere-smooth components')
    lattice=entry['block_lattice']
    if not lattice['complete']:raise ValueError('incomplete primitivity analysis')
    group_order=int(entry['group']['order_decimal'])
    if group_order>2**53-1:raise ValueError('group order exceeds exact JSON/JavaScript numeric range')
    blocks=[]
    for i,s in enumerate(lattice['proper_systems']):
        q=tuple(tuple(x-1 for x in p['images_1_based']) for p in s['quotient_permutations'])
        order=quotient_order(state_key(q));kernel,remainder=divmod(group_order,order)
        if remainder:raise ValueError('quotient order does not divide source group order')
        blocks.append(dict(block_count=s['blocks'],block_size=s['block_size'],blocks=s['partition'],
            tuple=[p['cycle_notation'] for p in s['quotient_permutations']],genus=s['genus'],group_order=order,
            primitive=i in lattice['primitive_quotient_indices'],kernel_order=kernel,
            indices=s['indices'],ramification_excess=s['ramification_excess']))
        if s['block_size']==2:
            if kernel&(kernel-1):raise ValueError('two-point block kernel must be elementary abelian')
            blocks[-1]['kernel_structure']='C2^'+str(kernel.bit_length()-1) if kernel>1 else 'trivial'
    blocks.sort(key=lambda b:(b['block_size'],b['block_count']))
    pairs=entry['representative']['adjacent_pairs']
    return public_component(dict(degree=entry['degree'],genus=entry['genus'],states=entry['quotient_state_count'],
        primitive=lattice['primitive'],group_order=group_order,
        **({'id':entry['component_id']} if 'component_id' in entry else {}),
        tuple=[p['cycle_notation'] for p in entry['representative']['permutations']],
        cycle_types=entry['representative']['passport_ordered'],indices=entry['representative']['indices'],
        defect_distribution=entry['defect_distribution'],
        representative_pair_orbit_sizes={f'g{i+1}_g{i+2}':p['orbit_sizes'] for i,p in enumerate(pairs)},
        representative_pair_product_indices={f'g{i+1}_g{i+2}':p['product_index'] for i,p in enumerate(pairs)},
        quotient_genera=sorted({b['genus'] for b in blocks}),block_sizes=sorted({b['block_size'] for b in blocks}),
        block_quotients=blocks))


def clean_catalog(source):
    entries=source['positive_genus_smooth_components']
    if len({e['component_id'] for e in entries})!=len(entries):raise ValueError('duplicate components')
    result=[]
    for i,e in enumerate(entries,1):
        result.append(clean_component(e))
        if i%1000==0:print('CLEANED',i,flush=True)
    return dict(component_count=len(result),components=result)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--file',type=Path,required=True)
    args=p.parse_args();path=args.file.resolve();source=json.loads(path.read_text());result=clean_catalog(source)
    temp=path.with_name(path.name+'.tmp')
    if temp.exists():raise ValueError('temporary output already exists')
    with temp.open('x') as stream:json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
    check=json.loads(temp.read_text())
    if check!=result:raise ValueError('JSON round trip failed; original retained')
    temp.replace(path)
    print('SAVED',path,'COMPONENTS',len(result['components']),'BYTES',path.stat().st_size,
          'QUOTIENT_ORDER_CACHE',quotient_order.cache_info(),flush=True)


if __name__=='__main__':main()
