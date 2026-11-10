"""Read-only arithmetic audit of published/ancestral source-genus-three data.

Uses saved complete block lattices; neither changes classifications nor infers
smooth closure from the prospective geometric criteria.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
from collections import Counter
import json
from pathlib import Path
from genus_three_fiber_geometry import partial_fiber_witness,full_fiber_witness,translation_decks
from run_genus_three import sha,write
from verify_smooth_geometry import parse_tuple


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('catalog',type=Path);p.add_argument('new',type=Path)
    p.add_argument('--out',type=Path)
    a=p.parse_args();catalog=json.loads(a.catalog.read_text());new=json.loads(a.new.read_text())
    counts=Counter();genera=Counter();active=Counter();conflicts=[];bindings={}
    for item in catalog['components']:
        if item['genus']!=3:continue
        counts['published']+=1;n=item['degree'];t=parse_tuple('['+', '.join(item['tuple'])+']',n)
        for b in item['block_quotients']:
            genera[str(b['genus'])]+=1
            q=parse_tuple('['+', '.join(b['tuple'])+']',b['block_count'])
            ident=tuple(range(b['block_count']))
            for special in range(4):
                if q[special]!=ident or b['genus']!=1:continue
                number=sum(any(t[special][x-1]!=x-1 for x in block)for block in b['blocks'])
                active[str(number)]+=1
                if number>=2:counts['published_multi_active_translation_deck_systems']+=bool(translation_decks(q))
                if partial_fiber_witness(t,b['blocks'],special):conflicts.append(item['id'])
                if full_fiber_witness(t,b['blocks'],special):conflicts.append(item['id'])
    for item in new:
        n=item['degree'];t=parse_tuple(item['representative'],n);path=Path(item['raw_folder'])/'GEOMETRY_REVIEW.json'
        digest=sha(path)
        if digest!=item['geometry_review_sha256']:raise ValueError('ancestral geometry binding changed')
        bindings[item['id']]=digest;profile=json.loads(path.read_text())
        if item['genus']!=3 or item['original_isotrivial'] is not True:raise ValueError('ancestral target population changed')
        counts['unpublished']+=1
        for b in profile['block_lattice']['proper_systems']:
            for entry in b['identity_entries']:
                if partial_fiber_witness(t,b['partition'],entry-1):conflicts.append(item['id'])
                if full_fiber_witness(t,b['partition'],entry-1):conflicts.append(item['id'])
    if conflicts:raise ValueError(f'conflicting criteria: {conflicts}')
    result=dict(passed=True,counts=counts,published_block_system_genera=genera,
        published_elliptic_identity_entry_active_blocks=active,conflicts=conflicts,
        catalog_sha256=sha(a.catalog),new_components_sha256=sha(a.new),geometry_bindings=bindings,
        scope='finite arithmetic consistency check; classifications unchanged; does not prove geometric statements')
    if a.out:write(a.out,result)
    print(json.dumps({k:v for k,v in result.items()if k!='geometry_bindings'},indent=2))


if __name__=='__main__':main()
