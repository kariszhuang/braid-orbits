"""Independent group/block profile and deduplicated catalog for completed pilots.

Never edits the prior RESULTS file or open TeX. Geometry is not classified here.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import sympy
from sympy.combinatorics import Permutation,PermutationGroup

from audit_guided_native import known_keys
from odd_order_spin import parse_seed,cycle_lengths
from run_smooth_degree9 import sha
from run_smooth_noncyclic import recover_blocks

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--validation',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    args=p.parse_args();source=args.validation.resolve();result=json.loads((source/'RESULTS.json').read_text())
    root=args.out.resolve();root.mkdir();catalog=[];total_states=0
    # The initial bank was independently compared to the complete 66 catalog.
    known=root/'known_states.tsv';known.write_bytes((source/'known_states.tsv').read_bytes())
    for item in result['distinct_new_components']:
        folder=Path(item['job'])/item['folder'];summary=item['summary'];n=item['degree'];t=parse_seed(summary['seed'],n)
        support=[{x for x in range(n) if p[x]!=x} for p in t]
        # Three transversal entries have identical type r^k; the localized
        # support is strictly smaller, so its position is unambiguous here.
        candidates=[]
        for r in range(3,n//2+1,2):
            if n%r or n//r<3:continue
            k=n//r
            local=[i for i,s in enumerate(support) if 0<len(s)<=k]
            if len(local)!=1:continue
            if not all(sorted(cycle_lengths(t[i]))==[r]*k for i in range(4) if i!=local[0]):continue
            blocks=recover_blocks(t,support[local[0]],r,k)
            membership={x:j for j,b in enumerate(blocks) for x in b}
            if len(membership)!=n:raise ValueError('incomplete block partition')
            for q in t:
                if any(len({membership[q[x]] for x in b})!=1 for b in blocks):raise ValueError('noninvariant blocks')
            candidates.append((r,k,blocks))
        if len(candidates)!=1:raise ValueError('ambiguous/missing transversal shape')
        r,k,blocks=candidates[0]
        order=int(PermutationGroup([Permutation(list(p)) for p in t]).order())
        if order!=int(summary['group_order']):raise ValueError('independent group order mismatch')
        entry=dict(degree=n,genus=int(summary['genus']),states=int(summary['states']),
                   representative=summary['seed'],passport=[sorted(cycle_lengths(p),reverse=True) for p in t],
                   group_order=order,full_group=False,imprimitive=True,blocks=r,block_size=k,
                   block_partition=[[x+1 for x in b] for b in blocks],
                   defect=0,smooth_everywhere=True,raw_folder=str(folder),
                   summary_sha256=sha(folder/'SUMMARY.tsv'),table_sha256=sha(folder/'STATES.tsv'),
                   independent_closure_check_sha256=sha(Path(item['job'])/'CHECKS.json'))
        catalog.append(entry);total_states+=entry['states']
        with known.open('a') as output, (folder/'STATES.tsv').open() as stream:
            for row in csv.DictReader(stream,delimiter='\t'):output.write(str(n)+'\t'+row['tuple']+'\n')
    if len(known_keys(known))!=25192+total_states:raise ValueError('new-state union overlaps old catalog')
    output=dict(validation_root=str(source),validation_sha256=sha(source/'RESULTS.json'),
                components=catalog,new_components=len(catalog),new_states=total_states,
                by_degree=dict(Counter(c['degree'] for c in catalog)),by_genus=dict(Counter(c['genus'] for c in catalog)),
                total_components=66+len(catalog),total_states=25192+total_states,
                independent_group_orders=True,sympy_version=sympy.__version__,
                geometry='not classified in this mission; existing theorem may apply, human review pending',
                known_bank_sha256=sha(known))
    (root/'CATALOG.json').write_text(json.dumps(output,indent=2)+'\n')
    # A new input-binding record, not an overwrite of the original measurement.
    launch=dict(result,known_bank_sha256=sha(known),catalog_sha256=sha(root/'CATALOG.json'),
                expanded_known_components=output['total_components'],expanded_known_states=output['total_states'])
    (root/'RESULTS.json').write_text(json.dumps(launch,indent=2)+'\n')
    (root/'controls.tsv').write_bytes((source/'controls.tsv').read_bytes())
    print(json.dumps({k:v for k,v in output.items() if k!='components'},indent=2))

if __name__=='__main__':main()
