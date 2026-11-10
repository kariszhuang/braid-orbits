"""Exploratory elliptic-quotient pure-braid action tests; no theorem assertion."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
from collections import Counter
import csv
import json
from pathlib import Path

from sympy.combinatorics import Permutation,PermutationGroup
from run_smooth_degree9 import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--theory',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    args=p.parse_args();data=json.loads(args.theory.read_text());entries=data['components'];results=[]
    for e in entries:
        if e['quotient_genus']!=1:continue
        path=Path(e['raw_folder'])/'STATES.tsv'
        if sha(path)!=e['table_sha256']:raise ValueError('pattern input changed')
        with path.open() as stream:adj=[list(map(int,row['edges'].split(','))) for row in csv.DictReader(stream,delimiter='\t')]
        words=[list(range(j-1,i,-1))+[i,i]+[m+3 for m in range(i+1,j)] for i in range(3) for j in range(i+1,4)]
        def walk(s,word):
            for m in word:s=adj[s][m]
            return s
        fiber={0};queue=[0]
        for s in queue:
            for word in words:
                h=walk(s,word)
                if h not in fiber:fiber.add(h);queue.append(h)
        points=sorted(fiber);labels={s:i for i,s in enumerate(points)}
        if len(points)!=e['action']['pure_braid_fiber_size']:raise ValueError('pure fiber disagreement')
        image=int(PermutationGroup([Permutation([labels[walk(s,w)] for s in points]) for w in words]).order())
        results.append(dict(degree=e['degree'],blocks=e['blocks'],block_size=e['block_size'],states=e['states'],
            pure_fiber=len(points),pure_image_order=image,quotient_order=e['quotient_order'],
            fiber_divides_blocks=e['blocks']%len(points)==0,
            pure_image_order_divides_quotient=e['quotient_order']%image==0,table_sha256=e['table_sha256']))
    summary=dict(elliptic_components=len(results),
        fiber_divisibility_violations=sum(not e['fiber_divides_blocks'] for e in results),
        image_order_divisibility_violations=sum(not e['pure_image_order_divides_quotient'] for e in results),
        by_pure_image_order=dict(Counter(e['pure_image_order'] for e in results)),
        higher_quotient_components=sum(e['quotient_genus']>1 for e in entries),
        higher_minimum_pure_fiber=min(e['action']['pure_braid_fiber_size'] for e in entries if e['quotient_genus']>1),
        theory_sha256=sha(args.theory),
        interpretation='Exploratory finite catalog pattern, NOT a universal proof or a marked-action factorization certificate. Test outside this one-local/unramified family too.')
    args.out.write_text(json.dumps(dict(summary=summary,components=results),indent=2)+'\n');print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
