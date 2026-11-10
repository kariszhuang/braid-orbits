"""Exact graph-derived passport and pure-braid fibers of audited components.

This verifies graph identities, not tuple edge correctness (already independently
checked in the complete certificate). No large permutation-group enumeration.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
from collections import Counter
import csv
from pathlib import Path

from verify_smooth_geometry import parse_tuple,cycles
from run_smooth_degree9 import sha


def action_profile(folder,expected_sha):
    path=Path(folder)/'STATES.tsv'
    if sha(path)!=expected_sha:raise ValueError('action certificate changed')
    with path.open() as stream:rows=list(csv.DictReader(stream,delimiter='\t'))
    if not rows:raise ValueError('empty action')
    n=len(parse_tuple(rows[0]['tuple'],int(rows[0].get('degree',0)))) if 'degree' in rows[0] else None
    if n is None:
        with (Path(folder)/'SUMMARY.tsv').open() as stream:n=int(next(csv.DictReader(stream,delimiter='\t'))['degree'])
    adj=[];passports=[]
    for i,row in enumerate(rows):
        if int(row['state'])!=i:raise ValueError('action row order')
        edges=list(map(int,row['edges'].split(',')))
        if len(edges)!=6 or any(not 0<=e<len(rows) for e in edges):raise ValueError('complete six-edge action required')
        adj.append(edges);t=parse_tuple(row['tuple'],n)
        passports.append(tuple(tuple(sorted(map(len,cycles(p)),reverse=True)) for p in t))
    words=[]
    for i in range(3):
        for j in range(i+1,4):words.append(list(range(j-1,i,-1))+[i,i]+[m+3 for m in range(i+1,j)])
    parent=list(range(len(rows)))
    def root(x):
        while x!=parent[x]:parent[x]=parent[parent[x]];x=parent[x]
        return x
    def join(a,b):
        a,b=root(a),root(b)
        if a!=b:parent[b]=a
    for s in range(len(rows)):
        for m in range(3):
            if adj[adj[s][m]][m+3]!=s or adj[adj[s][m+3]][m]!=s:raise ValueError('inverse graph action')
        for i in range(2):
            def walk(word):
                x=s
                for m in word:x=adj[x][m]
                return x
            if walk([i,i+1,i])!=walk([i+1,i,i+1]):raise ValueError('braid relation')
        if adj[adj[s][0]][2]!=adj[adj[s][2]][0]:raise ValueError('distant braid commute')
        for word in words:
            h=s
            for m in word:h=adj[h][m]
            if passports[s]!=passports[h]:raise ValueError('pure word changed ordered passport')
            join(s,h)
    sizes=Counter(root(s) for s in range(len(rows)))
    if len(set(sizes.values()))!=1:raise ValueError('normal pure braid fibers not equisized')
    passport_sizes=Counter(passports);pure_size=next(iter(sizes.values()))
    fibers=len(sizes)
    if 24%fibers:raise ValueError('B4/P4=S4 quotient has invalid transitive index')
    return dict(states=len(rows),ordered_passports=len(passport_sizes),
        passport_fiber_sizes=dict(Counter(passport_sizes.values())),
        pure_braid_orbits=fibers,pure_braid_fiber_size=pure_size,
        pure_braid_fixes_every_state=pure_size==1,
        graph_braid_relations_checked=True,table_sha256=expected_sha)
