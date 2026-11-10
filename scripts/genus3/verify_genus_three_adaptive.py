"""Independent provenance for negative-witness smooth-prefix templates.

Templates are proposal inputs only. This verifier never certifies a positive
component, classifies isotriviality, or admits a state into a skip bank.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import csv
from pathlib import Path
import time

from verify_smooth_geometry import parse_tuple, hurwitz, compose, cycles, orbits
from run_s18_linked_roots import state_key
from run_genus_three import sha


def smooth(t):
    n=len(t[0]);ix=[n-len(cycles(p))for p in t]
    return n-len(cycles(compose(t[0],t[1])))in(ix[0]+ix[1],ix[2]+ix[3])


def verify_adaptive(folder,deadline_epoch=None):
    folder=Path(folder)
    if not(folder/'ADAPTIVE.tsv').exists():return None
    with(folder/'COUNTS.tsv').open()as stream:counts=next(csv.DictReader(stream,delimiter='\t'))
    with(folder/'ADAPTIVE_COUNTS.tsv').open()as stream:budget=next(csv.DictReader(stream,delimiter='\t'))
    n=int(counts['degree']);capacity=int(budget['capacity']);draws=int(budget['adaptive_frame_draws'])
    if counts['lane']!='refactor' or not 1<=capacity<=1024 or draws<0:raise ValueError('adaptive capacity/lane')
    with(folder/'ADAPTIVE.tsv').open()as stream:records=list(csv.DictReader(stream,delimiter='\t'))
    if not len(records)==int(budget['admitted'])<=capacity:raise ValueError('adaptive mass')
    with(folder/'RESOLUTIONS.tsv').open()as stream:roots={row['candidate']:row for row in csv.DictReader(stream,delimiter='\t')}
    seen=set();keys=set()
    for row in records:
        if deadline_epoch is not None and time.time()>=deadline_epoch:raise TimeoutError('absolute research deadline reached during adaptive audit')
        candidate=row['candidate'];prefix=int(row['prefix'])
        if candidate in seen or candidate not in roots:raise ValueError('adaptive candidate binding')
        seen.add(candidate);root=roots[candidate]
        if root['status']!='nonsmooth':raise ValueError('adaptive source is not a negative witness')
        word=list(map(int,root['word'].split(',')))
        if not 0<prefix<len(word):raise ValueError('adaptive prefix length')
        t=parse_tuple(root['seed'],n)
        if not smooth(t):raise ValueError('adaptive root not initially smooth')
        for move in word[:prefix]:
            t=hurwitz(t,move)
            if not smooth(t):raise ValueError('adaptive prefix passes an earlier nonsmooth state')
        if t!=parse_tuple(row['seed'],n):raise ValueError('adaptive prefix tuple mismatch')
        if smooth(hurwitz(t,word[prefix])):raise ValueError('adaptive prefix is not immediately before first failure')
        ident=tuple(range(n))
        if ident in t or compose(compose(t[0],t[1]),compose(t[2],t[3]))!=ident or len(orbits(t))!=1:
            raise ValueError('adaptive tuple scope')
        if sum(n-len(cycles(p))for p in t)!=2*n+4:raise ValueError('adaptive source genus')
        key=state_key(t)
        if key in keys:raise ValueError('duplicate adaptive template state')
        keys.add(key)
    return dict(passed=True,templates=len(records),adaptive_frame_draws=draws,capacity=capacity,
                sha256={name:sha(folder/name)for name in('ADAPTIVE.tsv','ADAPTIVE_COUNTS.tsv')},
                qualification='literal negative-witness smooth-prefix provenance, NOT positive closure or a skip bank')
