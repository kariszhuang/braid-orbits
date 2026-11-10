"""Arithmetic audit of theorem-target exclusions, not literal braid negatives."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import csv
from pathlib import Path

from verify_smooth_geometry import parse_tuple,compose,cycles
from run_genus_three import sha


def verify_theorem_exclusions(folder):
    folder=Path(folder)
    with(folder/'COUNTS.tsv').open()as stream:counts=next(csv.DictReader(stream,delimiter='\t'))
    with(folder/'TARGET_THEOREM.tsv').open()as stream:meta=next(csv.DictReader(stream,delimiter='\t'))
    with(folder/'TARGET_EXCLUSIONS.tsv').open()as stream:rows=list(csv.DictReader(stream,delimiter='\t'))
    total=int(meta['excluded_occurrences']);sampled=int(meta['sampled_occurrences']);n=int(counts['degree'])
    if meta['theorem']!='no-positive-genus-all-smooth-transposition'or not 0<=sampled==len(rows)<=min(total,16):raise ValueError('theorem exclusion metadata')
    seen=set()
    for row in rows:
        ident=int(row['candidate']);t=parse_tuple(row['seed'],n)
        if ident in seen or not 0<ident<=int(counts['constructed']):raise ValueError('theorem exclusion candidate')
        seen.add(ident);indices=[n-len(cycles(p))for p in t]
        if 0 in indices or sum(indices)!=2*n+4 or 1 not in indices:raise ValueError('not an identity-free genus-three transposition tuple')
        if compose(compose(t[0],t[1]),compose(t[2],t[3]))!=tuple(range(n)):raise ValueError('theorem exclusion product')
    return dict(passed=True,excluded_occurrences=total,sampled_occurrences=sampled,
                avoided_sorted_profiles=int(meta['avoided_sorted_profiles']),
                sources={name:sha(folder/name)for name in('TARGET_THEOREM.tsv','TARGET_EXCLUSIONS.tsv')},
                qualification='goal/theorem exclusions, NOT saved nonsmooth braid witnesses or a degree census')
