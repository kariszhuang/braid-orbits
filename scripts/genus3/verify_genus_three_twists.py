#!/usr/bin/env python3
"""Constructive full-cycle twist obstructions, checked by literal braid replay.

These functions return a genuine nonsmooth word, or None (undecided). They do
not use geometric/isotriviality criteria or producer modules.
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

from audit_guided_native import audit_guided_native
from run_genus_three import write,sha
from verify_genus_three_finite import mul,profile,braid
from verify_smooth_geometry import parse_tuple


def smooth(t):
    n=len(t[0]);ix=[n-len(profile(p))for p in t]
    return n-len(profile(mul(t[0],t[1])))in(sum(ix[:2]),sum(ix[2:]))


def fixed_point_word(t):
    n=len(t[0]);ix=[n-len(profile(p))for p in t];p=mul(t[0],t[1])
    if profile(p)!=(n,)or ix[0]+ix[2]!=n-1 or ix[1]+ix[3]<=n-1:return None
    x=next((i for i in range(n)if t[0][i]==i),None)
    target=next((i for i in range(n)if t[2][i]==i),None)
    if x is None or target is None:return None
    for k in range(n):
        if x==target:
            word=(0,)*(2*k)if k<=n-k else(3,)*(2*(n-k))
            word+=(4,);u=t
            for m in word:u=braid(u,m)
            assert not smooth(u)
            return word
        x=p[x]
    raise AssertionError('invalid full cycle')


def audit_samples(folder,out,auditor):
    """Every selected seed is literally checked; not exhaustive type coverage."""
    out.mkdir();count=0;types=Counter();lengths=Counter();n=None
    with (out/'LEADS.tsv').open('w')as leads,(out/'RESOLUTIONS.tsv').open('w')as proofs:
        proofs.write('candidate\tstatus\tword\tfolder\tseed\n')
        with (folder/'COUNTS.tsv').open()as stream:source_n=int(next(csv.DictReader(stream,delimiter='\t'))['degree'])
        with (folder/'LEADS.tsv').open()as stream:
            for row in stream:
                _,seed=row.rstrip('\n').split('\t',1);t=tuple(tuple(p)for p in parse_tuple(seed,source_n));n=source_n
                word=fixed_point_word(t)
                if word is None:continue
                count+=1;types[profile(t[0])]+=1;lengths[len(word)]+=1
                leads.write(f'{count}\t{seed}\n');proofs.write(f'{count}\tnonsmooth\t'+','.join(map(str,word))+f'\t\t{seed}\n')
    assert count and n
    (out/'COUNTS.tsv').write_text('degree\tlane\tseed\tconstructed\tscope\ttransitive\tpassport\tshort_reject\tduplicate\tclosed\tunresolved\tseconds\n'+f'{n}\tsolved\t0\t{count}\t{count}\t{count}\t{count}\t0\t0\t0\t0\t0\n')
    native=audit_guided_native(out,auditor,native_components=True)
    report=dict(passed=True,records=count,special_types={str(k):v for k,v in types.items()},word_lengths=dict(lengths),native=native,
        input_path=str(folder),input_sha256=sha(folder/'LEADS.tsv'),scope='literal regression witnesses for sampled seeds,not a census; infinite index-family coverage is the separate elementary proof')
    write(out/'CHECKS.json',report);return report


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('folder',type=Path);p.add_argument('out',type=Path);p.add_argument('auditor',type=Path);a=p.parse_args()
    r=audit_samples(a.folder,a.out,a.auditor);print(json.dumps({k:r[k]for k in ('passed','records','special_types','word_lengths')}))


if __name__=='__main__':main()
