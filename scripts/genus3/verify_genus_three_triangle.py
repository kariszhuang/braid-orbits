#!/usr/bin/env python3
"""Independent exact coverage and literal certificate check. No Rust imports."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
from collections import Counter
import csv
import json
from math import factorial
from pathlib import Path
from audit_guided_native import audit_guided_native
from verify_genus_three_finite import connection,mul,profile,pair_key
from verify_smooth_geometry import parse_tuple
from run_genus_three import sha,write

def audit(folder,auditor):
    folder=Path(folder);c=json.loads((folder/'COVERAGE.json').read_text());n=c['degree'];assert n in (5,9,13)
    invol=(1,)+(2,)*((n-1)//2);expected=connection(invol,(n,))
    assert connection(invol,invol)==n
    assert c['involutions']==factorial(n)//(2**((n-1)//2)*factorial((n-1)//2))
    native=audit_guided_native(folder,auditor,native_components=True)
    assert native['independent_native']['records']==expected
    distinct=set();words=Counter();a=tuple((-i)%n for i in range(n));b=tuple((-i-1)%n for i in range(n));h=tuple((i+1)%n for i in range(n))
    assert mul(a,b)==h
    # The pair centralizer is trivial: it centralizes h, hence is h^k, and
    # commuting with a forces it to fix a's unique fixed point0, hence k=0.
    for k in range(1,n):
        power=tuple((i+k)%n for i in range(n));assert mul(power,a)!=mul(a,power)
    with (folder/'RESOLUTIONS.tsv').open() as stream:
        for row in csv.DictReader(stream,delimiter='\t'):
            t=parse_tuple(row['seed'],n)
            assert tuple(t[0])==a and tuple(t[1])==b
            assert all(profile(tuple(p))==invol for p in t[:3]) and profile(tuple(t[3]))==(n,)
            third=tuple(t[2]);assert third not in distinct;distinct.add(third)
            assert len(row['word'].split(','))<=n
            words[len(row['word'].split(','))]+=1
    assert len(distinct)==expected==c['necessary_Sn_states']==c['negative']
    assert c['genus']==(n-1)//4 and c['unresolved']==0
    return dict(status='complete_negative_passport',coverage=c,exact_character_count=expected,
        first_pair_frames=1,centralizer_order=1,word_length_histogram=dict(words),native=native,
        scope='one entire passport,not a degree census; degrees5/9 are regressions of completed censuses')

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('folder',type=Path);p.add_argument('auditor',type=Path)
    a=p.parse_args();r=audit(a.folder,a.auditor);write(a.folder/'EXACT_AUDIT.json',r);print(json.dumps(r['coverage']))

if __name__=='__main__':main()
