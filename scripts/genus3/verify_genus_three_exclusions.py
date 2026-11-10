"""Independent literal witnesses for the supplied single-local elliptic criterion.

These are isotriviality-target exclusions, NEVER non-smooth braid witnesses.
No production Rust/constructor helpers are imported.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import csv
import hashlib
from pathlib import Path

from audit_mathematical_patterns import quotient_check, orbit_sizes, cycles, compose
from verify_smooth_geometry import parse_tuple


def audit_exclusions(path, expected=None):
    path=Path(path); count=0; previous={}
    with path.open() as stream:
        reader=csv.DictReader(stream,delimiter='\t')
        if reader.fieldnames!=['degree','candidate','special','labels','seed']:raise ValueError('isotrivial witness header')
        for row in reader:
            n=int(row['degree']); k=int(row['special'])-1; rank=int(row['candidate'])
            if not 0<=k<4 or rank<=previous.get(n,0):raise ValueError('isotrivial witness order')
            previous[n]=rank
            t=parse_tuple(row['seed'],n); ident=tuple(range(n))
            indices=[n-len(cycles(p)) for p in t]
            if ident in t or orbit_sizes(t)!=[n] or compose(compose(compose(t[0],t[1]),t[2]),t[3])!=ident:
                raise ValueError('isotrivial witness tuple scope')
            if sum(indices)!=2*n+4:raise ValueError('isotrivial witness genus')
            labels=list(map(int,row['labels'].split(',')))
            if len(labels)!=n or set(labels)!=set(range(max(labels)+1)):raise ValueError('isotrivial block labels')
            blocks=[[x+1 for x in range(n) if labels[x]==b] for b in range(max(labels)+1)]
            if not 1<len(blocks)<n:raise ValueError('proper blocks required')
            q,qi,excess,h,branch,_=quotient_check(t,blocks)
            active={labels[x] for x in range(n) if t[k][x]!=x}
            if q[k]!=tuple(range(len(blocks))) or len(active)!=1 or h!=1 or excess[k]!=4 or any(v for j,v in enumerate(excess) if j!=k):
                raise ValueError('Proposition 1 hypotheses do not hold')
            count+=1
    if expected is not None and count!=expected:raise ValueError('isotrivial exclusion counter mismatch')
    return dict(records=count,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),hypotheses='Proposition1(a)-(c), quotient genus1; source genus3',
        scope='static sufficient isotriviality witnesses, not a full-orbit smoothness census')


if __name__=='__main__':
    import json,sys
    print(json.dumps(audit_exclusions(Path(sys.argv[1])),indent=2))
