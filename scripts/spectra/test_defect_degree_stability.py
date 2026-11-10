"""Literal Python audit + small exact population controls for128-point probe."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import csv
import itertools
import json
import subprocess
import tempfile
from collections import Counter
from pathlib import Path
from defect_degree_stability import (audit_tuples,analyze,coins,cycles,genus,mul,inv,orbits,braid,prob,tv,write_jobs)


def exact(d,g):
    ps=list(itertools.permutations(range(d)));identity=ps[0];hist=Counter();count=0
    for a,b,c in itertools.product(ps,repeat=3):
        t=(a,b,c,inv(mul(mul(a,b),c)))
        if identity in t or orbits(t[:3])!=1 or genus(t)!=g:continue
        count+=1;cc=cycles(mul(a,b));l=orbits(t[:2]);r=orbits(t[2:]);delta=cc-l-r+1;assert 0<=delta<=g
        hist[delta]+=1
        for m in range(6):assert braid(braid(t,m),m^1)==t and genus(braid(t,m))==g
    return hist,count


def main():
    p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);args=p.parse_args()
    for g in range(1,20):assert abs(sum(coins(g)[::2])-.5)<1e-12
    assert abs(tv([1,0,0],coins(2))-2/3)<1e-12
    assert tv([1],[0,.5,.5])==1 # Missing model tails must count.
    # Full-group parity family: a full cycle and b=a*(012).
    for d in range(3,129):
        a=tuple(list(range(1,d))+[0]);c=list(range(d));c[:3]=[1,2,0];b=mul(a,c)
        assert cycles(a)==cycles(b)==1 and mul(inv(a),b)==tuple(c)
        t=(a,inv(a),b,inv(b));assert genus(t)==d-1 and orbits(t[:3])==1
        for m in range(6):
            u=braid(t,m);assert all(cycles(p)==1 for p in u)
            assert (cycles(mul(*u[:2]))-1)%2==(d-1)%2
        if d<=6:
            identity=tuple(range(d));seen={identity};todo=[identity]
            while todo:
                q=todo.pop()
                for v in (a,b):
                    u=mul(q,v)
                    if u not in seen:seen.add(u);todo.append(u)
            import math
            assert len(seen)==math.factorial(d)//(2 if d%2 else 1)
    with tempfile.TemporaryDirectory(prefix='hurwitz-degree-control-') as raw:
        out=Path(raw);jobs=[]
        for g in (1,2):
            for rep in range(4):jobs.append(['conditional',4,g,'S',rep,89001+rep*37+g*101,1250,50000,12])
        for d in (60,90,120):jobs.append(['smooth_control',d,2,'special',0,89000+d,20,50,8])
        for d in (60,61,80):jobs.append(['parity_control',d,d-1,'special',0,91000+d,20,50,8])
        write_jobs(out/'JOBS.tsv',jobs)
        with (out/'HISTOGRAMS.tsv').open('w') as f:subprocess.run([str(args.binary),str(out/'JOBS.tsv'),str(out/'TUPLES.tsv'),str(out/'METRICS.tsv')],stdout=f,stderr=subprocess.DEVNULL,check=True)
        checked=audit_tuples(out,jobs);result=analyze(out)
        distances={}
        for g in (1,2):
            h,n=exact(4,g);c=next(c for c in result['cells'] if c['kind']=='conditional' and c['genus']==g)
            distance=tv(c['observed'],prob(h,g));assert distance<.05,(g,distance)
            distances[g]=dict(exact_tuples=n,tv=distance)
        print(json.dumps(dict(controls='PASS',literal_tuple_checks=checked,small_exact_conditional=distances,
                              fullcycle_literal_degrees=[3,128],fullcycle_exact_group_orders_degrees=[3,6])))


if __name__=='__main__':main()
