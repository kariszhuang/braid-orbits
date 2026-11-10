"""Fresh genus-three DOUBLE-transposition markings, without imposed blocks.

The two edges may share one ordinary cycle (a forest chain). Requiring four
different cycle colours would be incorrect. Coupled product rotations use
their actual common phase; different faces are never rotated independently.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
from itertools import combinations
from math import gcd

import genus_three_map_core as maps


def forest(edges,colours):
    """Exact incidence forest of TWO disjoint transpositions after contraction."""
    a,b=(tuple(sorted((colours[x],colours[y])))for x,y in edges)
    return a[0]!=a[1] and b[0]!=b[1] and a!=b


class DoubleScreen:
    """Reuse right-map colours and individual-edge phase traces in RAM."""
    def __init__(self,p,b,e):
        self.rows=maps.cycles(p)
        self.locations={x:(j,i)for j,row in enumerate(self.rows)for i,x in enumerate(row)}
        self.colours=(maps.labels(b),maps.labels(e));self.cache={};self.phases={}

    def trace(self,edge):
        key=tuple(sorted(edge))
        if key in self.cache:return self.cache[key]
        x,y=key;face,at=self.locations[x];other,j=self.locations[y]
        if face!=other:raise ValueError('an additive chord lies within one product face')
        row=self.rows[face];length=len(row);traces=[]
        for colour in self.colours:
            trace=[]
            for t in range(length):
                pair=tuple(sorted((colour[row[(at+t)%length]],colour[row[(j+t)%length]])))
                if pair[0]==pair[1]:self.cache[key]=None;return None
                trace.append(pair)
            traces.append(tuple(trace))
        value=tuple(traces);self.cache[key]=value;return value

    def palette(self,edge,trace,which,common):
        key=(tuple(sorted(edge)),which,common)
        if key not in self.phases:
            phases=[set()for _ in range(common)]
            for t,pair in enumerate(trace):phases[t%common].add(pair)
            self.phases[key]=phases
        return self.phases[key]

    def accepts(self,edges):
        traces=[self.trace(edge)for edge in edges]
        if any(trace is None for trace in traces):return False
        for which in(0,1):
            left,right=traces[0][which],traces[1][which]
            common=gcd(len(left),len(right))
            if len(left)==len(right):
                if any(x==y for x,y in zip(left,right)):return False
            else:
                a=self.palette(edges[0],left,which,common);b=self.palette(edges[1],right,which,common)
                # CRT compatibility, never unrestricted independent phases.
                if any(x&y for x,y in zip(a,b)):return False
        return True


def first_twists(edges,p,b,e):return DoubleScreen(p,b,e).accepts(edges)


def _one_face(row):
    length=len(row)
    for tail in combinations(range(1,length),3):
        points=(0,)+tail
        for local in(((points[0],points[1]),(points[2],points[3])),
                      ((points[0],points[3]),(points[1],points[2]))):
            key=tuple(sorted(tuple(sorted(edge))for edge in local))
            canonical=min(tuple(sorted(tuple(sorted(((x-s)%length,(y-s)%length)))for x,y in local))for s in points)
            if key==canonical:yield tuple((row[x],row[y])for x,y in local)


def _two_faces(first,second):
    length,other=len(first),len(second)
    for distance in range(1,length//2+1):
        # Only common product powers fixing the first UNORDERED edge remain.
        period=distance if 2*distance==length else length
        common=gcd(period,other)
        for x in range(common):
            for y in range(x+1,other):
                edge=(x,y)
                canonical=min(tuple(sorted(((x-(s//common)*common)%other,
                                            (y-(s//common)*common)%other)))for s in edge)
                if edge==canonical:yield((first[0],first[distance]),(second[x],second[y]))


def chord_pairs(p):
    rows=maps.cycles(p)
    for row in rows:
        if len(row)>=4:yield from _one_face(row)
    for i,first in enumerate(rows):
        if len(first)<2:continue
        for second in rows[i+1:]:
            if len(second)>=2:yield from _two_faces(first,second)


def candidates(b,e,limit=0):
    n=len(b);p=maps.inverse(maps.compose(b,e));attempted=0;screen=DoubleScreen(p,b,e)
    for edges in chord_pairs(p):
        attempted+=1
        if limit and attempted>limit:return
        if not screen.accepts(edges):continue
        a=list(range(n))
        for x,y in edges:a[x],a[y]=y,x
        a=tuple(a);c=maps.compose(a,p);t=(a,c,b,e)
        assert maps.index(a)==2
        assert maps.index(a)+maps.index(c)==maps.index(p)
        assert sum(map(maps.index,t))==2*n+4 and maps.transitive(t)
        yield t


def middle_word(t):
    a,c,b,e=t;n=len(a)
    if maps.index(b)+maps.index(e)<=n-1:return None
    edges=[row for row in maps.cycles(a)if len(row)>1]
    if len(edges)!=2 or any(len(row)!=2 for row in edges):raise ValueError('two transpositions required')
    colours=maps.labels(c)
    if not forest(tuple((e[x],e[y])for x,y in edges),colours):return(1,1)
    ai,ei=maps.inverse(a),maps.inverse(e)
    if not forest(tuple((ai[ei[x]],ai[ei[y]])for x,y in edges),colours):return(4,4)
    return None
