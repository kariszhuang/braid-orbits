"""Permutation and braid-orbit helpers shared by projects/geosmooth and projects/genus4_structure.

Permutations are tuples p with p[x] the image of x, 0-based. Composition is rightmost first:
compose(p, q)[x] = p[q[x]]. parse() reads 1-based cycle notation such as "(1 2 3)(4 5)".
move(t, m) applies one of the six braid moves; orbit(t) closes the braid orbit modulo
simultaneous relabelling by S_d.
"""
import re
def compose(p,q): return tuple(p[x] for x in q)
def inverse(p):
    q=[0]*len(p)
    for x,y in enumerate(p): q[y]=x
    return tuple(q)
def cycles(p):
    seen=[0]*len(p); out=[]
    for x in range(len(p)):
        if seen[x]: continue
        l=0
        while not seen[x]: seen[x]=1; x=p[x]; l+=1
        out.append(l)
    return out
def index(p): return len(p)-len(cycles(p))
def move(t,m):
    j=m%3; a,b=t[j:j+2]; u=list(t)
    u[j:j+2]=(compose(a,compose(b,inverse(a))),a) if m<3 else (b,compose(inverse(b),compose(a,b)))
    return tuple(u)
def canonical(t):
    n=len(t[0]); gens=t+tuple(inverse(p) for p in t); best=None
    for r in range(n):
        lab={r:0}; order=[r]
        for x in order:
            for p in gens:
                y=p[x]
                if y not in lab: lab[y]=len(order); order.append(y)
        key=tuple(lab[p[x]] for p in t[:3] for x in order)
        if best is None or key<best: best=key
    return best
def unpack(key,n):
    t=tuple(tuple(key[j*n:(j+1)*n]) for j in range(3))
    return t+(inverse(compose(t[0],compose(t[1],t[2]))),)
def smooth(t):
    z=index(compose(t[0],t[1])); ix=[index(p) for p in t]
    return z==ix[0]+ix[1] or z==ix[2]+ix[3]
def orbit(t,cap=10**6):
    n=len(t[0]); root=canonical(t); seen={root}; st=[root]; i=0
    while i<len(st):
        u=unpack(st[i],n); i+=1
        for m in range(6):
            k=canonical(move(u,m))
            if k not in seen:
                seen.add(k); st.append(k)
                if len(st)>cap: return None
    return [unpack(k,n) for k in st]
def parse(parts,n):
    out=[]
    for s in parts:
        p=list(range(n))
        for row in re.findall(r'\(([^)]*)\)',s):
            c=[int(x)-1 for x in row.split()]
            for a,b in zip(c,c[1:]+c[:1]): p[a]=b
        out.append(tuple(p))
    return tuple(out)
def genus(t):
    """Genus of the cover, by Riemann-Hurwitz: 2g - 2 = -2d + sum of the four indices."""
    return 1-len(t[0])+sum(index(p) for p in t)//2
