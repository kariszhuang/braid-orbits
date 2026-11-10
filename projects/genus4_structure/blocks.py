"""Block-system / quotient analysis of the genus-4 non-isotrivial all-smooth components.
For each block system with a genus-1 quotient and a block-trivial entry s, report block size k,
quotient degree m, moving relative branch types (s on each block) and fixed types (holonomy of the
other entries over each quotient cycle)."""
import json,sys,collections
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../common'))
from braid import parse, compose, inverse, cycles, index
def block_systems(t):
    n=len(t[0]); out=set()
    for j in range(1,n):
        par=list(range(n))
        def f(x):
            while par[x]!=x: par[x]=par[par[x]]; x=par[x]
            return x
        par[f(j)]=f(0); changed=True
        while changed:
            changed=False
            cls=collections.defaultdict(list)
            for x in range(n): cls[f(x)].append(x)
            for c in cls.values():
                for p in t:
                    a=f(p[c[0]])
                    for x in c[1:]:
                        b=f(p[x])
                        if a!=b: par[b]=a; a=f(a); changed=True
        cls=collections.defaultdict(list)
        for x in range(n): cls[f(x)].append(x)
        bl=tuple(sorted(tuple(sorted(c)) for c in cls.values()))
        if 1<len(bl)<n: out.add(bl)
    return out
def analyse(t):
    n=len(t[0]); res=[]
    for bl in block_systems(t):
        m=len(bl); k=n//m; bof={x:i for i,b in enumerate(bl) for x in b}
        q=[tuple(bof[p[b[0]]] for b in bl) for p in t]
        hq=1-m+sum(index(p) for p in q)//2
        triv=[i for i,p in enumerate(q) if index(p)==0]
        rel=[]
        for i,(p,pq) in enumerate(zip(t,q)):
            seen=set(); types=[]
            for c0 in range(m):
                if c0 in seen: continue
                cyc=[]; x=c0
                while x not in seen: seen.add(x); cyc.append(x); x=pq[x]
                L=len(cyc); pts=set(y for b in cyc for y in bl[b]); lens=[]; s2=set()
                for y in pts:
                    if y in s2: continue
                    l=0; z=y
                    while z not in s2: s2.add(z); z=p[z]; l+=1
                    lens.append(l//L)
                ty=tuple(sorted([a for a in lens if a>1],reverse=True))
                if ty: types.append(ty)
            rel.append(sorted(types,reverse=True))
        res.append(dict(k=k,m=m,qgenus=hq,trivial=triv,rel=rel))
    return res
if __name__=="__main__":
    # usage: python3 blocks.py FILE.jsonl  (one record per line with "tuple" and "degree")
    rows=[json.loads(l) for l in open(sys.argv[1])]
    summ=collections.Counter()
    for r in rows:
        t=parse(r['tuple'],r['degree'])
        for a in analyse(t):
            if a['qgenus']==1 and a['trivial']:
                s=a['trivial'][0]; mov=a['rel'][s]; fix=[x for i,x in enumerate(a['rel']) if i!=s]
                summ[(r['degree'],a['k'],a['m'],str(mov),str(fix))]+=1
    for k,v in sorted(summ.items()): print(v,k)
