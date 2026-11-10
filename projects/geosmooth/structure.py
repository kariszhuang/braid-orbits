"""Group order and block systems of the geometric-only orbits (generic, small d)."""
import sys, itertools
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../common')); sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../genus4_structure'))
from braid import parse, compose, inverse, index
from blocks import block_systems
def closure(gens):
    n=len(gens[0]); idp=tuple(range(n)); G={idp}; fr=[idp]
    while fr:
        nf=[]
        for g in fr:
            for s in gens:
                h=compose(s,g)
                if h not in G: G.add(h); nf.append(h)
        fr=nf
    return G
for d,tup in [(6,["(1 2 3 6 5 4)","(1 3 5)","(1 3 5)","(1 4)(2 3)(5 6)"]),(6,["(1 2 5)(3 4 6)","(1 2 5)","(1 3)(2 4)(5 6)","(1 4)(2 6)(3 5)"])]:
    t=parse(tup,d); G=closure(list(t))
    print(tup,'order',len(G),'blocks',sorted(block_systems(t)))
    for bl in block_systems(t):
        m=len(bl); bof={x:i for i,b in enumerate(bl) for x in b}
        q=[tuple(bof[p[b[0]]] for b in bl) for p in t]
        hq=1-m+sum(index(p) for p in q)//2
        print('   quotient on',m,'blocks genus',hq,'quotient indices',[index(p) for p in q])
