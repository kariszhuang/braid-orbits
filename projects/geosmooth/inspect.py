import sys
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../common'))
from braid import parse, compose, inverse, orbit, smooth, index
from geo_check import orbit_partition, cycle_list, local_genus
def show(t):
    n=len(t[0]); a,b,c,e=t; p=compose(a,b)
    la,na=orbit_partition([a,b],n); lb,nb=orbit_partition([c,e],n)
    pinv=inverse(p)
    A=[[x for x in range(n) if la[x]==k] for k in range(na)]
    B=[[x for x in range(n) if lb[x]==k] for k in range(nb)]
    print(' alpha orbits',[(blk,local_genus([a,b,pinv],blk)) for blk in A])
    print(' beta orbits ',[(blk,local_genus([c,e,p],blk)) for blk in B])
    print(' node cycles ',[ (cy,la[cy[0]],lb[cy[0]]) for cy in cycle_list(p)])
    print(' indices',[index(x) for x in t],'ind(g1g2)',index(p))
n=int(sys.argv[1]); t=parse(sys.argv[2:6],n)
O=orbit(t); bad=[u for u in O if not smooth(u)]
print('states',len(O),'combinatorially non-smooth states',len(bad))
show(bad[0])
