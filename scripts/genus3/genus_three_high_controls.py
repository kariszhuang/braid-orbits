"""Known localized-unramified positive-control constructions above degree25.

These are controls from the established isotrivial family, not new discoveries.
Independent literal permutation arithmetic; no production constructor imports.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
from audit_mathematical_patterns import cycles,compose,orbit_sizes


def inverse(p):
    q=[0]*len(p)
    for x,y in enumerate(p):q[y]=x
    return tuple(q)


def notation(p):
    return ''.join('('+' '.join(str(x+1) for x in c)+')' for c in cycles(p) if len(c)>1) or '()'


def controls():
    for r,k in [(3,9),(3,10),(3,11),(4,7),(4,8),(3,12),(3,16),(3,20),(3,21),(4,13)]:
        ident=tuple(range(k));u=tuple(range(1,k))+(0,);v=list(ident)
        v[0],v[1]=v[1],v[0];v[3],v[4]=v[4],v[3];v=tuple(v)
        if r==3:voltage=[u,v,inverse(compose(v,u))];order=[0,2,1]
        else:voltage=[compose(compose(v,inverse(u)),inverse(v)),inverse(v),u,v];order=[0,2]
        returned=ident
        for j in order:returned=compose(voltage[j],returned)
        local=inverse(returned);n=r*k
        a=local+tuple(range(k,n))
        b=tuple(((j+1)%r)*k+x for j in range(r) for x in range(k))
        c=tuple(((j+1)%r)*k+voltage[j][x] for j in range(r) for x in range(k))
        d=inverse(compose(compose(a,b),c));t=[a,b,c,d]
        assert orbit_sizes(t)==[n]
        assert sum(n-len(cycles(p)) for p in t)==2*n+4
        yield dict(degree=n,genus=3,tuple=t,text='['+', '.join(map(notation,t))+']',
                   reason='localized unramified elliptic lift; isotrivial positive control, not discovery')
