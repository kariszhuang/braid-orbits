"""Fresh genus-three map proposals, without prescribed invariant blocks.

This is a nonuniform construction law, NOT an exhaustive kernel census.
The source quadruple starts smooth; full braid closure is still required.
"""
import random


def compose(a,b):return tuple(a[x]for x in b)


def inverse(a):
    out=[0]*len(a)
    for x,y in enumerate(a):out[y]=x
    return tuple(out)


def cycles(a):
    unseen=set(range(len(a)));out=[]
    while unseen:
        x=min(unseen);row=[]
        while x in unseen:unseen.remove(x);row.append(x);x=a[x]
        out.append(tuple(row))
    return out


def index(a):return len(a)-len(cycles(a))


def transitive(generators):
    seen={0};queue=[0]
    for x in queue:
        for p in generators:
            y=p[x]
            if y not in seen:seen.add(y);queue.append(y)
    return len(seen)==len(generators[0])


def pair_genus(b,e):
    n=len(b)
    assert len(e)==n and sorted(b)==list(range(n)) and sorted(e)==list(range(n))
    value=index(b)+index(e)+index(compose(b,e))
    assert value%2==0
    return 1-n+value//2


def validate(b,e):
    if pair_genus(b,e)!=3 or not transitive((b,e)):raise ValueError('connected genus-three map required')


def seed(n):
    steps={7:(1,1),8:(1,2),9:(3,1),12:(3,4)}
    u,v=steps[n];b=tuple((x+u)%n for x in range(n));e=tuple((x+v)%n for x in range(n))
    validate(b,e);return b,e


def insert_after(p,corner):
    n=len(p);out=list(p)+[p[corner]];out[corner]=n;return tuple(out)


def leaf(b,e,colour,corner):
    n=len(b)
    if colour==0:return tuple(b)+(n,),insert_after(e,corner)
    if colour==1:return insert_after(b,corner),tuple(e)+(n,)
    raise ValueError('leaf colour')


def split_face(b,e,u,k):
    p=compose(b,e)
    if not any(u in row and k in row for row in cycles(p)):raise ValueError('same-face corners required')
    return insert_after(b,e[u]),insert_after(e,k)


def equal_cycle(p,x,y):return any(x in row and y in row for row in cycles(p))


def mix(b,e,x,y):
    """b*transposition, transposition*e, keeping product and genus exactly."""
    if x==y or equal_cycle(b,x,y)==equal_cycle(e,x,y):return None
    tau=list(range(len(b)));tau[x],tau[y]=y,x
    bb,ee=compose(b,tau),compose(tau,e)
    if not transitive((bb,ee)):return None
    assert compose(bb,ee)==compose(b,e) and pair_genus(bb,ee)==3
    return bb,ee


def grow(n,rng,face_split_probability=.5,mixes=2):
    """Leaf/split and exact genus-preserving mixes; no coverage assertion."""
    if n<7:raise ValueError('this seed family begins in degree7')
    b,e=seed(rng.choice([d for d in(7,8,9,12)if d<=n]))
    while len(b)<n:
        if rng.random()<face_split_probability:
            face=rng.choice(cycles(compose(b,e)));b,e=split_face(b,e,rng.choice(face),rng.choice(face))
        else:b,e=leaf(b,e,rng.randrange(2),rng.randrange(len(b)))
        for _ in range(mixes):
            x,y=rng.sample(range(len(b)),2);changed=mix(b,e,x,y)
            if changed is not None:b,e=changed
    validate(b,e);return b,e


def labels(p):
    out=[-1]*len(p)
    for j,row in enumerate(cycles(p)):
        for x in row:out[x]=j
    return out


def initial_rainbow_support(b,e,cycle,k,rng,tries=32):
    """A proposal restriction: every product shift hits distinct b/e cycles.

    In a STAR target this is necessary; outside STAR it is merely a bias.
    This function never marks a failed random clique attempt as exclusion.
    """
    length=len(cycle)
    if length<k:return None
    forbidden=set()
    for p in(b,e):
        colour=labels(p)
        for i in range(length):
            for j in range(i+1,length):
                if colour[cycle[i]]==colour[cycle[j]]:
                    forbidden.add((i-j)%length);forbidden.add((j-i)%length)
    allowed=set(range(1,length))-forbidden
    for _ in range(tries):
        positions=[0];available=set(allowed)
        for _ in range(k-1):
            if not available:break
            x=rng.choice(sorted(available));positions.append(x)
            available={y for y in available if y!=x and (y-x)%length in allowed}
        if len(positions)==k:
            positions.sort();return tuple(cycle[x]for x in positions)
    return None


def quadruple(b,e,k,rng,rainbow=True):
    """Split ONE cycle of p into k cycles: ind(a*c)=ind(a)+ind(c)."""
    n=len(b);p=inverse(compose(b,e));eligible=[row for row in cycles(p)if len(row)>=k]
    rng.shuffle(eligible)
    for row in eligible:
        support=initial_rainbow_support(b,e,row,k,rng)if rainbow else tuple(row[x]for x in sorted(rng.sample(range(len(row)),k)))
        if support is None:continue
        a=list(range(n))
        for x,y in zip(support,support[1:]+support[:1]):a[x]=y
        a=tuple(a);c=compose(inverse(a),p);t=(a,c,b,e)
        assert index(compose(a,c))==index(a)+index(c)
        assert compose(compose(a,c),compose(b,e))==tuple(range(n))
        assert sum(map(index,t))==2*n+4 and transitive(t)
        assert len(cycles(c))==len(cycles(p))+k-1
        return t
    return None


def braid(t,code):
    i,negative=divmod(code,2);a,b=t[i:i+2]
    pair=(b,compose(compose(inverse(b),a),b))if negative else(compose(compose(a,b),inverse(a)),a)
    return t[:i]+pair+t[i+2:]


def smooth(t):
    joint=index(compose(t[0],t[1]))
    return joint in(index(t[0])+index(t[1]),index(t[2])+index(t[3]))


def canonical(t):
    """Exact simultaneous-S_n key for a transitive tuple, no group library."""
    best=None;n=len(t[0])
    for root in range(n):
        lab={root:0};queue=[root]
        for x in queue:
            for p in t:
                y=p[x]
                if y not in lab:lab[y]=len(lab);queue.append(y)
        assert len(queue)==n
        key=bytes(lab[p[x]]for p in t for x in queue)
        if best is None or key<best:best=key
    return best
