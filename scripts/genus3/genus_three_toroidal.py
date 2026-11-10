"""Construct connected elliptic dessins by genus-preserving edge insertions.

This is a diverse proposal law, not a census or a claim that all maps occur.
No invariant block system of the quotient group is prescribed.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import random

from genus_three_intrinsic import cycles,compose


def insert_after(p,corner):
    q=len(p);out=list(p)+[p[corner]];out[corner]=q;return tuple(out)


def leaf(beta,gamma,colour,corner):
    q=len(beta)
    if colour=='black':return tuple(beta)+(q,),insert_after(gamma,corner)
    if colour=='white':return insert_after(beta,corner),tuple(gamma)+(q,)
    raise ValueError('leaf colour')


def split_face(beta,gamma,u,k):
    delta=compose(beta,gamma)
    if not any(u in row and k in row for row in cycles(delta)):raise ValueError('corners must share a face')
    return insert_after(beta,gamma[u]),insert_after(gamma,k)


def validate(beta,gamma):
    q=len(beta)
    if len(gamma)!=q or sorted(beta)!=list(range(q))or sorted(gamma)!=list(range(q)):raise ValueError('dessin permutation')
    if sum(len(cycles(p))for p in(beta,gamma,compose(beta,gamma)))!=q:raise ValueError('not elliptic')
    seen={0};queue=[0]
    for x in queue:
        for p in(beta,gamma):
            if p[x]not in seen:seen.add(p[x]);queue.append(p[x])
    if len(seen)!=q:raise ValueError('disconnected dessin')


def grow(q,rng,seed=None,need_two_unramified=False):
    if q<3:raise ValueError('elliptic quotient degree>=3')
    beta,gamma=seed or((1,2,0),(1,2,0))
    target=q-2 if need_two_unramified else q
    if target<len(beta):beta=gamma=(1,2,0)
    while len(beta)<target:
        counts=[len(cycles(beta)),len(cycles(gamma)),len(cycles(compose(beta,gamma)))];small=min(counts)
        operation=rng.choice([j for j,c in enumerate(counts)if c==small])if rng.random()<.8 else rng.randrange(3)
        if operation<2:beta,gamma=leaf(beta,gamma,'black'if operation==0 else'white',rng.randrange(len(beta)))
        else:
            row=rng.choice(cycles(compose(beta,gamma)));beta,gamma=split_face(beta,gamma,rng.choice(row),rng.choice(row))
    if need_two_unramified:
        while len(beta)<q:beta,gamma=leaf(beta,gamma,'white',rng.randrange(len(beta)))
    validate(beta,gamma)
    if need_two_unramified and sum(len(c)==1 for c in cycles(gamma))<2:raise ValueError('missing fixed unramified vertices')
    return beta,gamma


def marked_key(beta,gamma,mask,fixed):
    """Rooted canonical ordered quotient plus moving/fixed decorations."""
    q=len(beta);vertices={'b':cycles(beta),'g':cycles(gamma),'d':cycles(compose(beta,gamma))};best=None
    for start in range(q):
        labels={start:0};queue=[start]
        for x in queue:
            for p in(beta,gamma):
                y=p[x]
                if y not in labels:labels[y]=len(queue);queue.append(y)
        pair=bytes(labels[p[x]]for p in(beta,gamma)for x in queue)
        moving=sum(1<<labels[x]for x in range(q)if mask&(1<<x))
        decoration=sorted((kind,sum(1<<labels[x]for x in vertices[kind][j]))for kind,j in fixed)
        key=pair+moving.to_bytes(4,'little')+b''.join(kind.encode()+value.to_bytes(4,'little')for kind,value in decoration)
        if best is None or key<best:best=key
    return best
