"""Exact initial-twist clique markings and middle-square necessary tests.

This constructs a nonuniform subset of right maps; it is NOT a full census.
Word outputs use native order Q1,Q2,Q3,Q1^-1,Q2^-1,Q3^-1.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import genus_three_map_core as maps


def cliques(b,e,row,k):
    length=len(row)
    if length<k:return
    forbidden=set()
    for p in(b,e):
        colour=maps.labels(p)
        for i in range(length):
            for j in range(i+1,length):
                if colour[row[i]]==colour[row[j]]:
                    forbidden.update(((i-j)%length,(j-i)%length))
    allowed=set(range(1,length))-forbidden
    def extend(chosen,candidates):
        if len(chosen)==k:
            # One representative of the first-pair square-twist support orbit.
            canonical=min(tuple(sorted((x-s)%length for x in chosen))for s in chosen)
            if tuple(chosen)==canonical:yield tuple(row[x]for x in chosen)
            return
        while candidates:
            x=min(candidates);candidates.remove(x)
            following={y for y in candidates if(y-x)%length in allowed}
            yield from extend(chosen+[x],following)
    yield from extend([0],set(allowed))


def all_markings(b,e,k):
    n=len(b);p=maps.inverse(maps.compose(b,e))
    for row in maps.cycles(p):
        for support in cliques(b,e,row,k):
            a=list(range(n))
            for x,y in zip(support,support[1:]+support[:1]):a[x]=y
            a=tuple(a);c=maps.compose(maps.inverse(a),p)
            yield(a,c,b,e)


def middle_square_failures(t):
    """Indices forbid complement rescue: return literal +/- square failures."""
    a,c,b,e=t;n=len(a)
    if maps.index(b)+maps.index(e)<=n-1:return False,False
    support=[x for x in range(n)if a[x]!=x]
    if len(support)!=maps.index(a)+1:raise ValueError('a SINGLE nontrivial cycle required')
    colour=maps.labels(c);ai,ei=maps.inverse(a),maps.inverse(e)
    positive=[colour[e[x]]for x in support]
    negative=[colour[ai[ei[x]]]for x in support]
    return len(set(positive))<len(support),len(set(negative))<len(support)


def middle_negative_word(t):
    positive,negative=middle_square_failures(t)
    return(1,1)if positive else(4,4)if negative else None
