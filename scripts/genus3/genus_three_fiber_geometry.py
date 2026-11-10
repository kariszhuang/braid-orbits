"""Literal hypotheses of the fixed-linear-class full-fibre criterion.

This checks finite arithmetic, not the analytic proof or all-smoothness. See
evidence/genus_three_fiber_20261022/README.md for the separate proof argument.
No failure of this sufficient criterion is a non-isotriviality verdict.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
from audit_mathematical_patterns import quotient_check, cycles, compose, orbit_sizes


def translation_decks(q):
    """Deck transformations of a genus-one quotient that have no fixed point.

    A nonidentity automorphism of an elliptic curve is a translation iff it
    has no fixed point. A cover deck transformation can fix points only over
    branch values; there it fixes a point exactly when it preserves its inertia
    cycle setwise. This is literal finite monodromy arithmetic, not moduli proof.
    """
    r=len(q[0]);result=[]
    if orbit_sizes(q)!=[r]:raise ValueError('quotient action must be transitive')
    for target in range(1,r):
        image={0:target};used={target};todo=[0];valid=True
        while todo and valid:
            x=todo.pop();y=image[x]
            for p in q:
                a,b=p[x],p[y]
                if a in image:
                    if image[a]!=b:valid=False;break
                else:
                    if b in used:valid=False;break
                    image[a]=b;used.add(b);todo.append(a)
        if not valid:continue
        h=tuple(image[x] for x in range(r))
        fixed_over_branch=any(set(h[x] for x in cycle)==set(cycle)
                              for p in q for cycle in cycles(p))
        if not fixed_over_branch:result.append(h)
    return result


def partial_fiber_witness(t,blocks,special):
    """Sufficient argument: >=2 moving points and no quotient translation deck.

    Failure is unknown, never an isotriviality or nonsmoothness certificate.
    Independent full braid closure and separate proof review remain mandatory.
    """
    n=len(t[0]);ident=tuple(range(n))
    if not 0<=special<4 or ident in t:return None
    if compose(compose(compose(t[0],t[1]),t[2]),t[3])!=ident or orbit_sizes(t)!=[n]:return None
    if sum(n-len(cycles(p)) for p in t)!=2*n+4:return None
    q,_,excess,h,_,_=quotient_check(t,blocks)
    if h!=1 or q[special]!=tuple(range(len(blocks))):return None
    active=sum(any(t[special][x-1]!=x-1 for x in block) for block in blocks)
    if active<2 or any(e for i,e in enumerate(excess) if i!=special):return None
    translations=translation_decks(q)
    if translations:return None
    return dict(criterion='fixed-elliptic-partial-fibre-no-translation-decks',special_entry=special+1,
                quotient_genus=h,relative_excess=excess,active_blocks=active,partition=blocks,
                nonidentity_translation_decks=translations,arithmetic_checked=True,
                geometric_argument='separate proof review; no smoothness assertion')


def full_fiber_witness(t, blocks, special):
    n = len(t[0]); ident = tuple(range(n))
    if not 0 <= special < 4 or ident in t:
        return None
    if compose(compose(compose(t[0], t[1]), t[2]), t[3]) != ident or orbit_sizes(t) != [n]:
        return None
    if sum(n-len(cycles(p)) for p in t) != 2*n+4:
        return None
    q, _, excess, h, _, _ = quotient_check(t, blocks)
    if h != 1 or q[special] != tuple(range(len(blocks))):
        return None
    active = [any(t[special][x-1] != x-1 for x in block) for block in blocks]
    if not all(active) or any(e for i, e in enumerate(excess) if i != special):
        return None
    if excess[special] != 4:
        raise ValueError('relative Riemann-Hurwitz arithmetic does not reconcile')
    return dict(criterion='fixed-elliptic-full-fibre-linear-class', special_entry=special+1,
        quotient_genus=1, relative_excess=excess, active_blocks=len(blocks), partition=blocks,
        arithmetic_checked=True, geometric_argument='separate proof review; not implied by smoothness')
