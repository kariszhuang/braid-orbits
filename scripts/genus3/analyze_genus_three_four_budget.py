"""Independent finite necessary-budget audit, not a source-degree census.

Enumerate the marked four-block rational quotient with two pure transposition
lifts followed by a full-cycle lift. No closure or geometric classification is
inferred from feasibility. Multiplication acts rightmost first.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
from itertools import permutations, product
from collections import deque
import json
from audit_mathematical_patterns import cycles, compose


def inverse(p):
    result=[0]*len(p)
    for x,y in enumerate(p):result[y]=x
    return tuple(result)


def index(p):return len(p)-len(cycles(p))


def orbit_masks(a,d,dc):
    pending=set(range(len(a)));result=[]
    while pending:
        todo=[min(pending)];seen=set(todo)
        while todo:
            x=todo.pop()
            for p in (a,d):
                if p[x] not in seen:seen.add(p[x]);todo.append(p[x])
        pending-=seen
        result.append(sum(1<<i for i,c in enumerate(dc) if c[0] in seen))
    return result


def rows(a,b,d):
    dc=cycles(d);result=set();seen=set()
    while (a,b) not in seen:
        seen.add((a,b))
        result.update(orbit_masks(a,d,dc));result.update(orbit_masks(b,d,dc))
        a,b=compose(compose(a,b),inverse(a)),a
    return sorted(result)


def full_braid_rows(t):
    """Exact finite quotient orbit with D-cycle identities transported.

    A,B are unramified, C is the full-cycle entry. Therefore C cannot occur
    in the additive source pair: every other source index is >=k, while C's
    index is >=3k. The complementary pair is forced at every quotient state.
    """
    original=cycles(t[3]);tag=tuple(next(i for i,c in enumerate(original) if x in c)for x in range(4))
    start=(tuple(t),2,3,tag);seen={start};queue=deque([(start,[])])
    witnesses={}
    while queue:
        (q,cpos,dpos,tags),word=queue.popleft()
        selected=[0,1]if cpos>=2 else [2,3]
        if dpos in selected:
            # Each forced quotient pair must itself be an incidence forest.
            a,b=(q[i]for i in selected)
            assert index(compose(a,b))==index(a)+index(b)
            dc=cycles(q[dpos]);masks=orbit_masks(a,b,dc)
            for mask in masks:
                translated=0
                for i,c in enumerate(dc):
                    if mask&(1<<i):
                        assert len({tags[x]for x in c})==1
                        translated|=1<<tags[c[0]]
                if translated and translated not in witnesses:witnesses[translated]=word
        for move in range(6):
            i=move%3;a,b=q[i:i+2];v=list(q);cp,dp=cpos,dpos;newtags=tags
            if move<3:
                v[i:i+2]=[compose(compose(a,b),inverse(a)),a]
                relabel=a if dpos==i+1 else None
            else:
                v[i:i+2]=[b,compose(compose(inverse(b),a),b)]
                relabel=inverse(b)if dpos==i else None
            if relabel is not None:
                h=[0]*4
                for x in range(4):h[relabel[x]]=tags[x]
                newtags=tuple(h)
            if cp in(i,i+1):cp=2*i+1-cp
            if dp in(i,i+1):dp=2*i+1-dp
            state=(tuple(v),cp,dp,newtags)
            if state not in seen:seen.add(state);queue.append((state,word+[move]))
    masks=sorted(witnesses)
    maximal=[m for m in masks if not any(o!=m and m&o==m for o in masks)]
    return dict(states=len(seen),maximal_constraints=maximal,witnesses=witnesses)


def maximum_budget(constraints,cap):
    # Three D cycles. Exact two-dimensional enumeration; eliminate z using
    # all remaining inequalities. Pure integer arithmetic, no LP tolerance.
    best=-1;witness=None
    for x,y in product(range(cap+1),repeat=2):
        bound=cap
        for mask in constraints:
            used=(x if mask&1 else 0)+(y if mask&2 else 0)
            if used>cap:bound=-1;break
            if mask&4:bound=min(bound,cap-used)
        if bound>=0 and x+y+bound>best:best=x+y+bound;witness=[x,y,bound]
    return best,witness


def analyze():
    full=(1,2,3,0);transpositions=[p for p in permutations(range(4)) if index(p)==1]
    marked=[]
    for a,b in product(transpositions,repeat=2):
        ab=compose(a,b)
        if index(ab)!=2:continue
        c=full;d=inverse(compose(ab,c))
        if index(d)!=1:continue
        restricted=rows(a,b,d);braid=full_braid_rows([a,b,c,d]);constraints=sorted(map(int,braid['witnesses']))
        first=None
        for k in range(3,64):
            budget,witness=maximum_budget(constraints,k-1)
            # C has one return, at most k-1 excess; hence D needs >=k+5.
            if budget>=k+5:first=dict(inner_degree=k,source_degree=4*k,
                C_excess=k-1,D_excess=k+5,maximum_D_excess=budget,witness=witness);break
        maximal=[m for m in constraints if not any(other!=m and m&other==m for other in constraints)]
        marked.append(dict(quotient=[a,b,c,d],restricted_constraints=restricted,constraints=constraints,
            maximal_constraints=maximal,full_decorated_quotient_orbit=braid,first_budget_feasible=first))
    return dict(scope='all marked r=4 quotients with indices (1,1,3,1); pure first two lifts only',
        quotient_count=len(marked),constraint_patterns=sorted({tuple(q['constraints']) for q in marked}),
        least_budget_feasible_degree=min(q['first_budget_feasible']['source_degree'] for q in marked),
        qualification='necessary integer budgets only; not existence, smoothness, nonisotriviality or a degree lower bound outside this marked construction',
        marked=marked)


if __name__=='__main__':print(json.dumps(analyze(),indent=2))
