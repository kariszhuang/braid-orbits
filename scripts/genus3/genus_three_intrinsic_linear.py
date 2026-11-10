"""Projective160-fiber coordinates and exact inverse tree--cotree gauge.

No curve is certified all-smooth here. Matrix paths give literal necessary
tests, and surviving candidates require an independent full six-move audit.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import itertools
from functools import lru_cache

from genus_three_intrinsic import IntrinsicFrame,PERMS,RANK,MUL,INV,TRANS,cycles,compose,inverse
from verify_smooth_geometry import hurwitz

PUSH_WORDS=((0,0),(4,0,0,1))
PROJECTIVE=tuple(v for v in itertools.product(range(3),repeat=4)if any(v)and next(x for x in v if x)!=2)
assert len(PROJECTIVE)==40
PROJECTIVE_ID={v:i for i,v in enumerate(PROJECTIVE)}
SIGN_STATES=tuple(itertools.product(range(2),repeat=2))
HYPERPLANES={form:sum(1<<i for i,v in enumerate(PROJECTIVE)if sum(a*b for a,b in zip(form,v))%3==0)for form in itertools.product(range(3),repeat=4)}


def affine(sign,t):return RANK[tuple(((-1 if sign else 1)*x+t)%3 for x in range(3))]


def sign(p):return int((PERMS[p][1]-PERMS[p][0])%3==2)


def conjugate(h,p):return MUL[MUL[h][p]][INV[h]]


def apply_word(t,word):
    for m in word:t=hurwitz(t,m)
    return t


def project(v):
    first=next((x for x in v if x),0)
    if first==0:return None
    return tuple((x if first==1 else -x)%3 for x in v)


class LinearFrame:
    def __init__(self,beta,gamma):
        self.frame=IntrinsicFrame(beta,gamma);self.q=len(beta)
        self.black=cycles(beta);self.white=cycles(gamma)
        self.bl={x:j for j,c in enumerate(self.black)for x in c}
        self.wh={x:j for j,c in enumerate(self.white)for x in c}

    @lru_cache(maxsize=4096)
    def lift(self,mask,signs,vector):
        active=[j for j in range(self.q)if mask&(1<<j)]
        if len(active)!=4:raise ValueError('four active branch values required')
        u,v,c2,c3=vector;local=[0]*self.q
        local[active[0]]=affine(1,0);local[active[1]]=affine(1,c2);local[active[2]]=affine(1,c3)
        found=[]
        for c4 in range(3):
            local[active[3]]=affine(1,c4)
            t=self.frame.solve(local,[affine(signs[0],u),affine(signs[1],v)])
            if t is not None:found.append(t)
        if len(found)!=1:raise ValueError('root translation equation must have a unique fourth color')
        return found[0]

    def extract(self,t,check=True):
        q=self.q
        def sheet(p,j,target):
            if any(p[3*j+x]//3!=target for x in range(3)):raise ValueError('point push changed the quotient marking')
            return RANK[tuple(p[3*j+x]%3 for x in range(3))]
        local=[sheet(t[0],j,j)for j in range(q)]
        bv=[sheet(t[1],j,self.frame.beta[j])for j in range(q)]
        cv=[sheet(t[2],j,self.frame.gamma[j])for j in range(q)]
        # First trivialize beta lifts, separately on every beta-cycle.
        h=[None]*q
        for row in self.black:
            h[row[0]]=0
            for j in row:
                nxt=self.frame.beta[j];value=MUL[h[j]][INV[bv[j]]]
                if h[nxt]is None:h[nxt]=value
                elif h[nxt]!=value:raise ValueError('nontrivial beta return')
        local=[conjugate(h[j],local[j])for j in range(q)]
        cv=[MUL[MUL[h[self.frame.gamma[j]]][cv[j]]][INV[h[j]]]for j in range(q)]
        # Recover gamma potentials; one initial potential per white vertex.
        z=[None]*q
        for row in self.white:
            z[row[0]]=0
            for j in row:
                nxt=self.frame.gamma[j];value=MUL[z[j]][INV[cv[j]]]
                if z[nxt]is None:z[nxt]=value
                elif z[nxt]!=value:raise ValueError('nontrivial gamma return')
        # Residual black/white gauges set all primal-tree edge potentials to1.
        black=[None]*len(self.black);white=[None]*len(self.white);black[0]=0
        pending=set(self.frame.primal)
        while pending:
            progress=False
            for j in list(pending):
                b,w=self.bl[j],self.wh[j]
                if black[b]is not None:
                    value=MUL[black[b]][INV[z[j]]]
                    if white[w]is not None and white[w]!=value:raise ValueError('primal gauge inconsistency')
                    white[w]=value
                elif white[w]is not None:black[b]=MUL[white[w]][z[j]]
                else:continue
                pending.remove(j);progress=True
            if not progress:raise ValueError('disconnected primal gauge')
        local=[conjugate(black[self.bl[j]],local[j])for j in range(q)]
        z=[MUL[MUL[white[self.wh[j]]][z[j]]][INV[black[self.bl[j]]]]for j in range(q)]
        if any(z[j]!=0 for j in self.frame.primal):raise ValueError('tree potential not identity')
        active=[j for j,a in enumerate(local)if a!=0]
        if len(active)!=4 or any(local[j]not in TRANS for j in active):raise ValueError('wrong moving ramification')
        mask=sum(1<<j for j in active)
        translation=affine(0,PERMS[local[active[0]]][0])
        local=[conjugate(translation,p)for p in local];z=[conjugate(translation,p)for p in z]
        a,b=(z[j]for j in self.frame.handles)
        signs=(sign(a),sign(b));vector=(PERMS[a][0],PERMS[b][0],PERMS[local[active[1]]][0],PERMS[local[active[2]]][0])
        if check:
            normalized=self.lift(mask,signs,vector)
            gauges=[MUL[MUL[translation][black[self.bl[j]]]][h[j]]for j in range(q)]
            relabel=tuple(3*j+PERMS[gauges[j]][x]for j in range(q)for x in range(3));ri=inverse(relabel)
            actual=tuple(compose(compose(relabel,p),ri)for p in t)
            if actual!=normalized:raise ValueError('inverse chart is not literal simultaneous conjugation')
        return mask,signs,vector

    @lru_cache(maxsize=8192)
    def matrix(self,mask,signs,word):
        zero=self.extract(apply_word(self.lift(mask,signs,(0,0,0,0)),word))
        if zero[2]!=(0,0,0,0):raise ValueError('affine transition has a nonzero origin')
        columns=[]
        for j in range(4):
            v=tuple(int(i==j)for i in range(4));out=self.extract(apply_word(self.lift(mask,signs,v),word))
            if out[:2]!=zero[:2]:raise ValueError('sign/support transition depends on translation')
            columns.append(out[2])
        matrix=tuple(tuple(columns[j][i]for j in range(4))for i in range(4))
        if len(set(projective_images(matrix)))!=40:raise ValueError('singular projective transition')
        return zero[0],zero[1],matrix


def matvec(matrix,vector):return tuple(sum(a*b for a,b in zip(row,vector))%3 for row in matrix)


@lru_cache(maxsize=4096)
def projective_images(matrix):
    return tuple(PROJECTIVE_ID[project(matvec(matrix,v))]for v in PROJECTIVE)


def minimal_block_genus(beta,gamma,mask):
    """Minimal invariant equivalence identifying S; proposed geometric test."""
    q=len(beta);parent=list(range(q))
    def root(x):
        while parent[x]!=x:parent[x]=parent[parent[x]];x=parent[x]
        return x
    def join(x,y):
        x,y=root(x),root(y)
        if x==y:return False
        parent[y]=x;return True
    active=[j for j in range(q)if mask&(1<<j)]
    for x in active[1:]:join(active[0],x)
    changed=True
    while changed:
        changed=False
        for x in range(q):
            for y in range(x):
                if root(x)==root(y):
                    for p in(beta,gamma):changed=join(p[x],p[y])or changed
    blocks={}
    for x in range(q):blocks.setdefault(root(x),[]).append(x)
    partition=sorted(blocks.values());r=len(partition);labels={x:j for j,b in enumerate(partition)for x in b}
    if len({len(b)for b in partition})!=1:raise ValueError('invariant block sizes')
    quotient=[tuple(labels[p[b[0]]]for b in partition)for p in(beta,gamma)];quotient.append(compose(*quotient))
    numerator=2-2*r+sum(r-len(cycles(p))for p in quotient)
    if numerator not in(0,2):raise ValueError('genus-one quotient must have intermediate genus0or1')
    return numerator//2,tuple(tuple(b)for b in partition)


def double_isogeny_quotients(beta,gamma):
    """All three nonzero H1(E,F2) characters, with zero ordinary returns."""
    f=IntrinsicFrame(beta,gamma);q=len(beta);result=[]
    for handles in((0,1),(1,0),(1,1)):
        z=[None]*q
        for j in f.primal:z[j]=0
        for j,value in zip(f.handles,handles):z[j]=value
        for face,edge in f.leaves:
            occurrences=[i for kind,i,_ in f.words[face]if kind=='z']
            if occurrences.count(edge)!=1:raise ValueError('double-isogeny leaf orientation')
            known=[z[i]for i in occurrences if i!=edge]
            if None in known:raise ValueError('unassigned double-isogeny edge')
            z[edge]=sum(known)%2
        if sum(z[i]for kind,i,_ in f.words[0]if kind=='z')%2:raise ValueError('double-isogeny root relation')
        b=tuple(2*beta[j]+x for j in range(q)for x in range(2))
        c=tuple(2*gamma[j]+(x+z[j]+z[gamma[j]])%2 for j in range(q)for x in range(2))
        for p,old in((b,beta),(c,gamma),(compose(b,c),compose(beta,gamma))):
            if sorted(map(len,cycles(p)))!=sorted(k for row in cycles(old)for k in(len(row),len(row))):raise ValueError('double-isogeny ordinary return')
        result.append((b,c,handles))
    return result


def triangular_quotient(a,b,c):
    """Index-ac lattice <(a,0),(b,c)> in the Euclidean333 triangle group."""
    return euclidean_quotient(3,a,b,c)


def euclidean_quotient(order,a,b,c):
    """Right cosets of a translation sublattice in triangle333/244/236."""
    if min(a,c)<1 or not 0<=b<a:raise ValueError('lattice Hermite normal form')
    steps={3:((1,0),(0,1),(-1,-1)),4:((1,0),(0,1),(-1,0),(0,-1)),6:((1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1))}[order]
    gamma_shift=2 if order==6 else 1
    beta=[];gamma=[]
    for y in range(c):
        for x in range(a):
            for r,(dx,dy)in enumerate(steps):
                beta.append(order*(a*y+x)+(r+1)%order)
                quotient,remainder=divmod(y+dy,c);nx=(x+dx-quotient*b)%a
                gamma.append(order*(a*remainder+nx)+(r+gamma_shift)%order)
    return tuple(beta),tuple(gamma)


def _collision_colors(t,special):
    q=len(t[0])//3;which=1-special if special<2 else 5-special
    a=t[special];ordinary=t[which]
    local=[];voltage=[];quotient=[];active=set()
    for j in range(q):
        if any(a[3*j+x]//3!=j for x in range(3)):raise ValueError('special quotient not identity')
        local.append(RANK[tuple(a[3*j+x]%3 for x in range(3))])
        if local[-1]!=0:
            if local[-1]not in TRANS:raise ValueError('moving inertia not a reflection')
            active.add(j)
        target=ordinary[3*j]//3
        if any(ordinary[3*j+x]//3!=target for x in range(3)):raise ValueError('ordinary not block preserving')
        quotient.append(target);voltage.append(RANK[tuple(ordinary[3*j+x]%3 for x in range(3))])
    colors=[]
    for row in cycles(quotient):
        selected=[j for j in row if j in active]
        if len(selected)>2:return None
        transports={row[0]:0}
        for j in row:
            nxt=quotient[j];value=MUL[voltage[j]][transports[j]]
            if nxt not in transports:transports[nxt]=value
            elif transports[nxt]!=value:raise ValueError('ordinary return not pure')
        if len(selected)==2:
            pulled=[conjugate(INV[transports[j]],local[j])for j in selected]
            colors.append((PERMS[pulled[0]][0]-PERMS[pulled[1]][0])%3)
    return colors


LOCAL_WORDS=((),)+tuple((m,)for m in range(6))+((4,0),(4,0,0))


def bad_forms(frame,mask,signs,word):
    special=0
    for m in word:
        k=m%3
        if special in(k,k+1):special=2*k+1-special
    vectors=((0,0,0,0),)+tuple(tuple(int(i==j)for i in range(4))for j in range(4))
    values=[_collision_colors(apply_word(frame.lift(mask,signs,v),word),special)for v in vectors]
    if values[0]is None:return((0,0,0,0),)
    if any(values[0]):raise ValueError('collision form has nonzero origin')
    if any(row is None or len(row)!=len(values[0])for row in values):raise ValueError('collision support depends on translations')
    return tuple(tuple(values[j+1][i]for j in range(4))for i in range(len(values[0])))


def initial_smooth(t):
    n=len(t[0]);indices=[n-len(cycles(p))for p in t]
    return n-len(cycles(compose(t[0],t[1])))in(indices[0]+indices[1],indices[2]+indices[3])


class _StateSpace:
    def __init__(self,masks):self.masks=tuple(masks)
    def __len__(self):return 160*len(self.masks)
    def __getitem__(self,i):
        if not 0<=i<len(self):raise IndexError(i)
        return self.masks[i//160],SIGN_STATES[(i%160)//40],PROJECTIVE[i%40]


class _Subset:
    def __init__(self,states,ids):self.states=states;self.ids=ids
    def __len__(self):return len(self.ids)
    def __getitem__(self,i):return self.states[self.ids[i]]


def greatest_good(frame,orbit,deadline=None):
    """Necessary subgroup closure; every rejected smooth start gets a word."""
    from array import array
    import time
    masks=sorted(orbit);states=_StateSpace(masks);n=len(states);lookup={m:i for i,m in enumerate(masks)}
    incoming=[array('i',[-1])*n for _ in PUSH_WORDS];bad=array('i',[-2])*n;parent_word=bytearray(n);local=bytearray(n);matrices={};queue=array('I')
    for mi,mask in enumerate(masks):
        if deadline is not None and time.time()>=deadline:raise TimeoutError('absolute research deadline in intrinsic matrix graph')
        for si,signs in enumerate(SIGN_STATES):
            conditions=[(word,bad_forms(frame,mask,signs,word))for word in LOCAL_WORDS]
            edges=[(word,frame.matrix(mask,signs,word))for word in PUSH_WORDS]
            matrices[mask,signs]=edges
            marked=0
            for wi,(word,forms)in enumerate(conditions):
                bits=0
                for form in forms:bits|=HYPERPLANES[form]
                fresh=bits&~marked;marked|=bits
                while fresh:
                    bit=fresh&-fresh;fresh-=bit;i=160*mi+40*si+bit.bit_length()-1
                    bad[i]=-1;local[i]=wi;queue.append(i)
            for ei,(word,(target,ns,matrix))in enumerate(edges):
                target_base=160*lookup[target]+40*SIGN_STATES.index(ns)
                for vi,nv in enumerate(projective_images(matrix)):
                    j=target_base+nv;i=160*mi+40*si+vi
                    if incoming[ei][j]!=-1:raise ValueError('point-push transition is not a permutation')
                    incoming[ei][j]=i
    if any(-1 in edges for edges in incoming):raise ValueError('point-push transition misses a vertex')
    head=0
    while head<len(queue):
        if deadline is not None and head%1024==0 and time.time()>=deadline:raise TimeoutError('absolute research deadline in intrinsic bad propagation')
        dest=queue[head];head+=1
        for ei in range(len(PUSH_WORDS)):
            source=incoming[ei][dest]
            if bad[source]!=-2:continue
            bad[source]=dest;parent_word[source]=ei;queue.append(source)
    def witness(i):
        word=[]
        if bad[i]==-2:raise ValueError('survivor has no negative witness')
        while bad[i]>=0:word.extend(PUSH_WORDS[parent_word[i]]);i=bad[i]
        word.extend(LOCAL_WORDS[local[i]]);return tuple(word)
    survivors=array('I',(i for i,b in enumerate(bad)if b==-2))
    return dict(states=states,rejected_count=n-len(survivors),survivors=_Subset(states,survivors),witness=witness,matrices=matrices,
                compact_graph_bytes=sum(len(v)*v.itemsize for v in incoming)+len(bad)*bad.itemsize+len(parent_word)+len(local)+len(queue)*queue.itemsize)
