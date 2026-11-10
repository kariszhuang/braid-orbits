"""Native symmetric-group tree--cotree coordinates, with no affine model.

The first local transposition is fixed. Its centralizer is the remaining
simultaneous-conjugation gauge. This is a constructor, not a braid certificate.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import itertools
import time
from functools import lru_cache

from genus_three_intrinsic import IntrinsicFrame,compose,inverse,cycles


def reduced(word):
    out=[]
    for x in word:
        if out and out[-1]==-x:out.pop()
        else:out.append(x)
    return tuple(out)


def word_inverse(word):return tuple(-x for x in reversed(word))


class ProductRow:
    def __init__(self,product,a):self.product=product;self.a=a
    def __getitem__(self,b):return self.product(self.a,b)


class ConjugateRow:
    def __init__(self,frame,h):self.frame=frame;self.h=h
    def __getitem__(self,p):return self.frame.mul[self.frame.mul[self.h][p]][self.frame.inv[self.h]]


def partitions(total,minimum=1):
    if not total:yield();return
    for k in range(minimum,total+1):
        for rest in partitions(total-k,k):yield(k,)+rest


class NativeFrame:
    def __init__(self,beta,gamma,m=5,lazy=False,deadline_epoch=None):
        if not 3<=m<=8:raise ValueError('bounded native fibre degree3..8')
        self.frame=IntrinsicFrame(beta,gamma);self.q=len(beta);self.m=m;self.deadline_epoch=deadline_epoch
        self.perms=tuple(itertools.permutations(range(m)));self.rank={p:i for i,p in enumerate(self.perms)};self.size=len(self.perms)
        if lazy or m>=8:
            @lru_cache(maxsize=500000)
            def product(a,b):return self.rank[compose(self.perms[a],self.perms[b])]
            self.mul=tuple(ProductRow(product,a)for a in range(self.size))
        else:self.mul=tuple(tuple(self.rank[compose(a,b)]for b in self.perms)for a in self.perms)
        self.inv=tuple(self.rank[inverse(p)]for p in self.perms)
        self.trans=tuple(i for i,p in enumerate(self.perms)if sum(x!=p[x]for x in range(m))==2)
        self.trans_set=frozenset(self.trans)
        self.fixed=self.trans[0]
        self.centralizer=tuple(h for h in range(self.size)if self.mul[h][self.fixed]==self.mul[self.fixed][h])
        self.conjugate={h:ConjugateRow(self,h)if lazy or m>=8 else tuple(self.mul[self.mul[h][p]][self.inv[h]]for p in range(self.size))for h in self.centralizer}
        self.population=None

    def check_deadline(self):
        if self.deadline_epoch is not None and time.time()>=self.deadline_epoch:
            raise TimeoutError('absolute research drain in native fibre')

    @staticmethod
    def symbolic_value(word,z,local):
        out=()
        for kind,j,negative in word:
            value=z[j]if kind=='z'else local[j]
            if value is None:raise ValueError('unassigned symbolic edge')
            out=reduced((word_inverse(value)if negative else value)+out)
        return out

    def symbolic(self,mask):
        active=[j for j in range(self.q)if mask&(1<<j)]
        if len(active)!=4:raise ValueError('four simple branch values required')
        z=[None]*self.q;local=[()]*self.q
        for j in self.frame.primal:z[j]=()
        for k,j in enumerate(self.frame.handles):z[j]=(k+1,)
        for k,j in enumerate(active):local[j]=(k+3,)
        for face,edge in self.frame.leaves:
            z[edge]=(7,);word=self.symbolic_value(self.frame.words[face],z,local);cuts=[i for i,x in enumerate(word)if abs(x)==7]
            if len(cuts)!=1:raise ValueError('dual variable does not occur exactly once')
            cut=cuts[0];solved=reduced(word_inverse(word[:cut])+word_inverse(word[cut+1:]))
            z[edge]=solved if word[cut]>0 else word_inverse(solved)
        root=self.symbolic_value(self.frame.words[0],z,local)
        if sum(abs(x)==6 for x in root)!=1:raise ValueError('fourth puncture must occur once')
        return active,z,root

    def value(self,word,assignment):
        result=0
        for token in word:result=self.mul[result][assignment[token-1]if token>0 else self.inv[assignment[-token-1]]]
        return result

    def full_group_reference(self,data):
        seen={0};queue=[0]
        for x in queue:
            for g in data:
                y=self.mul[x][g]
                if y not in seen:seen.add(y);queue.append(y)
        return len(seen)==self.size

    def full_group(self,data):
        """The transposition normal-closure graph is connected iff G=S_m."""
        if len(data)!=6 or any(p not in self.trans_set for p in data[2:]):raise ValueError('two handles/four transpositions required')
        parent=list(range(self.m))
        def root(x):
            while parent[x]!=x:parent[x]=parent[parent[x]];x=parent[x]
            return x
        def join(x,y):
            x,y=root(x),root(y)
            if x==y:return False
            parent[y]=x;return True
        for rank in data[2:]:
            moved=[x for x,y in enumerate(self.perms[rank])if x!=y];join(*moved)
        changed=True
        while changed:
            changed=False;parts={}
            for x in range(self.m):parts.setdefault(root(x),[]).append(x)
            if len(parts)==1:return True
            for rank in data[:2]:
                p=self.perms[rank]
                for row in parts.values():
                    for x in row[1:]:changed=join(p[row[0]],p[x])or changed
        return len({root(x)for x in range(self.m)})==1

    def rooted_key(self,data):
        """Exact simultaneous gauge key using two roots on the first tau."""
        order=(data[2],*data[:2],*data[3:]);best=None
        for start in range(self.m):
            if self.perms[data[2]][start]==start:continue
            labels={start:0};queue=[start]
            for x in queue:
                for rank in order:
                    y=self.perms[rank][x]
                    if y not in labels:labels[y]=len(queue);queue.append(y)
            if len(queue)!=self.m:raise ValueError('rooted key requires transitivity')
            key=bytes(labels[self.perms[rank][x]]for rank in data for x in queue)
            if best is None or key<best:best=key
        return best

    def lift(self,active,zwords,data):
        m=self.m;q=self.q;z=[self.value(w,data)for w in zwords];local=[0]*q
        for k,j in enumerate(active):local[j]=data[k+2]
        a=tuple(m*j+self.perms[local[j]][x]for j in range(q)for x in range(m))
        b=tuple(m*self.frame.beta[j]+x for j in range(q)for x in range(m))
        c=tuple(m*self.frame.gamma[j]+self.perms[self.mul[self.inv[z[self.frame.gamma[j]]]][z[j]]][x]for j in range(q)for x in range(m))
        d=inverse(compose(compose(a,b),c));t=(a,b,c,d);n=m*q
        if sum(n-len(cycles(p))for p in t)!=2*n+4:raise ValueError('native source genus not three')
        for p,bar in((b,self.frame.beta),(c,self.frame.gamma),(d,inverse(self.frame.delta))):
            if sorted(map(len,cycles(p)))!=sorted(len(row)for row in cycles(bar)for _ in range(m)):raise ValueError('native ordinary return is not pure')
        return t

    def covers(self,mask):
        active,zwords,root=self.symbolic(mask);cut=next(i for i,x in enumerate(root)if abs(x)==6);seen=set();number=raw=valid=full=0
        for h1,h2 in itertools.product(range(self.size),repeat=2):
            for t2,t3 in itertools.product(self.trans,repeat=2):
                raw+=1
                if raw%4096==0:self.check_deadline()
                data=(h1,h2,self.fixed,t2,t3,0)
                before=self.value(root[:cut],data);after=self.value(root[cut+1:],data)
                last=self.mul[self.inv[before]][self.inv[after]]
                if root[cut]<0:last=self.inv[last]
                if last not in self.trans:continue
                valid+=1;data=(*data[:5],last)
                key=min(tuple(conj[p]for p in data)for conj in self.conjugate.values())
                if key in seen:continue
                seen.add(key)
                if not self.full_group(key):continue
                full+=1;number+=1
                if self.value(root,key):raise ValueError('native root relation failed')
                yield number,self.lift(active,zwords,key),dict(native_data=key)
        self.population=dict(raw_fixed_first=raw,valid_fixed_first=valid,canonical_all=len(seen),canonical_full=full,root_word=root,z_words=zwords,centralizer_size=len(self.centralizer))


class ConjugacyFrame(NativeFrame):
    """Solve A X A^-1=Y by exact cycle matchings, instead of looping A.

    Gauge order: fix first transposition; choose B modulo its centralizer;
    choose the other three transpositions modulo C(first,B); finally choose
    the remaining A modulo the stabilizer of all these already chosen data.
    Each step is an exact finite orbit quotient, not a random proposal.
    """
    @lru_cache(maxsize=1024)
    def conjugators(self,x,y):
        cx=cycles(self.perms[x]);cy=cycles(self.perms[y])
        if sorted(map(len,cx))!=sorted(map(len,cy)):return()
        lengths=sorted({len(row)for row in cx});partial=[tuple(range(self.m))]
        for length in lengths:
            left=[row for row in cx if len(row)==length];right=[row for row in cy if len(row)==length];expanded=[]
            for ordering in itertools.permutations(range(len(right))):
                for offsets in itertools.product(range(length),repeat=len(left)):
                    for old in partial:
                        h=list(old)
                        for j,row in enumerate(left):
                            target=right[ordering[j]]
                            for k,point in enumerate(row):h[point]=target[(k+offsets[j])%length]
                        expanded.append(tuple(h))
            partial=expanded
        result=tuple(self.rank[p]for p in partial)
        if len(result)!=len(set(result)):raise ValueError('duplicate cycle conjugator')
        if any(self.mul[self.mul[a][x]][self.inv[a]]!=y for a in result):raise ValueError('cycle conjugator orientation')
        return result

    def covers(self,mask):
        active,zwords,root=self.symbolic(mask)
        if root.count(1)!=1 or root.count(-1)!=1:raise ValueError('first handle not a geometric handle')
        pivot=root.index(1);rotated=root[pivot:]+root[:pivot];end=rotated.index(-1);xword=rotated[1:end];yword=word_inverse(rotated[end+1:])
        visited=set();number=trials=solutions=0;b_classes=0
        for b in range(self.size):
            if b in visited:continue
            orbit={c[b]for c in self.conjugate.values()};visited.update(orbit)
            if b!=min(orbit):raise ValueError('noncanonical second handle traversal')
            b_classes+=1;stabilizer=[c for c in self.conjugate.values()if c[b]==b];triple_seen=set()
            for triple in itertools.product(self.trans,repeat=3):
                trials+=1
                if trials%1024==0:self.check_deadline()
                if triple in triple_seen:continue
                orbit={tuple(c[p]for p in triple)for c in stabilizer};triple_seen.update(orbit)
                if triple!=min(orbit):raise ValueError('noncanonical local triple traversal')
                remaining=[c for c in stabilizer if all(c[p]==p for p in triple)]
                data=(0,b,self.fixed,*triple);x=self.value(xword,data);y=self.value(yword,data)
                for a in self.conjugators(x,y):
                    if a!=min(c[a]for c in remaining):continue
                    candidate=(a,*data[1:]);solutions+=1
                    if not self.full_group(candidate):continue
                    number+=1
                    if self.value(root,candidate):raise ValueError('conjugacy root relation failed')
                    yield number,self.lift(active,zwords,candidate),dict(native_data=candidate)
        self.population=dict(second_handle_classes=b_classes,local_triple_trials=trials,canonical_all=solutions,canonical_full=number,root_word=root,z_words=zwords,centralizer_size=len(self.centralizer),algorithm='exact conjugacy equation with stabilizer chain')


class MarkedCycleFrame(ConjugacyFrame):
    """Direct C(tau) second-handle classes; no factorial B-orbit traversal."""
    def second_handles(self):
        marked=[x for x,y in enumerate(self.perms[self.fixed])if x!=y];first,second=marked;unmarked=[x for x in range(self.m)if x not in marked];result=[]
        def permutation(rows):
            p=list(range(self.m))
            for row in rows:
                for x,y in zip(row,row[1:]+row[:1]):p[x]=y
            return self.rank[tuple(p)]
        def leftover(rows,used,shape):
            spare=unmarked[used:];at=0;out=list(rows)
            for length in shape:out.append(tuple(spare[at:at+length]));at+=length
            return permutation(out)
        # Marked endpoints in one oriented B-cycle. Exchanging them identifies
        # separation k with length-k, and no further oriented-cycle invariant.
        for length in range(2,self.m+1):
            for separation in range(1,length//2+1):
                row=(first,*unmarked[:separation-1],second,*unmarked[separation-1:length-2])
                for shape in partitions(self.m-length):result.append(leftover([row],length-2,shape))
        # Marked endpoints in different cycles; exchanging endpoints lets us
        # order their two lengths. The other B-cycles are unlabeled.
        for left in range(1,self.m):
            for right in range(left,self.m-left+1):
                row1=(first,*unmarked[:left-1]);row2=(second,*unmarked[left-1:left+right-2])
                for shape in partitions(self.m-left-right):result.append(leftover([row1,row2],left+right-2,shape))
        if len(result)!=len(set(result)):raise ValueError('duplicate marked-cycle handle')
        return result

    def covers(self,mask):
        active,zwords,root=self.symbolic(mask)
        if root.count(1)!=1 or root.count(-1)!=1:raise ValueError('first handle not geometric')
        pivot=root.index(1);rotated=root[pivot:]+root[:pivot];end=rotated.index(-1);xword=rotated[1:end];yword=word_inverse(rotated[end+1:]);number=trials=solutions=0;handles=self.second_handles();full_seen=set()
        for b in handles:
            stabilizer=[self.conjugate[h]for h in self.centralizer if self.mul[h][b]==self.mul[b][h]];triple_seen=set()
            for triple in itertools.product(self.trans,repeat=3):
                trials+=1
                if trials%1024==0:self.check_deadline()
                if triple in triple_seen:continue
                orbit={tuple(c[p]for p in triple)for c in stabilizer};triple_seen.update(orbit)
                if triple!=min(orbit):raise ValueError('noncanonical local triple')
                data=(0,b,self.fixed,*triple);x=self.value(xword,data);y=self.value(yword,data)
                for a in self.conjugators(x,y):
                    candidate=(a,*data[1:]);solutions+=1
                    if not self.full_group(candidate):continue
                    key=self.rooted_key(candidate)
                    if key in full_seen:continue
                    full_seen.add(key)
                    number+=1
                    if self.value(root,candidate):raise ValueError('marked-cycle root relation')
                    yield number,self.lift(active,zwords,candidate),dict(native_data=candidate)
        self.population=dict(second_handle_classes=len(handles),local_triple_trials=trials,solutions_before_A_gauge=solutions,canonical_full=number,root_word=root,z_words=zwords,centralizer_size=len(self.centralizer),algorithm='marked cycle layout + conjugacy equation + exact two-root gauge')
