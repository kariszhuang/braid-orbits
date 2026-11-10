"""Elliptic zero profile211: two moving simple values, one fixed3-cycle.

The fixed3-cycle is ONE branch value, not three separate transpositions.
All covers returned have full symmetric relative monodromy. The chart is
only a complete marked-fibre constructor, never an all-smooth certificate.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
from itertools import product

from genus_three_intrinsic import compose,inverse,cycles
from genus_three_intrinsic_native import MarkedCycleFrame,reduced,word_inverse


class ClusteredFrame(MarkedCycleFrame):
    def __init__(self,beta,gamma,fixed_vertex,m=5,deadline_epoch=None):
        super().__init__(beta,gamma,m=m,lazy=m>=6,deadline_epoch=deadline_epoch)
        self.vertices={'b':cycles(self.frame.beta),'g':cycles(self.frame.gamma),'d':self.frame.faces}
        self.fixed_vertex=(str(fixed_vertex[0]),int(fixed_vertex[1]))
        kind,vertex=self.fixed_vertex
        if kind not in self.vertices or not 0<=vertex<len(self.vertices[kind]):raise ValueError('unknown fixed elliptic value')
        self.threes=tuple(i for i,p in enumerate(self.perms)if sorted(map(len,cycles(p)))==[1]*(m-3)+[3])

    def symbolic(self,mask):
        moving=[j for j in range(self.q)if mask&(1<<j)]
        if len(moving)!=2:raise ValueError('two distinct moving simple values required')
        aa=[()]*self.q;bb=[()]*self.q;gg=[()]*self.q;dd=[()]*len(self.frame.faces);z=[None]*self.q
        for k,j in enumerate(moving):aa[j]=(k+3,)
        kind,vertex=self.fixed_vertex
        if kind=='d':dd[vertex]=(5,)
        else:(bb if kind=='b'else gg)[self.vertices[kind][vertex][-1]]=(5,)
        for j in self.frame.primal:z[j]=()
        for k,j in enumerate(self.frame.handles):z[j]=(k+1,)
        def face_word(face):
            out=()
            for j in self.frame.faces[face]:
                if z[j]is None or z[self.frame.gamma[j]]is None:raise ValueError('unassigned clustered edge')
                out=reduced(aa[self.frame.delta[j]]+bb[self.frame.gamma[j]]+word_inverse(z[self.frame.gamma[j]])+gg[j]+z[j]+out)
            return reduced(word_inverse(dd[face])+out)
        for face,edge in self.frame.leaves:
            z[edge]=(6,);word=face_word(face);cuts=[i for i,x in enumerate(word)if abs(x)==6]
            if len(cuts)!=1:raise ValueError('clustered dual variable must occur once')
            cut=cuts[0];solved=reduced(word_inverse(word[:cut])+word_inverse(word[cut+1:]));z[edge]=solved if word[cut]>0 else word_inverse(solved)
        root=face_word(0)
        if any(sum(abs(x)==j for x in root)!=1 for j in range(3,6)):raise ValueError('three actual punctures must each occur once')
        return moving,z,root

    def full_group(self,data):
        """Normal closure of the TWO actual local transpositions is Sm iff connected."""
        if len(data)!=5 or data[2]not in self.trans_set or data[3]not in self.trans_set or data[4]not in self.threes:
            raise ValueError('211 branch data required')
        parent=list(range(self.m))
        def root(x):
            while parent[x]!=x:parent[x]=parent[parent[x]];x=parent[x]
            return x
        def join(x,y):
            x,y=root(x),root(y)
            if x==y:return False
            parent[y]=x;return True
        for rank in data[2:4]:join(*(x for x,y in enumerate(self.perms[rank])if x!=y))
        changed=True
        while changed:
            changed=False;parts={}
            for x in range(self.m):parts.setdefault(root(x),[]).append(x)
            if len(parts)==1:return True
            for rank in(*data[:2],data[4]):
                p=self.perms[rank]
                for row in parts.values():
                    for x in row[1:]:changed=join(p[row[0]],p[x])or changed
        return len({root(x)for x in range(self.m)})==1

    def lift(self,moving,zwords,data):
        m=self.m;q=self.q;z=[self.value(w,data)for w in zwords];aa=[0]*q;bb=[0]*q;gg=[0]*q
        for k,j in enumerate(moving):aa[j]=data[k+2]
        kind,vertex=self.fixed_vertex
        if kind!='d':(bb if kind=='b'else gg)[self.vertices[kind][vertex][-1]]=data[4]
        a=tuple(m*j+self.perms[aa[j]][x]for j in range(q)for x in range(m))
        b=tuple(m*self.frame.beta[j]+self.perms[bb[j]][x]for j in range(q)for x in range(m))
        c=tuple(m*self.frame.gamma[j]+self.perms[self.mul[self.mul[self.inv[z[self.frame.gamma[j]]]][gg[j]]][z[j]]][x]for j in range(q)for x in range(m))
        d=inverse(compose(compose(a,b),c));t=(a,b,c,d);n=m*q
        if sum(n-len(cycles(p))for p in t)!=2*n+4 or n-len(cycles(a))!=2:raise ValueError('clustered source genus/index')
        for fine,name in((b,'b'),(c,'g'),(d,'d')):
            expected=[]
            for j,row in enumerate(self.vertices[name]):
                length=len(row);expected.extend(([3*length]+[length]*(m-3))if(name,j)==self.fixed_vertex else[length]*m)
            if sorted(map(len,cycles(fine)))!=sorted(expected):raise ValueError('clustered ordinary return type')
        return t

    def covers(self,mask):
        moving,zwords,root=self.symbolic(mask)
        if root.count(1)!=1 or root.count(-1)!=1:raise ValueError('clustered first handle not geometric')
        pivot=root.index(1);rotated=root[pivot:]+root[:pivot];end=rotated.index(-1);xword=rotated[1:end];yword=word_inverse(rotated[end+1:])
        number=trials=solutions=0;seen=set();handles=self.second_handles()
        for b in handles:
            stabilizer=[self.conjugate[h]for h in self.centralizer if self.mul[h][b]==self.mul[b][h]];local_seen=set()
            for locals in product(self.trans,self.threes):
                trials+=1
                if trials%1024==0:self.check_deadline()
                if locals in local_seen:continue
                orbit={tuple(conj[p]for p in locals)for conj in stabilizer};local_seen.update(orbit)
                if locals!=min(orbit):raise ValueError('noncanonical clustered local pair')
                data=(0,b,self.fixed,*locals);x=self.value(xword,data);y=self.value(yword,data)
                for a in self.conjugators(x,y):
                    candidate=(a,*data[1:]);solutions+=1
                    if not self.full_group(candidate):continue
                    key=self.rooted_key(candidate)
                    if key in seen:continue
                    seen.add(key);number+=1
                    if self.value(root,candidate):raise ValueError('clustered root relation')
                    yield number,self.lift(moving,zwords,candidate),dict(native_data=candidate)
        self.population=dict(second_handle_classes=len(handles),local_pair_trials=trials,solutions_before_A_gauge=solutions,canonical_full=number,
            root_word=root,z_words=zwords,centralizer_size=len(self.centralizer),algorithm='three-puncture211 exact conjugacy and marked-cycle gauge')
