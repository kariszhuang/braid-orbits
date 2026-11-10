"""Four simple elliptic branch values: two/three moving, the rest fixed.

This does not split a 2^2 inertia at one elliptic value into two punctures.
The two fixed values are DISTINCT black/white/face vertices. Any fixed plus
moving branch values on fixed E make the source family non-isotrivial.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
from genus_three_intrinsic import cycles,compose,inverse
from genus_three_intrinsic_native import MarkedCycleFrame,reduced,word_inverse


class MixedFrame(MarkedCycleFrame):
    def __init__(self,beta,gamma,fixed_vertices,m=3,deadline_epoch=None):
        super().__init__(beta,gamma,m=m,lazy=m>=6,deadline_epoch=deadline_epoch)
        self.vertices={'b':cycles(self.frame.beta),'g':cycles(self.frame.gamma),'d':self.frame.faces}
        self.fixed_vertices=tuple((str(kind),int(j))for kind,j in fixed_vertices)
        if len(self.fixed_vertices)not in(1,2) or len(set(self.fixed_vertices))!=len(self.fixed_vertices):raise ValueError('one or two distinct fixed elliptic values required')
        self.moving_count=4-len(self.fixed_vertices)
        if any(kind not in self.vertices or not 0<=j<len(self.vertices[kind])for kind,j in self.fixed_vertices):raise ValueError('unknown fixed vertex')

    def symbolic(self,mask):
        moving=[j for j in range(self.q)if mask&(1<<j)]
        if len(moving)!=self.moving_count:raise ValueError('moving/fixed simple values must total four')
        a=[()]*self.q;b=[()]*self.q;g=[()]*self.q;d=[()]*len(self.frame.faces);z=[None]*self.q
        for k,j in enumerate(moving):a[j]=(k+3,)
        for k,(kind,vertex)in enumerate(self.fixed_vertices):
            if kind=='d':d[vertex]=(k+self.moving_count+3,)
            else:(b if kind=='b'else g)[self.vertices[kind][vertex][-1]]=(k+self.moving_count+3,)
        for j in self.frame.primal:z[j]=()
        for k,j in enumerate(self.frame.handles):z[j]=(k+1,)
        def face_word(face):
            out=()
            for j in self.frame.faces[face]:
                if z[j]is None or z[self.frame.gamma[j]]is None:raise ValueError('unassigned mixed edge')
                out=reduced(a[self.frame.delta[j]]+b[self.frame.gamma[j]]+word_inverse(z[self.frame.gamma[j]])+g[j]+z[j]+out)
            return reduced(word_inverse(d[face])+out)
        for face,edge in self.frame.leaves:
            z[edge]=(7,);word=face_word(face);cuts=[i for i,x in enumerate(word)if abs(x)==7]
            if len(cuts)!=1:raise ValueError('mixed dual edge must occur once')
            cut=cuts[0];solved=reduced(word_inverse(word[:cut])+word_inverse(word[cut+1:]));z[edge]=solved if word[cut]>0 else word_inverse(solved)
        root=face_word(0)
        if any(sum(abs(x)==j for x in root)!=1 for j in range(3,7)):raise ValueError('mixed four punctures must each occur once')
        return moving,z,root

    def lift(self,moving,zwords,data):
        q=self.q;m=self.m;z=[self.value(w,data)for w in zwords];a=[0]*q;b=[0]*q;g=[0]*q
        for k,j in enumerate(moving):a[j]=data[k+2]
        for k,(kind,vertex)in enumerate(self.fixed_vertices):
            if kind!='d':(b if kind=='b'else g)[self.vertices[kind][vertex][-1]]=data[k+self.moving_count+2]
        aa=tuple(m*j+self.perms[a[j]][x]for j in range(q)for x in range(m))
        bb=tuple(m*self.frame.beta[j]+self.perms[b[j]][x]for j in range(q)for x in range(m))
        cc=tuple(m*self.frame.gamma[j]+self.perms[self.mul[self.mul[self.inv[z[self.frame.gamma[j]]]][g[j]]][z[j]]][x]for j in range(q)for x in range(m))
        dd=inverse(compose(compose(aa,bb),cc));t=(aa,bb,cc,dd);n=m*q
        if sum(n-len(cycles(p))for p in t)!=2*n+4 or n-len(cycles(aa))!=self.moving_count:raise ValueError('mixed source genus/moving index')
        for fine,kind in((bb,'b'),(cc,'g'),(dd,'d')):
            expected=[]
            for vertex,row in enumerate(self.vertices[kind]):
                length=len(row)
                expected.extend(([2*length]+[length]*(m-2))if(kind,vertex)in self.fixed_vertices else[length]*m)
            if sorted(map(len,cycles(fine)))!=sorted(expected):raise ValueError('mixed ordinary local return mismatch')
        seen={0};queue=[0]
        for x in queue:
            for p in t[:3]:
                if p[x]not in seen:seen.add(p[x]);queue.append(p[x])
        if len(seen)!=n:raise ValueError('full relative monodromy did not induce connected fine cover')
        return t


def fixed_vertex_placements(beta,gamma,count=2):
    from itertools import combinations
    if count not in(1,2):raise ValueError('one or two fixed vertices required')
    vertices=[(kind,j)for kind,p in(('b',beta),('g',gamma),('d',compose(beta,gamma)))for j in range(len(cycles(p)))]
    yield from combinations(vertices,count)


def fixed_vertex_pairs(beta,gamma):
    yield from fixed_vertex_placements(beta,gamma,2)
