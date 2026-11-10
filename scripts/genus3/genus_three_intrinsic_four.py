"""Full-S4 elliptic covers:160 cubic resolvents times15 V4 extensions.

No cubic smoothness test is used: a bad cubic resolvent can have a smooth
quartic lift, since disjoint quartic transpositions share a cubic color.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
from functools import lru_cache

from genus_three_intrinsic import PERMS,MUL,cycles,compose,inverse
from genus_three_intrinsic_linear import LinearFrame,PROJECTIVE,SIGN_STATES,affine

LINEAR=tuple((0,p[0]+1,p[1]+1,p[2]+1)for p in PERMS)
FOUR=tuple(tuple(LINEAR[p][x]^u for x in range(4))for p in range(6)for u in range(4))
MUL4=tuple(tuple(4*MUL[a//4][b//4]+((a%4)^LINEAR[a//4][b%4])for b in range(24))for a in range(24))
INV4=tuple(next(b for b in range(24)if MUL4[a][b]==0)for a in range(24))
RANK3={p:i for i,p in enumerate(PERMS)}


def relative_order(generators):
    seen={0};queue=[0]
    for h in queue:
        for g in generators:
            k=MUL4[g][h]
            if k not in seen:seen.add(k);queue.append(k)
    return len(seen)


class QuarticFrame:
    def __init__(self,beta,gamma):self.cubic=LinearFrame(beta,gamma);self.frame=self.cubic.frame;self.q=len(beta)

    @staticmethod
    def value(word,z,local):
        result=0
        for kind,j,negative in word:
            v=z[j]if kind=='z'else local[j]
            if v is None:raise ValueError('unassigned quartic potential')
            if negative:v=INV4[v]
            result=MUL4[v][result]
        return result

    def algebra(self,local,handles):
        z=[None]*self.q
        for j in self.frame.primal:z[j]=0
        for j,v in zip(self.frame.handles,handles):z[j]=v
        for face,edge in self.frame.leaves:
            word=self.frame.words[face];positions=[k for k,(kind,j,_)in enumerate(word)if kind=='z'and j==edge]
            if len(positions)!=1:raise ValueError('quartic dual leaf')
            cut=positions[0];before=self.value(word[:cut],z,local);after=self.value(word[cut+1:],z,local)
            v=MUL4[INV4[after]][INV4[before]]
            if word[cut][2]:v=INV4[v]
            z[edge]=v
            if self.value(word,z,local)!=0:raise ValueError('quartic leaf equation')
        return z,self.value(self.frame.words[0],z,local)

    def cubic_parameters(self,mask,signs,vector):
        t=self.cubic.lift(mask,signs,vector)
        local=tuple(RANK3[tuple(t[0][3*j+x]%3 for x in range(3))]for j in range(self.q))
        return local,(affine(signs[0],vector[0]),affine(signs[1],vector[1]))

    def extension_space(self,mask,signs,vector,check_all=False):
        local3,handle3=self.cubic_parameters(mask,signs,vector);active=[j for j in range(self.q)if mask&(1<<j)]
        axes=[next(x for x in(1,2,3)if LINEAR[local3[j]][x]==x)for j in active]
        def values(code):
            local=[4*p for p in local3]
            for k,(j,axis)in enumerate(zip(active,axes)):
                if code&(1<<(4+k)):local[j]^=axis
            return local,(4*handle3[0]+(code&3),4*handle3[1]+((code>>2)&3))
        columns=[]
        for j in range(8):
            _,root=self.algebra(*values(1<<j))
            if root//4:raise ValueError('quartic root has nonidentity cubic projection')
            columns.append(root%4)
        def relation(code):
            result=0
            for j in range(8):
                if code&(1<<j):result^=columns[j]
            return result
        kernel=[code for code in range(256)if relation(code)==0]
        if len(kernel)!=64:raise ValueError('binary root rank is not two')
        gauge=[]
        for c in(1,2):
            code=(c^LINEAR[handle3[0]][c])|((c^LINEAR[handle3[1]][c])<<2)
            for k,(j,axis)in enumerate(zip(active,axes)):
                delta=c^LINEAR[local3[j]][c]
                if delta not in(0,axis):raise ValueError('local translation gauge not reflection compatible')
                if delta:code|=1<<(4+k)
            gauge.append(code)
        translations=(0,gauge[0],gauge[1],gauge[0]^gauge[1])
        if len(set(translations))!=4 or any(relation(g)for g in translations):raise ValueError('binary gauge rank/kernel')
        reps=sorted({min(code^g for g in translations)for code in kernel})
        if len(reps)!=16 or reps[0]!=0:raise ValueError('quartic extension quotient is not F2^4')
        if check_all:
            for code in range(256):
                _,root=self.algebra(*values(code))
                if root//4 or root%4!=relation(code):raise ValueError('binary root map not linear')
        return reps[1:],values,dict(columns=columns,gauge=gauge,kernel_size=len(kernel),classes=16)

    def lift(self,local,handles):
        z,root=self.algebra(local,handles)
        if root:raise ValueError('quartic root relation')
        if relative_order((*handles,*(p for p in local if p)))!=24:raise ValueError('quartic extension is not full S4')
        q=self.q;n=4*q
        a=tuple(4*j+FOUR[local[j]][x]for j in range(q)for x in range(4))
        b=tuple(4*self.frame.beta[j]+x for j in range(q)for x in range(4))
        c=tuple(4*self.frame.gamma[j]+FOUR[MUL4[INV4[z[self.frame.gamma[j]]]][z[j]]][x]for j in range(q)for x in range(4))
        d=inverse(compose(compose(a,b),c));t=(a,b,c,d)
        if sum(n-len(cycles(p))for p in t)!=2*n+4:raise ValueError('quartic source genus not three')
        for p,bar in((b,self.frame.beta),(c,self.frame.gamma),(d,inverse(self.frame.delta))):
            if sorted(map(len,cycles(p)))!=sorted(k for row in cycles(bar)for k in[len(row)]*4):raise ValueError('quartic ordinary return not pure')
        return t

    def covers(self,mask,check_all=False):
        number=0
        for signs in SIGN_STATES:
            for vector in PROJECTIVE:
                reps,values,meta=self.extension_space(mask,signs,vector,check_all=check_all)
                for code in reps:
                    number+=1;local,handles=values(code);yield number,self.lift(local,handles),dict(signs=signs,cubic_vector=vector,binary_code=code)
        if number!=2400:raise ValueError('full quartic fiber not2400')
