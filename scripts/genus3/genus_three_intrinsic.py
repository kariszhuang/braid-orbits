"""Independent tree--cotree elliptic S3-cover coordinates.

This is a mathematical constructor/reference, not a braid-closure certificate.
"""
import itertools


PERMS=tuple(itertools.permutations(range(3)))
RANK={p:i for i,p in enumerate(PERMS)}
IDENTITY=RANK[(0,1,2)]
MUL=tuple(tuple(RANK[tuple(a[x]for x in b)]for b in PERMS)for a in PERMS)
INV=tuple(next(j for j in range(6)if MUL[i][j]==IDENTITY)for i in range(6))
TRANS=tuple(i for i,p in enumerate(PERMS)if sum(p[x]!=x for x in range(3))==2)


def compose(a,b):return tuple(a[x]for x in b)


def inverse(p):return tuple(p.index(i)for i in range(len(p)))


def cycles(p):
    seen=set();result=[]
    for start in range(len(p)):
        if start in seen:continue
        row=[];at=start
        while at not in seen:seen.add(at);row.append(at);at=p[at]
        result.append(tuple(row))
    return tuple(result)


class DSU:
    def __init__(self,n):self.parent=list(range(n))
    def root(self,x):
        while self.parent[x]!=x:self.parent[x]=self.parent[self.parent[x]];x=self.parent[x]
        return x
    def join(self,a,b):
        a,b=self.root(a),self.root(b)
        if a==b:return False
        self.parent[b]=a;return True


class IntrinsicFrame:
    def __init__(self,beta,gamma):
        self.beta=tuple(beta);self.gamma=tuple(gamma);self.q=q=len(beta)
        self.delta=compose(beta,gamma)
        black=cycles(beta);white=cycles(gamma);self.faces=faces=cycles(self.delta)
        if len(black)+len(white)+len(faces)!=q:raise ValueError('quotient must have genus one')
        bl={x:j for j,row in enumerate(black)for x in row}
        wh={x:j for j,row in enumerate(white)for x in row}
        fa={x:j for j,row in enumerate(faces)for x in row}
        primal=DSU(len(black)+len(white));self.primal=[]
        for e in range(q):
            if primal.join(bl[e],len(black)+wh[e]):self.primal.append(e)
        if len(self.primal)!=len(black)+len(white)-1:raise ValueError('disconnected quotient dessin')
        gamma_inverse=inverse(gamma);dual=DSU(len(faces));self.dual=[];adj=[[]for _ in faces]
        for e in range(q):
            if e in self.primal:continue
            a,b=fa[e],fa[gamma_inverse[e]]
            if dual.join(a,b):self.dual.append(e);adj[a].append((b,e));adj[b].append((a,e))
        if len(self.dual)!=len(faces)-1:raise ValueError('dual complement is not connected')
        self.handles=tuple(e for e in range(q)if e not in self.primal and e not in self.dual)
        if len(self.handles)!=2:raise ValueError('elliptic handle count')
        todo=[0];parent={0:(0,-1)}
        for at in todo:
            for nxt,edge in adj[at]:
                if nxt not in parent:parent[nxt]=(at,edge);todo.append(nxt)
        self.leaves=tuple((face,parent[face][1])for face in reversed(todo[1:]))
        self.words=[]
        for face in faces:
            word=[]
            for j in face:word.extend((('z',j,False),('z',gamma[j],True),('a',self.delta[j],False)))
            self.words.append(tuple(word))

    @staticmethod
    def value(word,z,local):
        result=IDENTITY
        for kind,j,negative in word:
            value=z[j]if kind=='z'else local[j]
            if value is None:raise ValueError('unassigned word variable')
            if negative:value=INV[value]
            result=MUL[value][result]
        return result

    def solve(self,local,handles):
        if len(local)!=self.q or len(handles)!=2:raise ValueError('intrinsic input shape')
        z=[None]*self.q
        for j in self.primal:z[j]=IDENTITY
        for j,v in zip(self.handles,handles):z[j]=v
        for face,edge in self.leaves:
            word=self.words[face]
            where=[j for j,(kind,i,_)in enumerate(word)if kind=='z'and i==edge]
            if len(where)!=1:raise ValueError('dual leaf variable must occur once')
            cut=where[0]
            before=self.value(word[:cut],z,local);after=self.value(word[cut+1:],z,local)
            unknown=MUL[INV[after]][INV[before]]
            if word[cut][2]:unknown=INV[unknown]
            z[edge]=unknown
            if self.value(word,z,local)!=IDENTITY:raise AssertionError('leaf equation orientation')
        if self.value(self.words[0],z,local)!=IDENTITY:return None
        voltage=[MUL[INV[z[self.gamma[j]]]][z[j]]for j in range(self.q)]
        n=3*self.q
        a=tuple(3*j+PERMS[local[j]][x]for j in range(self.q)for x in range(3))
        b=tuple(3*self.beta[j]+x for j in range(self.q)for x in range(3))
        c=tuple(3*self.gamma[j]+PERMS[voltage[j]][x]for j in range(self.q)for x in range(3))
        d=inverse(compose(compose(a,b),c))
        assert compose(compose(a,b),compose(c,d))==tuple(range(n))
        for p,quotient in ((b,self.beta),(c,self.gamma),(d,inverse(self.delta))):
            assert sorted(map(len,cycles(p)))==sorted(length for row in cycles(quotient)for length in [len(row)]*3)
        return a,b,c,d


def uniform3_quotient(q):
    if q%3:raise ValueError('degree divisible by three required')
    return tuple((x//3)*3+(x+1)%3 for x in range(q))
