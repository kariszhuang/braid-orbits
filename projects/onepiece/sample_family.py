import random, subprocess, json, sys
import os
B=os.environ.get('ONEPIECE', os.path.join(os.path.dirname(os.path.abspath(__file__)),'target/release/onepiece'))
def mul(a,b): return [a[x] for x in b]
def inv(a):
    r=[0]*len(a)
    for i,x in enumerate(a): r[x]=i
    return r
def cyc(a):
    seen=set(); out=[]
    for s in range(len(a)):
        if s in seen: continue
        c=[];x=s
        while x not in seen: seen.add(x); c.append(x); x=a[x]
        out.append(c)
    return out
def fmt(a): return ''.join('('+' '.join(str(v+1) for v in c)+')' for c in cyc(a) if len(c)>1) or '()'
random.seed(int(sys.argv[2]))
d=int(sys.argv[1]); best={}
for trial in range(int(sys.argv[3])):
    def dt():
        s=random.sample(range(d),4); p=list(range(d)); p[s[0]],p[s[1]]=s[1],s[0]; p[s[2]],p[s[3]]=s[3],s[2]; return p
    x,y=dt(),dt()
    perm=random.sample(range(d),d); z=list(range(d))
    for i in range(d): z[perm[i]]=perm[(i+1)%d]
    w=inv(mul(x,mul(y,z)))
    if len(cyc(w))!=1: continue
    t=';'.join(fmt(p) for p in (x,y,z,w))
    o=json.loads(subprocess.run([B,'--orbit',t,'--cap','20000'],capture_output=True,text=True).stdout)
    if o.get('capped'): continue
    k=(o['states'],o['clusters'],o['bad_clusters'])
    if k not in best: best[k]=t
for k in sorted(best,key=lambda k:k[2])[:5]: print(d,'states,clusters,bad',k,best[k])
