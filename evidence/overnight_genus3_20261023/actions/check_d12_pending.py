#!/usr/bin/env python3
"""Independent literal/full-Sn-orbit audit of the retained S12 starts."""
import argparse
from collections import deque, Counter
import hashlib
import json
from pathlib import Path
import time

def compose(a,b):return tuple(a[j] for j in b)
def inverse(a):
    r=[0]*len(a)
    for i,j in enumerate(a):r[j]=i
    return tuple(r)
def cycle_count(a):
    seen=set();n=0
    for i in range(len(a)):
        if i in seen:continue
        n+=1
        while i not in seen:seen.add(i);i=a[i]
    return n
def index(a):return len(a)-cycle_count(a)
def smooth(t):return index(compose(t[0],t[1])) in (index(t[0])+index(t[1]),index(t[2])+index(t[3]))
def move(t,j):
    k,negative=divmod(j,2);a,b=t[k:k+2]
    pair=(b,compose(compose(inverse(b),a),b)) if negative else (compose(compose(a,b),inverse(a)),a)
    return t[:k]+pair+t[k+2:]
def replay(t,w):
    for j in w:t=move(t,j)
    return t
def key(t):
    n=len(t[0]);best=None
    for root in range(n):
        labels={root:0};order=[root]
        for x in order:
            for a in t:
                y=a[x]
                if y not in labels:labels[y]=len(order);order.append(y)
        assert len(order)==n,'full-product frame must be transitive'
        out=tuple(labels[a[x]] for a in t for x in order)
        if best is None or out<best:best=out
    return best
def search(t,words,cap):
    for w in words:
        if not smooth(replay(t,w)):return dict(status='negative',word=w,states_seen=0)
    q=deque([(t,())]);seen={key(t)}
    while q:
        u,w=q.popleft()
        if not smooth(u):
            words.append(w);return dict(status='negative',word=w,states_seen=len(seen))
        for j in range(6):
            v=move(u,j);z=key(v)
            if z not in seen:
                seen.add(z)
                if len(seen)>cap:return dict(status='unresolved',states_seen=len(seen))
                q.append((v,w+(j,)))
    return dict(status='closed',states=len(seen))
def main():
    p=argparse.ArgumentParser();p.add_argument('input',type=Path);p.add_argument('output',type=Path);p.add_argument('--cap',type=int,default=10000);a=p.parse_args()
    assert not a.output.exists();a.output.mkdir(parents=True)
    start=time.monotonic();words=[];counts=Counter();rows=[]
    with (a.output/'WITNESSES.jsonl').open('w') as out:
        for line in a.input.read_text().splitlines():
            row=json.loads(line);t=tuple(tuple(x)for x in row['tuple']);identity=tuple(range(12))
            assert compose(compose(compose(*t[:2]),t[2]),t[3])==identity
            assert tuple(map(index,t))==(2,9,8,9) and smooth(t)
            r=search(t,words,a.cap);counts[r['status']]+=1
            if r['status']=='negative':assert not smooth(replay(t,r['word']))
            z=dict(tuple=t,**r);rows.append(z);out.write(json.dumps(z)+'\n');out.flush()
    result=dict(counts=dict(counts),input_sha256=hashlib.sha256(a.input.read_bytes()).hexdigest(),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),elapsed_seconds=time.monotonic()-start,cap=a.cap,maximum_word_length=max((len(r.get('word',[]))for r in rows),default=0),literal_words=words)
    (a.output/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
