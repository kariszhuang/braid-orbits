#!/usr/bin/env python3
"""Exact necessary marked S12 slice: indices(2,9,8,9), special a a3-cycle.

Coverage uses C209's saturated-single-cycle theorem and all initial twists.
No block system is assumed. Caps are unresolved and cannot yield exclusion.
"""
import argparse
from collections import Counter, deque
import hashlib
import itertools
import json
from pathlib import Path
import time

N = 12
ID = tuple(range(N))


def mul(a,b): return tuple(a[b[i]] for i in range(N))
def inv(a):
    z = [0]*N
    for i,j in enumerate(a): z[j] = i
    return tuple(z)
def cycles(a):
    seen, out = 0, []
    for i in range(N):
        if seen>>i&1: continue
        c = []
        while not seen>>i&1:
            seen |= 1<<i; c.append(i); i=a[i]
        out.append(c)
    return out
def ix(a): return N-len(cycles(a))
def rainbow(a): return all(len({x%4 for x in c})==len(c) for c in cycles(a))
def smooth(t):
    prod_ix = ix(mul(t[0],t[1]))
    return prod_ix in (ix(t[0])+ix(t[1]),ix(t[2])+ix(t[3]))
def braid(t,code):
    j,neg=divmod(code,2);a,b=t[j:j+2]
    pair=(b,mul(mul(inv(b),a),b)) if neg else (mul(mul(a,b),inv(a)),a)
    return t[:j]+pair+t[j+2:]
def replay(t,w):
    for code in w: t=braid(t,code)
    return t


def proof(start,known,cap):
    for w in known:
        if not smooth(replay(start,w)): return 'negative',w,0
    seen={start}; q=deque([(start,())])
    while q:
        t,w=q.popleft()
        if not smooth(t):
            if w not in known: known.append(w)
            return 'negative',w,len(seen)
        for code in range(6):
            u=braid(t,code)
            if u not in seen:
                seen.add(u)
                if len(seen)>cap: return 'unresolved',(),len(seen)
                q.append((u,w+(code,)))
    return 'labelled_closed',(),len(seen)


def enumerate_e():
    assignments=list(itertools.permutations(range(3)))
    orders=list(itertools.permutations((1,2,3)))
    for maps in itertools.product(assignments,repeat=3):
        blocks=[(4*i,4*maps[0][i]+1,4*maps[1][i]+2,4*maps[2][i]+3) for i in range(3)]
        for os in itertools.product(orders,repeat=3):
            e=list(ID)
            for block,order in zip(blocks,os):
                cyc=(block[0],)+tuple(block[i] for i in order)
                for x,y in zip(cyc,cyc[1:]+cyc[:1]):e[x]=y
            yield tuple(e)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--cap',type=int,default=5000)
    parser.add_argument('--seconds',type=float,default=120)
    args=parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=False)
    p=tuple(list(range(1,N))+[0]);pi=inv(p)
    a=list(ID);a[0],a[4],a[8]=4,8,0;a=tuple(a)
    c=mul(inv(a),p)
    assert [len(x) for x in cycles(c)]==[4,4,4]
    started=time.monotonic();counts=Counter();known=[];digest=hashlib.sha256()
    exhausted=True;unresolved=[]
    with (args.out/'WITNESSES.jsonl').open('w') as output:
        for e in enumerate_e():
            if time.monotonic()-started>args.seconds:
                exhausted=False;break
            counts['e_constructed']+=1;digest.update(bytes(e))
            assert rainbow(e)
            b=mul(pi,inv(e))
            if ix(b)!=8:counts['wrong_b_index']+=1;continue
            if not rainbow(b):counts['b_twist_failure']+=1;continue
            t=(a,c,b,e)
            assert mul(mul(mul(a,c),b),e)==ID
            assert tuple(map(ix,t))==(2,9,8,9) and smooth(t)
            counts['qualified']+=1
            status,word,states=proof(t,known,args.cap)
            counts[status]+=1
            if status=='negative':assert not smooth(replay(t,word))
            else:unresolved.append(dict(e=e,tuple=t,status=status,states=states))
            output.write(json.dumps(dict(e=e,status=status,word=word,states_seen=states))+'\n')
    if exhausted:assert counts['e_constructed']==46656
    result=dict(scope='Necessary marked slice for ALL source-genus3 degree12 single3-cycle candidates, NOT double-transposition case',
                counts=dict(counts),exhausted=exhausted,unresolved=unresolved,
                cap=args.cap,elapsed_seconds=time.monotonic()-started,
                enumeration_sha256=digest.hexdigest(),a=a,c=c,p=p,
                word_codes='0=Q1,1=Q1^-1,2=Q2,3=Q2^-1,4=Q3,5=Q3^-1',
                literal_words=known,negative_certificate_complete=exhausted and not unresolved,
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (args.out/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
