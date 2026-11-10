#!/usr/bin/env python3
"""Independent RGS-colouring coverage, full tuple reconstruction and negatives."""
import argparse
from collections import Counter
from itertools import combinations, permutations, product
import hashlib
import json
from pathlib import Path
import time
N=13
I=tuple(range(N))
def mul(a,b):return tuple(a[b[i]] for i in range(N))
def inverse(a):
    b=[0]*N
    for i,j in enumerate(a):b[j]=i
    return tuple(b)
def cycles(a):
    used=set();cs=[]
    for i in range(N):
        if i in used:continue
        c=[]
        while i not in used:used.add(i);c.append(i);i=a[i]
        cs.append(c)
    return cs
def ix(a):return N-len(cycles(a))
def canonical_edges(es):
    # The frozen Rust file minimizes the permutation array, not the edge list.
    # Both mark the same rotation orbits, but raw tuple-set comparison requires
    # the SAME marking. Keep this audit's independent edge-based enumeration.
    rotations=(tuple(sorted(tuple(sorted(((x+k)%N,(y+k)%N)))for x,y in es))for k in range(N))
    return min(rotations,key=match)
def match(es):
    a=list(I)
    for x,y in es:a[x],a[y]=y,x
    return tuple(a)
def rgs(word,top,sizes,out):
    if len(word)==N:
        if top==3 and min(sizes)>=1:out.append(tuple(word))
        return
    for c in range(min(top+1,3)+1):
        if sizes[c]>=6:continue
        sizes[c]+=1;word.append(c);rgs(word,max(top,c),sizes,out);word.pop();sizes[c]-=1
def allowed(labels,es):
    for k in range(N):
        u,v=es[0];x,y=es[1];a,b,c,d=(labels[(t+k)%N]for t in (u,v,x,y))
        if a==b or c==d or {a,b}=={c,d}:return False
    return True
def step(t,code):
    j,neg=divmod(code,2);a,b=t[j:j+2]
    pair=(b,mul(mul(inverse(b),a),b))if neg else(mul(mul(a,b),inverse(a)),a)
    return t[:j]+pair+t[j+2:]
def bad(t,w):
    for code in w:t=step(t,code)
    r=ix(mul(t[0],t[1]));return r!=ix(t[0])+ix(t[1])and r!=ix(t[2])+ix(t[3])
def main():
    ap=argparse.ArgumentParser();ap.add_argument('producer',type=Path);ap.add_argument('output',type=Path);a=ap.parse_args();assert not a.output.exists();a.output.mkdir(parents=True);started=time.monotonic()
    p=tuple(range(1,N))+(0,);pi=inverse(p);frames=set();all_count=full=bounded=0
    for es in combinations(combinations(range(N),2),2):
        if len(set(es[0]+es[1]))!=4:continue
        all_count+=1;c=mul(match(es),p);cc=cycles(c)
        if len(cc)!=3:continue
        full+=1
        if not all(len(v)>=3 for v in cc):continue
        bounded+=1;frames.add(canonical_edges(es))
    assert(all_count,full,bounded,len(frames))==(2145,1430,325,25)
    labels=[];rgs([0],0,[1,0,0,0],labels);assert len(labels)==2341339
    raw=a.producer.joinpath('ENUMERATION.bin').read_bytes();record_size=2*N;assert len(raw)%record_size==0
    produced={raw[i:i+record_size]for i in range(0,len(raw),record_size)};assert len(produced)*record_size==len(raw)
    counts=Counter();words=([3,3],[2,2],[0,3,3],[2,1,2],[1,2,2],[0,0,3,3],[2,2,2],[3,3,3],[3,0,3],[0,3,0,3],[1,2,1,2],[0,3,3,3],[1,2,2,2],[3,0,0,3],[2,1,1,2]);pending=[];seen=set()
    for es in sorted(frames):
        first=match(es);second=mul(first,p)
        twists=[match(tuple(((x+k)%N,(y+k)%N)for x,y in es))for k in range(N)]
        for lab in labels:
            if not allowed(lab,es):continue
            counts['supports']+=1;blocks=[tuple(i for i in range(N)if lab[i]==c)for c in range(4)]
            orders=[[(b[0],)+z for z in permutations(b[1:])]for b in blocks]
            for cs in product(*orders):
                e=list(I)
                for c in cs:
                    for i in range(len(c)):e[c[i]]=c[(i+1)%len(c)]
                e=tuple(e);record=bytes(first)+bytes(e);assert record in produced and record not in seen;seen.add(record);counts['generated']+=1
                b=mul(pi,inverse(e))
                if ix(b)!=9:counts['wrong_b_index']+=1;continue
                t=(first,second,b,e);assert mul(mul(mul(first,second),b),e)==I
                for k,twist in enumerate(twists):
                    if ix(mul(twist,b))!=11:
                        assert bad(t,[0]*(2*k)+[3]);counts['twist_negatives']+=1;break
                else:
                    counts['late_starts']+=1
                    for w in words:
                        if bad(t,w):counts['late_literal_negatives']+=1;break
                    else:pending.append(t)
    assert seen==produced
    report=dict(scope='ALL degree13 genus3 remaining double-transposition budget2,9,9,10, using proved support bounds',a_class_labelled=all_count,full_pair_frames_labelled=full,bounded_frames_labelled=bounded,rotation_classes=len(frames),counts=dict(counts),unresolved=len(pending),pending=pending,sorted_enumeration_sha256=hashlib.sha256(b''.join(sorted(seen))).hexdigest(),producer_enumeration_sha256=hashlib.sha256(raw).hexdigest(),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),elapsed_seconds=time.monotonic()-started)
    (a.output/'RESULT.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
