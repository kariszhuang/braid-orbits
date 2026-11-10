#!/usr/bin/env python3
"""Independent S13 support-conjugacy and RGS-colouring coverage audit."""
import argparse
from itertools import combinations, permutations, product
from collections import Counter
import hashlib
import json
from pathlib import Path
import time
N=13
I=tuple(range(N))
def mul(a,b):return tuple(a[b[i]]for i in range(N))
def inverse(a):
    r=[0]*N
    for i,j in enumerate(a):r[j]=i
    return tuple(r)
def cycles(a):
    used=set();out=[]
    for i in range(len(a)):
        if i in used:continue
        c=[]
        while i not in used:used.add(i);c.append(i);i=a[i]
        out.append(c)
    return out
def ix(a):return N-len(cycles(a))
def canon(a):return min(tuple((a[(i-k)%N]+k)%N for i in range(N))for k in range(N))
def step(t,code):
    j,neg=divmod(code,2);a,b=t[j:j+2]
    pair=(b,mul(mul(inverse(b),a),b))if neg else(mul(mul(a,b),inverse(a)),a)
    return t[:j]+pair+t[j+2:]
def bad(t,w):
    for code in w:t=step(t,code)
    r=ix(mul(t[0],t[1]));return r!=ix(t[0])+ix(t[1])and r!=ix(t[2])+ix(t[3])
def labels_rec(word,counts,q,profile,out):
    if len(word)==N:
        if sorted(counts)==profile:out.append(tuple(word))
        return
    maxlen=max(profile)
    for j in range(min(len(counts)+1,q)):
        if j==len(counts):counts.append(1);word.append(j);labels_rec(word,counts,q,profile,out);word.pop();counts.pop()
        elif counts[j]<maxlen:
            counts[j]+=1;word.append(j);labels_rec(word,counts,q,profile,out);word.pop();counts[j]-=1
def allowed(labels,acs):
    for k in range(N):
        # Independent graph test by explicit component sets, not producer DSU.
        components=[{j}for j in range(max(labels)+1)]
        for cy in acs:
            for point in cy[1:]:
                x,y=labels[(cy[0]+k)%N],labels[(point+k)%N]
                cx=next(c for c in components if x in c);dy=next(c for c in components if y in c)
                if cx is dy:return False
                cx.update(dy);components.remove(dy)
    return True
def main():
    p=argparse.ArgumentParser();p.add_argument('kind',choices=['double','three-two','triple']);p.add_argument('producer',type=Path);p.add_argument('output',type=Path);a=p.parse_args();assert not a.output.exists();a.output.mkdir(parents=True);start=time.monotonic()
    shape={'double':[2,2],'three-two':[2,3],'triple':[2,2,2]}[a.kind];s=sum(v-1 for v in shape);profile=[4,4,5]if s==2 else[3,3,3,4];support_size=sum(shape)
    # Enumerate the complete small active class first, then all support subsets.
    local=[v for v in permutations(range(support_size))if sorted(map(len,cycles(v)))==shape]
    target=tuple(range(1,N))+(0,);frames=set();counts=Counter()
    for support in combinations(range(N),support_size):
        for v in local:
            first=list(I)
            for i,j in enumerate(v):first[support[i]]=support[j]
            first=tuple(first);counts['all_a']+=1;second=mul(inverse(first),target);cs=cycles(second)
            if len(cs)!=s+1:continue
            counts['full_labelled_frames']+=1
            if sorted(map(len,cs))!=profile:continue
            counts['bounded_labelled_frames']+=1;frames.add(canon(first))
    producer=json.loads((a.producer/'RESULT.json').read_text());assert len(frames)==producer['rotation_frames'];assert frames=={tuple(v['a'])for v in producer['frames']}
    labels=[];labels_rec([0],[1],s+1,profile,labels);assert len(labels)==producer['support_partitions']
    raw=(a.producer/'ENUMERATION.bin').read_bytes();assert len(raw)%(2*N)==0;expected={raw[i:i+2*N]for i in range(0,len(raw),2*N)};assert len(expected)*2*N==len(raw);seen=set()
    for first in sorted(frames):
        second=mul(inverse(first),target);acs=[c for c in cycles(first)if len(c)>1];twists=[tuple((first[(i-k)%N]+k)%N for i in range(N))for k in range(N)]
        for lab in labels:
            if not allowed(lab,acs):continue
            counts['accepted_supports']+=1;blocks=[tuple(i for i in range(N)if lab[i]==j)for j in range(s+1)];os=[[(b[0],)+v for v in permutations(b[1:])]for b in blocks]
            for cys in product(*os):
                e=list(I)
                for c in cys:
                    for i in range(len(c)):e[c[i]]=c[(i+1)%len(c)]
                e=tuple(e);record=bytes(first)+bytes(e);assert record in expected and record not in seen;seen.add(record);counts['generated_e']+=1;b=mul(inverse(target),inverse(e))
                if ix(b)!=6+s:counts['wrong_b_index']+=1;continue
                t=(first,second,b,e);assert mul(mul(mul(first,second),b),e)==I
                for k,twist in enumerate(twists):
                    if ix(mul(twist,b))!=6+2*s:assert bad(t,[0]*(2*k)+[3]);counts['literal_negative']+=1;break
                else:counts['unresolved']+=1
    assert seen==expected
    for key in ['accepted_supports','generated_e','wrong_b_index','literal_negative','unresolved']:assert counts[key]==producer[key]
    result=dict(passed=True,kind=a.kind,counts=dict(counts),rotation_frames=len(frames),support_partitions=len(labels),sorted_record_sha256=hashlib.sha256(b''.join(sorted(seen))).hexdigest(),producer_binary_sha256=hashlib.sha256(raw).hexdigest(),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),elapsed_seconds=time.monotonic()-start)
    (a.output/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
