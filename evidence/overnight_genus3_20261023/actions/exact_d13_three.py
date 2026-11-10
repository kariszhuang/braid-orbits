#!/usr/bin/env python3
"""Single3-cycle in the unsaturated S13 genus3 budget; exact necessary slice."""
import argparse
from collections import Counter
from itertools import permutations, product
import hashlib
import json
from pathlib import Path
import time
N=13
I=tuple(range(N))
def mul(a,b):return tuple(a[b[i]]for i in range(N))
def inv(a):
    b=[0]*N
    for i,j in enumerate(a):b[j]=i
    return tuple(b)
def cycles(a):
    used=set();out=[]
    for i in range(N):
        if i in used:continue
        c=[]
        while i not in used:used.add(i);c.append(i);i=a[i]
        out.append(c)
    return out
def ix(a):return N-len(cycles(a))
def braid(t,code):
    j,neg=divmod(code,2);a,b=t[j:j+2]
    pair=(b,mul(mul(inv(b),a),b))if neg else(mul(mul(a,b),inv(a)),a)
    return t[:j]+pair+t[j+2:]
def bad(t,word):
    for code in word:t=braid(t,code)
    r=ix(mul(t[0],t[1]));return r!=ix(t[0])+ix(t[1])and r!=ix(t[2])+ix(t[3])
def main():
    ap=argparse.ArgumentParser();ap.add_argument('out',type=Path);args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    p=tuple(range(1,N))+(0,);a=list(I);a[0],a[4],a[8]=4,8,0;a=tuple(a);c=mul(inv(a),p);assert sorted(map(len,cycles(c)))==[4,4,5]
    counts=Counter();digest=hashlib.sha256();unresolved=[];words=[]
    # A fixed point removes one vertex of C13 squared in the p^4 order.
    # The other12 colours repeat with period3; cycle orientations stay free.
    with (args.out/'WITNESSES.jsonl').open('w')as output, (args.out/'ENUMERATION.bin').open('wb')as binary:
        for fixed in range(N):
            blocks=[tuple((fixed+4*j)%N for j in range(1,N)if(j-1)%3==colour)for colour in range(3)]
            os=[[(b[0],)+v for v in permutations(b[1:])]for b in blocks]
            for cs in product(*os):
                e=list(I)
                for cy in cs:
                    for j in range(4):e[cy[j]]=cy[(j+1)%4]
                e=tuple(e);digest.update(bytes(e));binary.write(bytes(e));counts['generated_e']+=1;b=mul(inv(p),inv(e))
                if ix(b)!=9:counts['wrong_b_index']+=1;continue
                t=(a,c,b,e);assert mul(mul(mul(a,c),b),e)==I
                for word in words:
                    if bad(t,word):break
                else:
                    word=None
                    for length in range(1,5):
                        for w in product(range(6),repeat=length):
                            if bad(t,w):word=w;break
                        if word is not None:break
                    if word is not None:words.append(word)
                if word is None:unresolved.append(t);counts['unresolved']+=1
                else:
                    assert bad(t,word);counts['literal_negative']+=1
                    output.write(json.dumps(dict(tuple=t,word=word))+'\n')
    assert counts['generated_e']==2808
    result=dict(scope='Complete necessary single3-cycle degree13 indexbudget2,9,9,10, not double-transposition',counts=dict(counts),unresolved=unresolved,enumeration_sha256=digest.hexdigest(),words=words,maximum_word_length=max(map(len,words),default=0),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),elapsed_seconds=time.monotonic()-start)
    (args.out/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
