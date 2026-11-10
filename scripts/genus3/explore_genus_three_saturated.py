#!/usr/bin/env python3
"""Exact remaining pure-3^5 genus-three passports, one coordinated CPU slot.

Reuses a literal/character-count validated 2800-element completion library.
No geometry rejection. Short-test survivors are saved, NEVER called negative.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
from collections import Counter,deque
import csv
from itertools import combinations,permutations
import json
from math import factorial
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

from audit_guided_native import audit_guided_native
from run_genus_three import sha,write
from verify_genus_three_finite import mul,inv,profile,braid,connection

N=15
IDENTITY=tuple(range(N))
P=tuple((x+1)%N for x in range(N))
ORDINARY=(3,)*5
PROFILES=((3,3),(3,2,2),(2,2,2,2))


def cycles_on(points,lengths):
    if not lengths:
        yield ()
        return
    root=points[0]
    for length in sorted(set(lengths)):
        remaining=list(lengths);remaining.remove(length)
        for rest in combinations(points[1:],length-1):
            support={root,*rest};unused=tuple(x for x in points if x not in support)
            for tail in permutations(rest):
                for later in cycles_on(unused,tuple(remaining)):
                    yield ((root,)+tail,)+later


def complete_frames(special):
    found=set();raw=0;tested=0
    for support in combinations(range(N),sum(special)):
        for cycles in cycles_on(support,special):
            tested+=1;a=list(IDENTITY)
            for cycle in cycles:
                for x,y in zip(cycle,cycle[1:]+cycle[:1]):a[x]=y
            a=tuple(a)
            if profile(mul(inv(a),P))!=ORDINARY:continue
            raw+=1;rotated=[]
            for shift in range(N):
                u=[0]*N
                for x in range(N):u[(x+shift)%N]=(a[x]+shift)%N
                rotated.append(tuple(u))
            assert len(set(rotated))==N,'nontrivial pair centralizer needs different state accounting'
            found.add(min(rotated))
    alpha=(1,)*(N-sum(special))+tuple(sorted(special))
    expected=connection(alpha,ORDINARY)
    assert raw==expected==N*len(found)
    denominator=factorial(N-sum(special))
    for length,count in Counter(special).items():denominator*=length**count*factorial(count)
    assert tested==factorial(N)//denominator
    return sorted(found),dict(class_permutations=tested,full_product_pairs=raw,frames=len(found),character_count=expected,centralizer_order=1)


def audit_frames(frames,special):
    """Coverage by exact characters plus independently literal distinct frames."""
    alpha=tuple(sorted((1,)*(N-sum(special))+special))
    assert len(frames)*N==connection(alpha,ORDINARY)
    keys=set()
    for a in frames:
        assert sorted(a)==list(IDENTITY) and profile(a)==alpha
        assert profile(mul(inv(a),P))==ORDINARY
        rotated=[]
        for shift in range(N):
            shift_p=tuple((x+shift)%N for x in range(N))
            rotated.append(mul(mul(shift_p,a),inv(shift_p)))
        assert len(set(rotated))==N
        key=min(rotated);assert key not in keys;keys.add(key)


def load_completions(path):
    with Path(path).open()as stream:rows=list(csv.DictReader(stream,delimiter='\t'))
    library=[tuple(map(int,row['c'].split(','))) for row in rows]
    assert len(library)==len(set(library))==connection(ORDINARY,ORDINARY)==2800
    assert all(sorted(c)==list(IDENTITY) and profile(c)==ORDINARY and profile(inv(mul(P,c)))==ORDINARY for c in library)
    return library


def is_smooth(t):
    k=N-len(profile(mul(t[0],t[1])));ix=[N-len(profile(p))for p in t]
    return k in (sum(ix[:2]),sum(ix[2:]))


def failure(t,bank):
    for w in bank:
        u=t
        for m in w:u=braid(u,m)
        if not is_smooth(u):return w
    queue=deque([(t,())])
    while queue:
        u,word=queue.popleft()
        if len(word)==4:continue
        for m in range(6):
            if word and (word[-1]+3)%6==m:continue
            v=braid(u,m);w=word+(m,)
            if not is_smooth(v):return w
            queue.append((v,w))
    return None


def notation(p):
    seen=set();cycles=[]
    for root in range(N):
        if root in seen or p[root]==root:continue
        cycle=[];x=root
        while x not in seen:seen.add(x);cycle.append(str(x+1));x=p[x]
        cycles.append('('+' '.join(cycle)+')')
    return ''.join(cycles)or'()'


def run(root):
    library=load_completions(root/'COMPLETIONS.tsv')
    reports=[]
    for special in PROFILES:
        began=time.perf_counter();frames,coverage=complete_frames(special);audit_frames(frames,special)
        folder=root/('special_'+'_'.join(map(str,special)));folder.mkdir()
        write(folder/'FRAMES.json',[list(a)for a in frames]);negative=pending=0;bank=[];hist=Counter()
        with (folder/'LEADS.tsv').open('w')as leads,(folder/'RESOLUTIONS.tsv').open('w')as proofs,(folder/'SURVIVORS.tsv').open('w')as survivors:
            proofs.write('candidate\tstatus\tword\tfolder\tseed\n');survivors.write('frame\tcompletion\tseed\n')
            for frame,a in enumerate(frames):
                b=mul(inv(a),P)
                for number,c in enumerate(library):
                    d=inv(mul(P,c));t=(a,b,c,d);text='['+', '.join(map(notation,t))+']';word=failure(t,bank)
                    if word is None:pending+=1;survivors.write(f'{frame}\t{number}\t{text}\n');continue
                    negative+=1;hist[len(word)]+=1
                    if len(bank)<64 and word not in bank:bank.append(word)
                    leads.write(f'{negative}\t{text}\n');proofs.write(f'{negative}\tnonsmooth\t'+','.join(map(str,word))+f'\t\t{text}\n')
        elapsed=time.perf_counter()-began
        (folder/'COUNTS.tsv').write_text('degree\tlane\tseed\tconstructed\tscope\ttransitive\tpassport\tshort_reject\tduplicate\tclosed\tunresolved\tseconds\n'+f'15\tsolved\t0\t{negative}\t{negative}\t{negative}\t{negative}\t0\t0\t0\t0\t{elapsed}\n')
        # Literal state membership, uniqueness, exact frame/completion counts,
        # and native witness replay form the independent coverage audit.
        allowed_frames=set(frames);allowed_c=set(library);seen=set()
        for file in ('RESOLUTIONS.tsv','SURVIVORS.tsv'):
            from verify_smooth_geometry import parse_tuple
            with (folder/file).open()as stream:
                for row in csv.DictReader(stream,delimiter='\t'):
                    a,b,c,d=tuple(tuple(p)for p in parse_tuple(row['seed'],N))
                    assert a in allowed_frames and c in allowed_c and mul(a,b)==P and d==inv(mul(P,c))
                    assert (a,c)not in seen;seen.add((a,c))
        assert len(seen)==negative+pending==len(frames)*len(library)
        native=audit_guided_native(folder,root/'independent_auditor',native_components=True)
        assert native['independent_native']['records']==negative
        report=dict(special=list(special),coverage=coverage,candidate_states=len(seen),negative=negative,pending=pending,word_length_histogram=dict(hist),producer_seconds=elapsed,native=native,
            status='complete_negative_passport' if not pending else 'needs_full_closure',files_sha256={q.name:sha(q)for q in folder.iterdir()if q.is_file()})
        write(folder/'AUDIT.json',report);reports.append(report)
        print(special,len(seen),negative,pending,flush=True)
    tests=subprocess.run([sys.executable,'-m','unittest','discover','-s',str(Path(__file__).resolve().parent),'-p','test_genus_three_saturated.py','-v'],check=True,capture_output=True,text=True,env=dict(os.environ,G3_SATURATED_ROOT=str(root)))
    write(root/'TESTS.json',dict(passed=True,stdout=tests.stdout,stderr=tests.stderr))
    write(root/'RESULT.json',dict(passed=True,cases=reports,binary_sha256=sha(root/'independent_auditor'),completion_library_sha256=sha(root/'COMPLETIONS.tsv'),tests_sha256=sha(root/'TESTS.json'),
        sources={str(q.relative_to(root/'source')):sha(q)for q in (root/'source').rglob('*.py')},
        scope='three complete pure ordinary-3^5 passports,not all degree15 profiles; same-code frame discovery is checked by independent exact characters and distinct rotation classes; literal negatives independently native-replayed; no geometry exclusion'))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('out',type=Path);p.add_argument('library',type=Path,nargs='?');p.add_argument('auditor',type=Path,nargs='?');p.add_argument('--worker',action='store_true');a=p.parse_args();root=a.out.resolve()
    if a.worker:run(root);return
    assert a.library and a.auditor
    root.mkdir();src=root/'source/scripts';src.mkdir(parents=True)
    for q in Path(__file__).resolve().parents[1].rglob('*.py'):shutil.copyfile(q,src/q.name)
    shutil.copyfile(a.library,root/'COMPLETIONS.tsv');shutil.copyfile(a.auditor,root/'independent_auditor');(root/'independent_auditor').chmod(0o700)
    write(root/'INPUTS.json',dict(library=dict(path=str(a.library.resolve()),sha256=sha(a.library)),auditor=dict(path=str(a.auditor.resolve()),sha256=sha(a.auditor)),cpu_workers=1,gpu=False))
    subprocess.run([sys.executable,str(src/Path(__file__).name),str(root),'--worker'],check=True)


if __name__=='__main__':main()
