"""Bounded, ONE auxiliary CPU; plan-only unless --run. Never alters smooth jobs."""
import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
EVIDENCE=ROOT/'evidence/defect_degree_stability_20261012'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def coins(g):
    p=[1.0]
    for j in range(1,g+1):
        q=[0.0]*(len(p)+1)
        for k,v in enumerate(p):
            q[k]+=v*j/(j+1);q[k+1]+=v/(j+1)
        p=q
    return p


def tv(p,q):
    return sum(abs((p[k] if k<len(p) else 0)-(q[k] if k<len(q) else 0))
               for k in range(max(len(p),len(q))))/2


def prob(h,g):
    n=sum(h.values())
    assert n>0 and all(0<=k<=g for k in h)
    return [h[k]/n for k in range(g+1)]


def schedule():
    rows=[]
    for d in (30,45,59,60,61,80,100,120):
        for a in ('S','A'):
            for rep in range(3):
                rows.append(['iid',d,0,a,rep,20261012100+d*101+rep*19+(a=='A')*10000,80000,0,1])
    for g in (1,2,4,8,16,24,32):
        for d in sorted({max(g+2,2*g-2),2*g,2*g+2,60,80,120}):
            if not 3<=d<=128 or g>=d:continue
            for rep in range(4):
                rows.append(['conditional',d,g,'S',rep,20261012200+d*103+g*1009+rep*29,1250,50000,12])
    for kind,degrees in [('parity_control',(60,61,80)),('smooth_control',(60,90,120))]:
        for d in degrees:
            for rep in range(2):
                rows.append([kind,d,d-1 if kind=='parity_control' else 2,'special',rep,
                             20261012300+d*107+rep*31,1250,5000,8])
    # Randomized order avoids processing time/temperature confounding degree.
    import random
    random.Random(20261012400).shuffle(rows)
    return rows


def write_jobs(path,rows):
    with path.open('w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['kind','degree','genus','ambient','replicate','seed','samples','burn','stride']);w.writerows(rows)


def mul(a,b):return tuple(a[x] for x in b)
def inv(a):
    p=[0]*len(a)
    for x,y in enumerate(a):p[y]=x
    return tuple(p)
def cycles(p):
    seen=set();c=0
    for i in range(len(p)):
        if i not in seen:
            c+=1;x=i
            while x not in seen:seen.add(x);x=p[x]
    return c
def orbits(ps):
    remaining=set(range(len(ps[0])));count=0
    while remaining:
        stack=[remaining.pop()];count+=1
        while stack:
            x=stack.pop()
            for p in ps:
                y=p[x]
                if y in remaining:remaining.remove(y);stack.append(y)
    return count
def genus(t):
    d=len(t[0]);total=sum(d-cycles(p) for p in t);assert total%2==0
    return 1-d+total//2
def braid(t,m):
    u=list(t);j=m//2;a,b=t[j:j+2]
    u[j:j+2]=[b,mul(mul(inv(b),a),b)] if m%2==0 else [mul(mul(a,b),inv(a)),a]
    return tuple(u)


def audit_tuples(out,jobs):
    checked=0
    for row in csv.DictReader((out/'TUPLES.tsv').open(),delimiter='\t'):
        job=int(row['job']);kind,d,g,*_=jobs[job];flat=list(map(int,row['permutations'].split(',')))
        assert len(flat)==4*d
        t=tuple(tuple(flat[i*d:(i+1)*d]) for i in range(4))
        identity=tuple(range(d))
        assert all(sorted(p)==list(identity) and p!=identity for p in t)
        assert mul(mul(mul(*t[:2]),t[2]),t[3])==identity and orbits(t[:3])==1
        actual_g=genus(t);c=cycles(mul(*t[:2]));l=orbits(t[:2]);r=orbits(t[2:]);delta=c-l-r+1
        assert 0<=delta<=actual_g
        if kind!='iid':assert actual_g==g
        for m in range(6):
            u=braid(t,m);assert braid(u,m^1)==t and genus(u)==actual_g
        if kind=='smooth_control':
            assert delta==0
            for m in range(6):
                u=braid(t,m);assert cycles(mul(*u[:2]))-orbits(u[:2])-orbits(u[2:])+1==0
        if kind=='parity_control':assert all(cycles(p)==1 for p in t) and delta%2==(d-1)%2
        checked+=1
    return checked


def analyze(out):
    pooled=defaultdict(Counter);chains=defaultdict(lambda:defaultdict(Counter));halves=defaultdict(lambda:defaultdict(Counter));joint=defaultdict(Counter)
    for r in csv.DictReader((out/'HISTOGRAMS.tsv').open(),delimiter='\t'):
        kind,d,a,g=r['kind'],int(r['degree']),r['ambient'],int(r['genus']);rep,block=int(r['replicate']),int(r['block'])
        c,l,rr,n=[int(r[k]) for k in ('cycles','left','right','count')];delta=c-l-rr+1;assert 0<=delta<=g
        key=(kind,d,a,g);pooled[key][delta]+=n;chains[key][rep][delta]+=n;halves[key][block//4][delta]+=n;joint[key][(c,l,rr)]+=n
    cells=[]
    for key,h in sorted(pooled.items()):
        kind,d,a,g=key;p=prob(h,g);q=coins(g);rep_probs=[prob(v,g) for v in chains[key].values()]
        cell=dict(kind=kind,degree=d,ambient=a,genus=g,n=sum(h.values()),observed=p,coin=q,
                  tv_coin=tv(p,q),mean=sum(k*v for k,v in enumerate(p)),coin_mean=sum(k*v for k,v in enumerate(q)),
                  replicate_tv_coin=[tv(v,q) for v in rep_probs],
                  max_replicate_tv=max(tv(x,y) for x in rep_probs for y in rep_probs),
                  even_probability=sum(p[::2]),both_pairs_transitive_fraction=sum(n for(c,l,r),n in joint[key].items() if l==r==1)/sum(h.values()),
                  joint=[dict(cycles=c,left=l,right=r,count=n) for(c,l,r),n in sorted(joint[key].items())])
        if kind!='iid':cell['half_tv']=tv(*[prob(halves[key][v],g) for v in (0,1)])
        cells.append(cell)
    cross=[]
    for g in sorted({c['genus'] for c in cells if c['kind']=='conditional'}):
        cc=[c for c in cells if c['kind']=='conditional' and c['genus']==g]
        ref=max(cc,key=lambda c:c['degree'])
        for c in cc:cross.append(dict(genus=g,degree=c['degree'],reference_degree=ref['degree'],tv_to_reference=tv(c['observed'],ref['observed']),
                                      above_2g=c['degree']>=2*g,coin_tv=c['tv_coin'],max_replicate_tv=c['max_replicate_tv'],half_tv=c['half_tv']))
    metrics=list(csv.DictReader((out/'METRICS.tsv').open(),delimiter='\t'))
    return dict(cells=cells,cross_degree=cross,metrics=metrics,
                limitations='iid labeled transitive ambient starts; conditional and braid samples correlated, large-degree mixing/irreducibility unproved; not full-group certification, orbit census or all-g proof')


def render_report(result):
    lines=['# Fixed-genus/large-degree defect pilot','',
           'TV is half the summed absolute probability difference. No universal verification is possible from this finite pilot. All cells, including sparse IID strata, are in RESULTS.json.',
           '', '## Independent labeled ambient starts (modal genus shown; ALL strata saved)','',
           '| d | ambient | modal g | n in stratum | coin TV | max replicate TV |','|---:|:---:|---:|---:|---:|---:|']
    iid=[c for c in result['cells'] if c['kind']=='iid']
    for d,a in sorted({(c['degree'],c['ambient']) for c in iid}):
        c=max((c for c in iid if (c['degree'],c['ambient'])==(d,a)),key=lambda c:c['n'])
        lines.append(f"|{d}|{a}|{c['genus']}|{c['n']}|{c['tv_coin']:.4f}|{c['max_replicate_tv']:.4f}|")
    lines+=['','## Fixed-genus conditional-MCMC pilot','',
            'Uniform conditional labeled tuples are stationary on each connected proposal class. Unknown large-degree connectivity/mixing means these estimates are exploratory, not an exact conditional census. Replicate and half discrepancies are diagnostics, NOT mixing proofs.',
            '', '| g | d | coin TV | TV to d=120 | max chain TV | half TV |','|---:|---:|---:|---:|---:|---:|']
    for c in result['cross_degree']:
        lines.append(f"|{c['genus']}|{c['degree']}|{c['coin_tv']:.4f}|{c['tv_to_reference']:.4f}|{c['max_replicate_tv']:.4f}|{c['half_tv']:.4f}|")
    lines+=['','## Special-passport braid controls','', '| family | d | g | coin TV |', '|---|---:|---:|---:|']
    for c in result['cells']:
        if c['kind'].endswith('control'):lines.append(f"|{c['kind']}|{c['degree']}|{c['genus']}|{c['tv_coin']:.4f}|")
    lines+=['','No full symmetric/alternating generation, new smooth component, geometric novelty, whole-degree coverage, or universal asymptotic theorem is claimed.',
            '', 'Design/statistical methods credit: Timothy Kassis, Vinayak Agarwal, Yuhuan He, Darshil Patel, Aubrey M. Brueckner (2026), *Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents*, [arXiv:2609.00065](https://doi.org/10.48550/arXiv.2609.00065); current record checked October12. This is not mathematical evidence.']
    return '\n'.join(lines)+'\n'


def main():
    p=argparse.ArgumentParser();p.add_argument('--run',action='store_true');p.add_argument('--confirm',type=Path);args=p.parse_args()
    if not args.run:print((EVIDENCE/'PLAN.md').read_text());return
    base=Path('data')
    if not Path('data').is_dir() or shutil.disk_usage(base).free<30*2**30:raise RuntimeError('mounted T7 with30GiB required')
    out=base/('defect_degree_stability_'+datetime.now().strftime('%Y%m%d_%H%M%S_%f'));out.mkdir();source=out/'source';(source/'src/bin').mkdir(parents=True);(source/'scripts/spectra').mkdir(parents=True)
    paths=['Cargo.toml','Cargo.lock','src/bin/defect_degree_probe.rs','scripts/spectra/defect_degree_stability.py','scripts/spectra/test_defect_degree_stability.py','evidence/defect_degree_stability_20261012/PLAN.md']
    hashes={}
    for name in paths:
        dest=source/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,dest);hashes[name]=sha(dest)
    # Standalone dependency-free target; no production default-run/target changes.
    binary=out/'defect_degree_probe';test_binary=out/'probe_tests'
    try:
        with (out/'BUILD.log').open('w') as log:
            rust_source=str(source/'src/bin/defect_degree_probe.rs')
            subprocess.run(['rustc','--edition=2021','--test',rust_source,'-o',str(test_binary)],stdout=log,stderr=subprocess.STDOUT,check=True)
            subprocess.run([str(test_binary),'--test-threads=1'],stdout=log,stderr=subprocess.STDOUT,check=True)
            subprocess.run(['rustc','--edition=2021','-O','-C','target-cpu=native','-C','codegen-units=1',rust_source,'-o',str(binary)],stdout=log,stderr=subprocess.STDOUT,check=True)
    except BaseException as e:
        (out/'FAILURE.json').write_text(json.dumps(dict(error=repr(e),phase='build',partial_data_preserved=True)));raise
    subprocess.run([sys.executable,str(source/'scripts/spectra/test_defect_degree_stability.py'),'--binary',str(binary)],check=True)
    if args.confirm:
        if json.loads((args.confirm/'FINISHED.json').read_text())['status']!='complete':raise ValueError('confirmation requires completed pilot')
        jobs=[['conditional',d,g,'S',rep,20261012500+d*109+g*1013+rep*41,12500,500000,24]
              for g in (4,8,16) for d in (60,120) for rep in range(4)]
    else:jobs=schedule()
    write_jobs(out/'JOBS.tsv',jobs);print(f'OUTPUT={out} JOBS={len(jobs)} one auxiliary CPU; smooth unchanged',flush=True)
    command=[str(binary),str(out/'JOBS.tsv'),str(out/'TUPLES.tsv'),str(out/'METRICS.tsv')]
    (out/'PLAN.json').write_text(json.dumps(dict(jobs=len(jobs),source_sha256=hashes,binary_sha256=sha(binary),command=command,timeout_seconds=600,
          exploratory_confirmation_of=str(args.confirm) if args.confirm else None,
          confirmation_rationale='adaptive replication of stable-but-mismatching genera4,8,16,10x burn/sample count,2x stride; new seeds; not preplanned confirmatory significance' if args.confirm else None),indent=2))
    try:
        with (out/'HISTOGRAMS.tsv').open('w') as h,(out/'PROGRESS.log').open('w') as log:subprocess.run(command,stdout=h,stderr=log,check=True,timeout=600)
        result=analyze(out);result['independent_literal_tuples_checked']=audit_tuples(out,jobs)
        result['raw_root']=str(out);result['raw_sha256']={n:sha(out/n) for n in ('JOBS.tsv','TUPLES.tsv','METRICS.tsv','HISTOGRAMS.tsv')}
        if args.confirm:
            old=json.loads((args.confirm/'RESULTS.json').read_text())
            result['pilot_comparison']=[]
            for c in result['cells']:
                previous=next(v for v in old['cells'] if all(c[k]==v[k] for k in ('kind','degree','ambient','genus')))
                result['pilot_comparison'].append(dict(degree=c['degree'],genus=c['genus'],tv_to_pilot=tv(c['observed'],previous['observed']),pilot_coin_tv=previous['tv_coin'],long_coin_tv=c['tv_coin']))
        (out/'RESULTS.json').write_text(json.dumps(result,indent=2));(out/'REPORT.md').write_text(render_report(result))
        (out/'FINISHED.json').write_text(json.dumps(dict(status='complete',jobs=len(jobs),result_sha256=sha(out/'RESULTS.json')),indent=2))
        print(f'COMPLETE={out} tuple_checks={result["independent_literal_tuples_checked"]}',flush=True)
    except BaseException as e:
        (out/'FAILURE.json').write_text(json.dumps(dict(error=repr(e),partial_data_preserved=True)));raise


if __name__=='__main__':main()
