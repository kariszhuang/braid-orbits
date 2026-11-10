"""Bounded one-CPU coin-model follow-up; plan-only unless --run is given."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import csv
import hashlib
import json
import math
import random
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from statistics import median

from test_defect_coin_model import ROOT, coins


def decode(text, d):
    perms = []
    for entry in text.strip('[]').split(','):
        p = list(range(d))
        for chunk in re.findall(r'\(([^()]*)\)', entry):
            cyc = [int(x)-1 for x in chunk.split()]
            for a, b in zip(cyc, cyc[1:]+cyc[:1]):
                p[a] = b
        assert sorted(p) == list(range(d))
        perms.append(p)
    assert len(perms) == 4
    return perms


def cycle_type(p):
    seen = set()
    sizes = []
    for i in range(len(p)):
        if i not in seen:
            x, size = i, 0
            while x not in seen:
                seen.add(x); size += 1; x = p[x]
            sizes.append(size)
    return tuple(sorted(sizes, reverse=True))


def passport(t):
    return tuple(sorted(cycle_type(p) for p in t))


def coin_float(g):
    return [float(p) for p in coins(g)]


def compare(hist, g, d, ambient='S'):
    n = sum(hist.values())
    observed = [hist.get(k,0)/n for k in range(d)]
    user = coin_float(g)
    fair = [math.comb(g,k)/2**g for k in range(g+1)]
    degree = coin_float(d-1)
    if ambient == 'A':
        degree = [2*p if k%2 == (d-1)%2 else 0 for k,p in enumerate(degree)]
    def tv(model):
        model = model+[0]*(d-len(model))
        return sum(abs(a-b) for a,b in zip(observed,model))/2
    return {'n': n, 'observed': observed[:g+1], 'coin': user,
            'tv_coin': tv(user), 'tv_fair': tv(fair), 'tv_degree_cycles': tv(degree),
            'mean': sum(k*p for k,p in enumerate(observed)),
            'coin_mean': sum(k*p for k,p in enumerate(user)),
            'even_probability': sum(observed[::2])}


def saved_cases():
    cases, sources, hashes = {}, {}, {}
    for path in (ROOT/'runs/defect/20260928/DEFECT_DISTRIBUTIONS.tsv',
                 ROOT/'runs/defect/20260928_S7_no_identity/DEFECT_DISTRIBUTIONS.tsv'):
        hashes[str(path.relative_to(ROOT))] = sha(path)
        with path.open() as f:
            for r in csv.DictReader(f, delimiter='\t'):
                if not r['case_id'].startswith(('S4:', 'S5:', 'S6:', 'S7:')):
                    continue
                item = cases.setdefault(r['case_id'], {'row': r, 'hist': Counter()})
                item['hist'][int(r['defect'])] += int(r['count'])
    answer = []
    for name,item in cases.items():
        r = item['row']; d,g,n = (int(r[k]) for k in ('degree','genus','states'))
        if g == 0: continue
        path = ROOT/r['source']
        if path not in sources:
            with path.open() as f: sources[path] = {r['component']:r for r in csv.DictReader(f)}
            hashes[str(path.relative_to(ROOT))] = sha(path)
        source = sources[path][name.rsplit('-',1)[1]]
        if '()' in source['seed']: continue
        t = decode(source['seed'],d)
        assert 1-d+sum(d-len(cycle_type(p)) for p in t)//2 == g
        assert n == sum(item['hist'].values()) == int(source['quotient_states'])
        answer.append({'case':name,'d':d,'g':g,'passport':passport(t),'hist':item['hist'],
                       **compare(item['hist'],g,d)})
    assert len(answer) == 1261
    return answer,hashes


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def seed_table():
    candidates = defaultdict(list); hashes = {}
    for shard in range(16):
        path = ROOT/f'runs/exact/S8_shards/shard_{shard}/generating_components.csv'
        hashes[str(path.relative_to(ROOT))] = sha(path)
        with path.open() as f:
            for r in csv.DictReader(f):
                if '()' in r['seed']: continue
                t = decode(r['seed'],8)
                g = 1-8+sum(8-len(cycle_type(p)) for p in t)//2
                if g > 0:
                    candidates[g].append((int(r['quotient_states']),f'S8_s{shard}_c{r["component"]}',t))
    rng = random.Random(2026101203); selected=[]
    for g, rows in sorted(candidates.items()):
        rows.sort(key=lambda x:(-x[0],x[1]))
        picked=rows[:2]+rng.sample(rows[2:],min(3,len(rows)-2))
        for n,name,t in picked:
            selected.append({'id':name,'degree':8,'genus':g,'states':n,'perms':t,
                             'passport':passport(t),'selection':'largest_two_or_seeded_random'})
    return selected,hashes


def analysis(raw, selected):
    cells=defaultdict(Counter); both={}; walkblocks=defaultdict(dict)
    for r in csv.DictReader(raw.splitlines(),delimiter='\t'):
        key=(r['kind'],r['id'],int(r['degree']),r['ambient'],int(r['genus']),int(r['replicate']),int(r['block']))
        cells[key][int(r['defect'])]+=int(r['count'])
        both[key]=int(r['both_pairs_transitive'])
    pooled=defaultdict(Counter); pair_counts=Counter(); reps=defaultdict(list)
    for key,h in cells.items():
        kind,name,d,a,g,rep,block=key;base=(kind,name,d,a,g)
        pooled[base].update(h);pair_counts[base]+=both[key]
        if kind=='iid':reps[base].append(compare(h,g,d,a))
        else:walkblocks[base][(rep,block)] = h
    iid=[];walk=[];selection={s['id']:s for s in selected}
    for key,h in sorted(pooled.items()):
        kind,name,d,a,g=key
        item={'id':name,'degree':d,'ambient':a,'genus':g,**compare(h,g,d,a),
              'both_pairs_transitive_fraction':pair_counts[key]/sum(h.values())}
        if kind=='iid':
            item['replicates']=reps[key];iid.append(item)
        else:
            chains=[];halves=[]
            for rep in range(4):
                ch=Counter()
                for block in range(8):ch.update(walkblocks[key][(rep,block)])
                chains.append(compare(ch,g,d,a))
            for half in range(2):
                hh=Counter()
                for rep in range(4):
                    for block in range(4*half,4*half+4):hh.update(walkblocks[key][(rep,block)])
                halves.append(compare(hh,g,d,a))
            def diff(p,q):return sum(abs(a-b) for a,b in zip(p['observed'],q['observed']))/2
            item.update({'exact_component_states':selection[name]['states'],
                         'max_chain_tv':max(diff(p,q) for p in chains for q in chains),
                         'half_tv':diff(*halves),'chain_summaries':chains,
                         'passport':selection[name]['passport']})
            walk.append(item)
    return iid,walk


def exact_analysis(cases):
    bins=[]
    for g in sorted({c['g'] for c in cases}):
        for low,high in ((0,1000),(1000,10000),(10000,10**20)):
            rows=[c for c in cases if c['g']==g and low<=c['n']<high]
            if rows:bins.append({'genus':g,'size_low':low,'size_high':None if high>10**10 else high,
                                 'components':len(rows),'median_tv':median(c['tv_coin'] for c in rows)})
    groups=defaultdict(list)
    for c in cases:groups[(c['d'],c['g'],c['passport'])].append(c)
    passports=[]
    for (d,g,pp),rows in sorted(groups.items()):
        h=Counter()
        for c in rows:h.update(c['hist'])
        passports.append({'degree':d,'genus':g,'passport':pp,'components':len(rows),
                          'component_tv_min':min(c['tv_coin'] for c in rows),
                          'component_tv_max':max(c['tv_coin'] for c in rows),**compare(h,g,d)})
    fair_wins=sum(c['tv_fair']<c['tv_coin']-1e-12 for c in cases)
    return {'components':len(cases),'fair_binomial_better_components':fair_wins,
            'coin_better_components':sum(c['tv_coin']<c['tv_fair']-1e-12 for c in cases),
            'parity_obstruction_components':sum(abs(c['even_probability']-.5)>=.25 for c in cases),
            'size_bins_by_genus':bins,'passport_cells':passports,
            'cases':[{k:v for k,v in c.items() if k!='hist'} for c in cases]}


def main():
    p=argparse.ArgumentParser();p.add_argument('--run',action='store_true');args=p.parse_args()
    if not args.run:
        print((ROOT/'evidence/defect_coin_followup_20261012/PLAN.md').read_text());return
    out=Path('data')/('defect_coin_followup_'+datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    out.mkdir();cases,hashes=saved_cases();selected,more_hashes=seed_table();hashes.update(more_hashes)
    seedtext='id\tdegree\tpermutations\n'+''.join(s['id']+'\t8\t'+','.join(str(x) for p in s['perms'] for x in p)+'\n' for s in selected)
    (out/'SEEDS.tsv').write_text(seedtext)
    (out/'SELECTION.json').write_text(json.dumps(selected,indent=2))
    # Freeze the executable and every Rust module before bounded acquisition.
    shutil.copytree(ROOT/'src',out/'source'/'src')
    shutil.copy2(ROOT/'Cargo.toml',out/'source'/'Cargo.toml')
    (out/'source'/'scripts').mkdir()
    for name in ('defect_coin_followup.py','test_defect_coin_model.py'):
        shutil.copy2(_scripts_path.find(name),out/'source'/'scripts'/name)
    shutil.copy2(ROOT/'target/release/defect_coin_probe',out/'defect_coin_probe')
    command=[str(out/'defect_coin_probe'),'--trials','100000','--seed','2026101202','--seeds',str(out/'SEEDS.tsv')]
    print(f'One-CPU bounded probe: {out}',flush=True)
    with (out/'HISTOGRAMS.tsv').open('w') as dest:
        subprocess.run(command,stdout=dest,check=True,timeout=300)
    raw=(out/'HISTOGRAMS.tsv').read_text();iid,walk=analysis(raw,selected)
    result={'raw_root':str(out),'command':command,'source_sha256':hashes,
            'binary_sha256':sha(ROOT/'target/release/defect_coin_probe'),
            'code_sha256':{str(f.relative_to(ROOT)):sha(f) for f in [Path(__file__),ROOT/'src/defect_coin_probe.rs',ROOT/'src/bin/defect_coin_probe.rs',ROOT/'src/main.rs',ROOT/'Cargo.toml']},
            'raw_histogram_sha256':sha(out/'HISTOGRAMS.tsv'),'seed_sha256':sha(out/'SEEDS.tsv'),
            'exact':exact_analysis(cases),'iid':iid,'walks':walk,
            'limitations':'iid transitive ambient tuples, not certified full-group generators; walk samples are correlated/non-exhaustive, mixing not proved; no iid walk p-values; no new discovery'}
    (out/'RESULTS.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({'raw_root':str(out),'iid_retained':sum(i['n'] for i in iid),'walks':len(walk),
                      'walk_observations':sum(w['n'] for w in walk),
                      'exact_model_wins':{k:result['exact'][k] for k in ('fair_binomial_better_components','coin_better_components')},
                      'walk_median_tv':median(w['tv_coin'] for w in walk),
                      'walk_max_chain_tv':max(w['max_chain_tv'] for w in walk)},indent=2))


if __name__=='__main__':main()
