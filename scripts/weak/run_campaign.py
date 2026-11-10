#!/usr/bin/env python3
"""Build, verify, and run resumable research campaigns. Stdlib only."""
from __future__ import annotations

import argparse
import csv
from collections import deque, OrderedDict
from contextlib import contextmanager, ExitStack
from functools import lru_cache
import datetime as dt
import fcntl
import hashlib
import heapq
import itertools
import json
import math
import os
from pathlib import Path
import random
import platform
import traceback
import signal
import re
import shutil
import subprocess
import sys
import time
from zoneinfo import ZoneInfo

PROJECT = Path(__file__).resolve().parents[2]
LEGACY_PROJECT = Path('legacy-project')  # old project folder; only used to rewrite saved paths
RUNS = PROJECT / "runs"
RESULTS = PROJECT / 'data'  # search output; may be a symlink to a larger disk
CAMPAIGNS = RUNS / "campaigns"
CONTROL = PROJECT / '.run-control'
LATEST = CONTROL / "LATEST.json"
CHILD = None
PAIR_CHILDREN = set()
STOPPING = False
CORE_TIMEZONE = ZoneInfo("America/Chicago")


def relocated_path(value):
    """Resolve only this project's historical root; never rewrite saved evidence."""
    text=str(value);legacy=str(LEGACY_PROJECT)
    if text==legacy or text.startswith(legacy+'/'):
        target=PROJECT/Path(text).relative_to(LEGACY_PROJECT)
        return str(target)
    return text


def require_results_volume(root):
    if Path(root).resolve().is_relative_to(RESULTS.resolve()) and not RESULTS.is_dir():
        raise RuntimeError('data/ does not exist; create it (or symlink it to a larger disk) before storing results')


def free_result_space(root):
    require_results_volume(root)
    return shutil.disk_usage(root).free


def run_event(root,event,**fields):
    record=dict(time_utc=dt.datetime.now(dt.timezone.utc).isoformat(),event=event,**fields)
    with (root/'EVENTS.jsonl').open('a') as stream:
        stream.write(json.dumps(record,sort_keys=True)+'\n')
        stream.flush();os.fsync(stream.fileno())


def scheduled_threads(now=None):
    now = now or dt.datetime.now(CORE_TIMEZONE)
    return 4 if 1 <= now.hour < 16 else 8


def seconds_until_thread_change(now=None):
    now = now or dt.datetime.now(CORE_TIMEZONE)
    if now.hour < 1:
        change = now.replace(hour=1, minute=0, second=0, microsecond=0)
    elif now.hour < 16:
        change = now.replace(hour=16, minute=0, second=0, microsecond=0)
    else:
        change = (now + dt.timedelta(days=1)).replace(hour=1, minute=0, second=0, microsecond=0)
    return max(1, math.ceil((change-now).total_seconds()))


def certificate_header(path):
    # Certificates can contain millions of states; reports only need the header.
    header = {}
    orders = []
    with path.open() as stream:
        for line in stream:
            if line.startswith("STATE\t"):
                break
            if line.startswith("ORDER\t"):
                orders.append(int(line.split()[1]))
            elif "=" in line:
                key, value = line.rstrip().split("=", 1)
                header[key] = value
    header["max_pair_order"] = max(orders, default=0)
    return header


def seed_profile(seed, degree):
    types = []
    total_index = 0
    for part in seed.strip()[1:-1].split(','):
        lengths = [len(c.split()) for c in re.findall(r"\(([^()]*)\)", part) if c.strip()]
        total_index += sum(n-1 for n in lengths)
        types.append('.'.join(map(str, sorted(n for n in lengths if n>1))) or '1')
    if len(types) != 4 or total_index % 2 or total_index < 2*degree-2:
        raise ValueError("Invalid generating tuple profile")
    return '/'.join(sorted(types)), 1-degree+total_index//2


def certificate_identity(path):
    """Exact component identity, independent of certificate filename/seed."""
    header = certificate_header(path)
    representative = header.get('representative')
    if representative:
        key = tuple(int(word, 16) for word in representative.split('/'))
    else:
        # Old certificates lack this header. Stream their state rows once;
        # do not infer discovery from a possibly arbitrary filename.
        with path.open() as stream:
            key = min(tuple(int(word, 16) for word in line.split()[1:])
                      for line in stream if line.startswith('STATE\t'))
    return header['group'], int(header['d']), key


def research_jobs(degrees, known, round_no, resolved_a9=()):
    # This is a discovery-priority heuristic, not a measured success probability.
    resolved_a9 = {'/'.join(sorted(p.split('/'))) for p in resolved_a9}
    rows = []
    for d in degrees:
        for group in ('A', 'S'):
            if group == 'A' and d == 8:
                continue
            eligible = any(int(h['d']) <= d <= int(h['d'])+2 and
                           (group == 'S' or h['group'] == 'A') for h in known)
            profiles = {seed_profile(h['seed'],d)[0] for h in known
                        if h['group'] == group and int(h['d']) == d}
            if group == 'A' and d >= 9:
                p = '.'.join(['3']*(d//3))
                profiles.add('/'.join([p]*4))
                # Every order-three genus-zero passport: if entry i consists
                # of k_i disjoint 3-cycles, Riemann-Hurwitz gives sum k_i=d-1.
                # Enumerating these few profiles is exact; tuples inside each
                # profile are still sampled, so no census claim follows.
                for counts in itertools.combinations_with_replacement(range(1,d//3+1),4):
                    if sum(counts)==d-1:
                        profiles.add('/'.join('.'.join(['3']*k) for k in counts))
            if group == 'A' and d == 9:
                profiles -= resolved_a9
            if eligible:
                rows.append((group,d,'deform',False,None))
            rows.extend((group,d,family,weak,None) for family,weak in
                        [('low-index',False),('low-index',True),('intransitive',True),('mixed',False)])
            # Keep a uniform control periodically; its throughput is not the
            # metric used to judge targeted discovery yield.
            if round_no % 4 == 0:
                rows.append((group,d,'uniform',False,None))
            for passport in sorted(profiles):
                rows.append((group,d,'passport',False,passport))
                rows.append((group,d,'passport',True,passport))
    def priority(row):
        group,d,family,weak,passport = row
        indices = [sum(int(n)-1 for n in p.split('.')) for p in passport.split('/')] if passport else [0]
        genus = 1-d+sum(indices)//2 if passport else 0
        focused = d <= 11 or (group == 'A' and d == 13 and passport == '/'.join(['3.3.3']*4))
        stage = {'deform':0,'passport':1 if genus==0 else 3,
                 'intransitive':2,'low-index':4,'mixed':5,'uniform':6}[family]
        degree_order = {10:0,11:1,9:2,8:3,13:4}.get(d,5+degrees.index(d))
        return (not focused,stage,degree_order,group,genus,max(indices)-min(indices),passport or '',weak)
    return sorted(rows,key=priority)


def class_size(cycle_type, degree):
    lengths=[int(n) for n in cycle_type.split('.')]
    lengths += [1]*(degree-sum(lengths))
    denominator=1
    for n in set(lengths):
        count=lengths.count(n)
        denominator *= n**count*math.factorial(count)
    return math.factorial(degree)//denominator


def genus_zero_profiles(degree):
    def partitions(n, smallest=1):
        if n==0:
            yield ()
        for k in range(smallest,n+1):
            for tail in partitions(n-k,k):
                yield (k,)+tail
    # An identity entry makes the strong statement immediate: the other three
    # have product one, so any two of them generate the full ambient group.
    types=[p for p in partitions(degree) if 0<degree-len(p) and (degree-len(p))%2==0]
    profiles=set()
    for four in itertools.combinations_with_replacement(types,4):
        if sum(degree-len(p) for p in four)==2*degree-2:
            profiles.add('/'.join('.'.join(str(n) for n in p if n>1) for p in four))
    return profiles


def exact_slices(known, all_genus_zero=False):
    profiles={row[4] for row in research_jobs([9],known,0) if row[0]=='A' and row[2]=='passport'}
    if all_genus_zero:
        profiles.update(genus_zero_profiles(9))
    profiles={'/'.join(sorted(p.split('/'))) for p in profiles}
    slices=[]
    for profile in profiles:
        # Braid moves can reorder the four cycle types. Put the two smallest
        # classes in positions 2/3, the only entries we enumerate independently.
        types=sorted(profile.split('/'),key=lambda p:(class_size(p,9),p))
        ordered='/'.join([types[2],types[0],types[1],types[3]])
        total=class_size(types[0],9)*class_size(types[1],9)
        if total<= (60_000_000 if all_genus_zero else 20_000_000):
            slices.append((ordered,total))
        elif all_genus_zero:
            raise ValueError(f'Genus-zero slice exceeds implemented budget: {ordered}; no silent omission')
    return sorted(slices,key=lambda row:(row[1],row[0]))


def exact_priority(item):
    # Keep stored slice IDs stable: change execution order, never renumber work.
    _, (passport,total) = item
    structured = set(passport.split('/')) <= {'3','3.3','3.3.3','2.2.2.2'}
    return (not structured,total,passport)


def degree11_involution_slices():
    """Remaining genus-zero equal-involution passports; not all A11/S11."""
    rows = []
    for group, m in [('A', 4), ('S', 5)]:
        involution = '.'.join(['2'] * m)
        cycles = 3*m - 11 + 2
        for lengths in itertools.combinations_with_replacement(range(1, 12), cycles):
            if sum(lengths) != 11:
                continue
            # Already completed, with zero unresolved cases: see C026.
            if group == 'A' and lengths in ((3, 4, 4), (3, 3, 5)):
                continue
            fourth = '.'.join(str(n) for n in lengths if n > 1)
            rows.append((group, '/'.join([involution]*3 + [fourth]), class_size(involution, 11)**2))
    return rows


def braid_fixed_jobs():
    """Biased symmetry construction, including positive genus; NOT a census."""
    rows = [('A', 25, 10, 0)]  # known example lies in this construction
    for d in range(4, 26):
        for genus in (0, 1):
            if genus == 0 and d <= 11:
                continue  # completed equal-involution genus-zero exclusions
            for m in range(1, d//2+1):
                cycles = 3*m-d+2-2*genus
                if 3 <= cycles <= d:
                    row = ('S' if m % 2 else 'A', d, m, genus)
                    if row not in rows:
                        rows.append(row)
    return rows


def s18_family_jobs():
    # n=6k, three 2^(3k) involutions and fourth 4^k 1^(2k): genus 1.
    return [('S' if k%2 else 'A',6*k,3*k,'.'.join(['4']*k)) for k in (3,5,4,2)]


def block_degree(block):
    return int(str(block).removeprefix('t').removeprefix('o').removeprefix('b').split(':')[0])


def mechanism_jobs(third_pairs=None):
    """Exact, explicitly bounded symmetry slices; never a degree census."""
    jobs=[dict(name='s18_position_spectrum',group='S',d=18,position_spectrum=True),
          dict(name='a25_unrestricted',group='A',d=25,rigid_blocks=[9,16],order=0)]
    small=[]
    def partitions(total,minimum=1):
        if total==0:
            yield ()
        for first in range(minimum,total+1):
            for rest in partitions(total-first,first):
                yield (first,)+rest
    for d in range(12,18):
        for parts in partitions(d):
            paths=[x for x in parts if x%2]
            cycles=[x for x in parts if not x%2]
            m=(d-len(paths))//2
            if not paths or len(parts)>16 or 2*m+(m if third_pairs is None else third_pairs)<d-1:
                continue
            small.append(dict(name='small_'+str(d)+'_'+'_'.join(map(str,parts)),
                              group='S' if m%2 else 'A',d=d,pair_blocks=cycles,pair_paths=paths))
    # Low-component profiles first; round-robin degrees prevents one degree
    # monopolizing the finite budget. Even-path pairs are a separate slice.
    small.sort(key=lambda j:(len(j['pair_blocks'])+len(j['pair_paths']),j['d'],j['name']))
    new=[]
    for blocks in rigid_glue_jobs(0):
        total=sum(blocks)+6
        if 26<=total<=34:
            new.append(dict(name='octahedral_o6_'+'_'.join(map(str,blocks)),group='A',d=total,
                            rigid_blocks=['o6',*blocks],order=0))
    new[:0]=[dict(name='octahedral_'+label,group='A',d=sum(map(block_degree,blocks)),
                  rigid_blocks=blocks,order=0) for label,blocks in
             [('o30_4',['o30',4]),('o30_3',['o30',3]),('o30',['o30'])]]
    # Existing direct sweep left these incomplete/unstarted. The completed
    # 2+4+4+4+4 slice is not duplicated. Old partial evidence is preserved.
    s18=[dict(name='s18_pair_'+'_'.join(map(str,p)),group='S',d=18,pair_blocks=p)
         for p in ([2,2,2,4,4,4],[2,2,2,2,2,4,4],[2,2,2,2,2,2,6],[2,2,2,2,2,2,2,4])]
    for i in range(max(len(small),len(new),len(s18))):
        for rows in (small,new,s18):
            if i<len(rows):jobs.append(rows[i])
    return jobs


def structure_jobs(third_pairs=None):
    """One-state algebra/lifts, balanced even paths, and obstruction audits."""
    jobs=[dict(name='one_state_symmetry',group='A',d=25,operation='braid-symmetry'),
          dict(name='independent_sympy',group='A',d=25,operation='independent'),
          dict(name='one_state_lifts',group='A',d=25,operation='research'),
          dict(name='one_state_rayleigh',group='A',d=25,operation='rayleigh')]
    # A colored pair graph is a multiset of even cycles, odd paths, and pairs
    # of even paths with opposite colors at each fixed length.
    kinds=sorted([(n,'cycle') for n in range(2,18,2)]+
                 [(n,'odd') for n in range(1,18,2)]+
                 [(2*n,'even') for n in range(2,9,2)])
    smaller=[]
    def extend(start,total,cycles,odd,even):
        fixed=len(odd)+2*len(even)
        m=(total-fixed)//2
        if 12<=total<=17 and even and 2*m+(m if third_pairs is None else third_pairs)>=total-1:
            label=lambda xs:'_'.join(map(str,xs)) or 'none'
            smaller.append(dict(name=f'even_{total}_c{label(cycles)}_p{label(odd)}_e{label(even)}',
                                group='S' if m%2 else 'A',d=total,pair_blocks=cycles,
                                pair_paths=odd,even_paths=even,resumable=True))
        for i in range(start,len(kinds)):
            weight,kind=kinds[i]
            if total+weight>17:continue
            extend(i,total+weight,cycles+([weight] if kind=='cycle' else []),
                   odd+([weight] if kind=='odd' else []),even+([weight//2] if kind=='even' else []))
    extend(0,0,[],[],[])
    smaller.sort(key=lambda j:(len(j['pair_blocks'])+len(j['pair_paths'])+2*len(j['even_paths']),j['d'],j['name']))
    audits=[dict(j,name='audit_'+j['name'],diagnose=True) for j in mechanism_jobs() if j['name'].startswith('octahedral')]
    jobs += smaller[:1]+audits[:1]
    jobs.append(dict(name='a17_remaining',group='A',d=17,pair_blocks=[2,2,2,2,2,2,4],pair_paths=[1],resumable=True))
    for i in range(1,max(len(smaller),len(audits))):
        for rows in (smaller,audits):
            if i<len(rows):jobs.append(rows[i])
    return jobs


def counterexample_pair_jobs():
    # New fixed-pair graphs: none repeats the completed S18 hexagon slice.
    profiles = ([6]*4, [4]*5, [8]*3, [6]*5, [4]*6, [10]*3,
                [4]*7, [8]*4, [4]*8,
                [6,6,8], [6,6,10], [6,8,8], [6,6,6,4], [6,8,10],
                [6,8,12], [6,6,6,6,2], [4,6,8,10], [6,10,12])
    return [dict(name='pair_'+'_'.join(map(str,p)),group='A' if sum(p)%4==0 else 'S',
                 d=sum(p),pair_blocks=p,resumable=True) for p in profiles]


def small_symmetric_slices(max_pairs=None):
    """Exact S9/S10 slices suggested by the S18 two-two and A25 three-one passports."""
    return structured_slices('S',(9,10),(0,1),max_pairs)


def structured_slices(group,degrees,genera,max_pairs=None,quotient=False):
    def partitions(total, minimum=1):
        if total == 0:
            yield ()
        for first in range(minimum, total+1):
            for rest in partitions(total-first, first):
                yield (first,)+rest

    rows = {}
    for d in degrees:
        types = []
        for partition in partitions(d):
            index = d-len(partition)
            if index and (group=='S' or index%2==0):
                text = '.'.join(map(str, (x for x in partition if x > 1)))
                types.append(dict(text=text, index=index, size=class_size(text, d)))
        by_text = {t['text']: t for t in types}

        def add(entries, mechanism, genus, always=False):
            key = (d, tuple(sorted(entries)))
            choices = []
            for fixed in set(entries):
                remaining = list(entries)
                remaining.remove(fixed)
                remaining.sort(key=lambda t: (by_text[t]['size'], t))
                total = by_text[remaining[0]]['size'] * by_text[remaining[1]]['size']
                # Approximate centralizer orbit cost; stabilizers can make the
                # actual count larger. The Rust scanner computes the exact count.
                score=max(1,(total*by_text[fixed]['size']+math.factorial(d)-1)//math.factorial(d)) if quotient else total
                choices.append((score,total, '/'.join([fixed, *remaining])))
            score, total, passport = min(choices)
            if max_pairs is not None and total > max_pairs and not always:
                return
            row = rows.setdefault(key, dict(group=group, d=d, passport=passport, total=total,
                                             genus=genus, mechanisms=set(), anchor=False))
            if quotient:row['estimated_pairs']=score
            row['mechanisms'].add(mechanism)
            row['anchor'] |= always

        # Nearest S18 analogue: two maximal odd involutions and two equal
        # complementary-index entries. These anchors run even if the broader
        # workload threshold is later reduced.
        m = d//2
        if m % 2 == 0:
            m -= 1
        involution = '.'.join(['2']*m)
        involution_index = m
        for genus in genera if group=='S' else ():
            other_index = d-1+genus-involution_index
            for other in types:
                if other['index'] == other_index:
                    add([involution]*2+[other['text']]*2, 's18-two-two', genus, True)

        # Full practical two-two and three-plus-exceptional strata. The latter
        # is the cycle-type shadow of the A25 mechanism. Only odd common types
        # can generate S_d in a three-one passport.
        for genus in genera:
            for i, left in enumerate(types):
                for right in types[i:]:
                    if (left['index']+right['index'] == d-1+genus
                            and (group=='A' or left['index'] % 2 or right['index'] % 2)):
                        add([left['text']]*2+[right['text']]*2, 's18-two-two', genus)
            for common in types:
                if group=='S' and common['index'] % 2 == 0:
                    continue
                exceptional_index = 2*d-2+2*genus-3*common['index']
                for exceptional in types:
                    if exceptional['index'] == exceptional_index:
                        add([common['text']]*3+[exceptional['text']], 'a25-three-one', genus)

    answer = []
    for row in rows.values():
        row['mechanisms'] = '+'.join(sorted(row['mechanisms']))
        row['tier'] = 0 if row.pop('anchor') else 1 if 'a25-three-one' in row['mechanisms'] else 2
        answer.append(row)
    return sorted(answer, key=lambda r: (r['tier'] != 0, r['d'], r.get('estimated_pairs',r['total']),r['total'], r['genus'], r['tier'], r['passport']))


def small_passport_plan(excluded=None,population=None):
    """All necessary-filter passports in the declared degree-nine population.

    Four arbitrary types are included, unlike the historical XXYY/XXXY queue.
    Class-size estimates select an ordered anchor and queue priority only.
    The Rust scanner measures the actual centralizer-reduced enumeration.
    """
    excluded=finite_passport_exclusions() if excluded is None else excluded
    def partitions(n,start=1):
        if not n:yield ()
        for first in range(start,n+1):
            for tail in partitions(n-first,first):yield (first,)+tail
    rows=[];omitted=[]
    for group,d,genera in (population or (('A',9,(1,2,3,4)),('S',9,(0,1)))):
        types=[]
        for p in partitions(d):
            index=d-len(p)
            if index and len(p)>=3 and (group=='S' or index%2==0):
                text='.'.join(str(x) for x in p if x>1)
                types.append((text,index,class_size(text,d)))
        for entries in itertools.combinations_with_replacement(types,4):
            index=sum(t[1] for t in entries)
            if index%2 or (group=='S' and not any(t[1]%2 for t in entries)):continue
            genus=1-d+index//2
            if genus not in genera:continue
            key=(group,d,tuple(sorted(t[0] for t in entries)))
            if key in excluded:
                omitted.append(dict(group=group,d=d,genus=genus,passport='/'.join(key[2]),evidence=excluded[key]))
                continue
            choices=[]
            for fixed in entries:
                remaining=list(entries);remaining.remove(fixed);remaining.sort(key=lambda t:(t[2],t[0]))
                raw=remaining[0][2]*remaining[1][2]
                estimate=max(1,(raw*fixed[2]+math.factorial(d)-1)//math.factorial(d))
                choices.append((estimate,raw,'/'.join([fixed[0],*[t[0] for t in remaining]])))
            estimate,raw,passport=min(choices)
            tag=hashlib.sha256('/'.join(key[2]).encode()).hexdigest()[:16]
            phase=(0 if (d==9 and genus==2) or (d==10 and group=='A' and genus<=1)
                   else 1 if d<=10 else 2) if population else (0 if group=='A' and genus==1 else 1 if group=='S' or genus==2 else 2)
            rows.append(dict(name=f'passport_v1_{group}{d}_{tag}',kind='passport',lane='exact',family='passport',
                group=group,d=d,genus=genus,passport=passport,phase=phase,
                estimated_pairs=estimate,raw_class_pairs=raw,cache_positive_witnesses=True))
    rows.sort(key=lambda r:(r['phase'],r['estimated_pairs'],r['group'],r['genus'],r['name']))
    return dict(rows=rows,excluded=omitted,
        scope=(f'All nonidentity necessary-filter passports in {population}; feasibility requires enumeration; no claim outside these strata' if population else 'All nonidentity necessary-filter A9 genus 1..4 and S9 genus 0..1 passports; feasibility requires enumeration; not other degrees or genera'),
        quotient='simultaneous S_d conjugation',pilot_seconds=15,chunk_seconds=120,
        estimates='class-size/centralizer heuristic, not actual candidate counts or an ETA')


def stability_passport_plan(excluded):
    return small_passport_plan(excluded,population=(('S',9,(2,3)),('A',10,(0,1,2,3)),
        ('S',10,(0,1)),('A',11,(0,1)),('S',11,(0,))))


def useful_passport_plan(excluded):
    return small_passport_plan(excluded,population=(('S',9,(2,3)),('A',10,(0,1,2,3)),
        ('S',10,tuple(range(6))),('A',11,(0,1)),('S',11,(0,))))


def covered_construction_targets(audit):
    # Only the declared parity-tight universe, only Factory's same nonidentity
    # min-three-cycle filter. This does not retire triples or other scanners.
    universes={('A',9):set(range(5)),('S',9):set(range(4)),('A',10):set(range(4))}
    covered=set()
    for target,genera in universes.items():
        rows=[r for r in audit['strata'] if (r['group'],r['degree'])==target]
        if {r['genus'] for r in rows}==genera and all(r['missing_passports']==0 for r in rows):covered.add(target)
    return covered


def small_degree_scope_audit(excluded):
    """Audit the entire parity-tight necessary-filter population, not a theorem.

    Input must be obtained through inherited_passport_exclusions; raw pinned
    intervals and historical evidence checks, not progress TSVs, authorize it.
    """
    population=(('A',9,tuple(range(5))),('S',9,tuple(range(4))),('A',10,tuple(range(4))))
    plan=small_passport_plan(excluded,population=population)
    strata=[]
    for group,d,genera in population:
        for genus in genera:
            missing=[r for r in plan['rows'] if (r['group'],r['d'],r['genus'])==(group,d,genus)]
            covered=[r for r in plan['excluded'] if (r['group'],r['d'],r['genus'])==(group,d,genus)]
            strata.append(dict(group=group,degree=d,genus=genus,total_passports=len(missing)+len(covered),
                excluded_with_provenance=len(covered),missing_passports=len(missing)))
    return dict(quotient='simultaneous S_d conjugation',strata=strata,
        necessary_filter_scope_complete=not plan['rows'],missing_rows=plan['rows'],
        exclusions=plan['excluded'],whole_degree_nonexistence_claim=False,
        proof_obligations=['review every-entry-at-least-three-cycles lemma',
            'review parity-tight genus bounds and identity-entry reduction',
            'independently review passport/scanner support and retained exclusions'],
        scope='Coverage audit of nonidentity min-three-cycle passports; not a new enumeration or blanket degree exclusion')


def inherited_passport_exclusions(parent):
    """Complete VERIFIED intervals only; any historical hard case prevents a skip.

    This is deliberately conservative even when a later replay resolved a lead.
    A progress table alone cannot authorize excluding a mathematical stratum.
    """
    # Invocation-local only: every root in discovery_chain is stopped. Repeated
    # ancestral prefixes share the strict validator, never progress-table claims.
    excluded=finite_passport_exclusions();exact_cache={}
    if not parent:return excluded
    for root,state,plan in discovery_chain(parent):
        for rows in plan['lanes'].values():
            for row in rows:
                if row['kind']!='passport':continue
                folder=root/'jobs'/row['name']
                coverage=discovery_coverage(row,folder,state['binary_sha256'],exact_cache=exact_cache)
                if not coverage['search_complete']:continue
                excluded[(row['group'],row['d'],tuple(sorted(row['passport'].split('/'))))]=str(folder)
    return excluded


class PassportQueue:
    """Cost-ranked pilots with resumable, bounded chunks and no dropped rows."""
    def __init__(self):self.heap=[];self.names=set()
    def push(self,row,attempts=0):
        if row['name'] in self.names:raise RuntimeError('Duplicate active passport queue row')
        priority=(row['phase'],row['estimated_pairs']*(1+attempts),attempts,row['name'])
        heapq.heappush(self.heap,(*priority,row));self.names.add(row['name'])
    def pop(self):
        if not self.heap:return None
        row=heapq.heappop(self.heap)[-1];self.names.remove(row['name']);return row


def structured_exclusions():
    """Validate the retained historical evidence before excluding any passport."""
    excluded={}
    # This full genus-zero A9 census is independently recorded in the project.
    progress=CAMPAIGNS/'smart_20260907_200313_392975/CENSUS_PROGRESS.tsv'
    with progress.open() as stream:
        prior=list(csv.DictReader(stream,delimiter='\t'))
    if len(prior)!=380 or not all(r['complete']=='True' for r in prior):
        raise RuntimeError('A9 exclusion requires the completed 380-slice census')
    for row in prior:
        excluded[('A',9,tuple(sorted(row['passport'].split('/'))))]=str(progress.relative_to(PROJECT))
    paths=list((RUNS/'hunts/20260925').glob('A10_three_involutions_*_exact/SUMMARY.txt'))
    paths+=list((RUNS/'hunts/20260925').glob('A11_three_involutions_*_exact/SUMMARY.txt'))
    prior11=CAMPAIGNS/'campaign_20260926_002343_265520'
    paths+=list(prior11.glob('jobs/*/attempt_*/SUMMARY.txt'))
    for path in paths:
        values=dict(l.split(' = ',1) for l in path.read_text().splitlines() if ' = ' in l)
        if (values.get('search mode')!='passport-exact' or values.get('scan range resolved')!='true'
                or int(values['scan start'])!=0 or values['scan next']!=values['scan total candidate pairs']
                or int(values.get('hard unsolved candidates',0))!=0):
            raise RuntimeError(f'Incomplete historical exclusion: {path}')
        if prior11 in path.parents and not (path.parent/'VERIFIED').exists():
            raise RuntimeError(f'Unverified historical exclusion: {path}')
        group,d=values['group'].split('_')
        excluded[(group,int(d),tuple(sorted(values['passport'].split('/'))))]=str(path.relative_to(PROJECT))
    if len(paths)!=26:raise RuntimeError('Expected nine A10 and seventeen A11/S11 completed exclusions')
    return excluded


def finite_passport_exclusions():
    """Only verified complete strata may be skipped by the finite S3 search."""
    excluded=structured_exclusions()
    record=json.loads((PROJECT/'evidence/structured_followup_20260929/ALTERNATING_COMPLETED.json').read_text())
    source=Path(relocated_path(PROJECT/record['source']))
    if len(record['passports'])!=140:raise RuntimeError('Expected 140 completed follow-up passports')
    for row in record['passports']:
        folder=source/'jobs'/row['job']
        coverage=exact_coverage(folder,record['binary_sha256'],row['passport'])
        if not coverage['search_complete']:raise RuntimeError(f'Unresolved finite exclusion: {folder}')
        excluded[(row['group'],row['degree'],tuple(sorted(row['passport'].split('/'))))]=str(folder)
    return excluded


def completed_s3_populations(root,parent=None):
    """Completed finite construction shards; saved closure leads remain replayable."""
    roots=[Path(root)]+[p for p,_,_ in discovery_chain(parent)]
    completed=set()
    for path in (p for r in roots for pattern in ('s3_exact_v1_*','fiber_exact_v1_*')
                 for p in r.glob(f'jobs/{pattern}/attempt_*/SUMMARY.txt')):
        if not (path.parent/'VERIFIED').exists():continue
        values=dict(l.split(' = ',1) for l in path.read_text().splitlines() if ' = ' in l)
        if values.get('family') in ('s3-exact','pair-fiber-exact') and values.get('finite construction complete')=='true':completed.add(path.parent.parent.name)
    return completed


def completed_fiber_validations(index=None):
    """Reuse curated, hash-checked finite harness evidence without rerunning it."""
    index=Path(index) if index is not None else PROJECT/'evidence/weak_refinement_20261004/FINITE_BRIDGES.json'
    if not index.exists():return set()
    completed=set()
    for row in json.loads(index.read_text())['shards']:
        summary=Path(row['summary'])
        if not summary.is_file() or hashlib.sha256(summary.read_bytes()).hexdigest()!=row['summary_sha256']:
            raise RuntimeError(f'Finite bridge evidence missing/changed: {summary}')
        values=dict(l.split(' = ',1) for l in summary.read_text().splitlines() if ' = ' in l)
        if values.get('family')!='pair-fiber-exact' or values.get('finite construction complete')!='true' or values.get('finite population resolved')!='true':
            raise RuntimeError(f'Finite bridge shard is not resolved: {summary}')
        source=Path(values['source certificate'])
        if not source.is_file() or hashlib.sha256(source.read_bytes()).hexdigest()!=row['source_sha256']:
            raise RuntimeError(f'Finite bridge source missing/changed: {source}')
        if row['name']!=summary.parent.name:raise RuntimeError('Finite bridge identity mismatch')
        completed.add(row['name'])
    return completed


def structured_followup(excluded=None):
    """Exact strata and symmetry slices; optional exclusions support isolated tests."""
    excluded=dict(structured_exclusions() if excluded is None else excluded)
    rows=[]
    for phase,group,d,genera in [('alternating','A',9,(1,)),('alternating','A',10,(0,1)),
                                 ('degree11','A',11,(0,1)),('degree11','S',11,(0,1))]:
        for row in structured_slices(group,(d,),genera,quotient=True):
            key=(group,d,tuple(sorted(row['passport'].split('/'))))
            if key in excluded:continue
            if any(d-sum(int(l)-1 for l in t.split('.'))<=2 for t in key[2]):
                # In a generating tuple some other entry joins the two cycles;
                # that pair is transitive. A one-cycle entry is already transitive.
                excluded[key]='theorem: an entry has at most two cycles; a transitive pair must exist'
                continue
            row=dict(row,phase=phase,centralizer=True)
            row['name']=f"structured_{group}{d}_"+hashlib.sha256('/'.join(key[2]).encode()).hexdigest()[:12]
            rows.append(row)
    mechanisms={}
    for third in (2,4,6,8):
        for row in mechanism_jobs(third)+structure_jobs(third):
            if row['group']!='A' or row['d'] not in (15,16,17) or row.get('pair_blocks') is None or third>row['d']//2:continue
            fixed=len(row.get('pair_paths',[]))+2*len(row.get('even_paths',[]))
            m=(row['d']-fixed)//2
            if third==m:continue  # all equal-involution profiles are already exhausted (C053/C055)
            key=(row['d'],third,tuple(row['pair_blocks']),tuple(row.get('pair_paths',[])),tuple(row.get('even_paths',[])))
            mechanisms[key]=dict(row,name=row['name']+f'_third{third}',third_pairs=third,resumable=True)
    forced_good=[]
    for key,row in list(mechanisms.items()):
        m=(row['d']-len(row.get('pair_paths',[]))-2*len(row.get('even_paths',[])))//2
        components=len(row['pair_blocks'])+len(row.get('pair_paths',[]))+2*len(row.get('even_paths',[]))
        reason=None
        if 2*m+row['third_pairs']<row['d']+1:
            reason='nonnegative genus forces the fourth entry to have at most two cycles'
        elif components==1:reason='the fixed first pair is already transitive'
        elif row['third_pairs']<components-1:reason='too few third-entry transpositions to join all fixed pair components'
        if reason:
            forced_good.append(dict(row,reason=reason))
            del mechanisms[key]
    return dict(exact=rows,mechanisms=sorted(mechanisms.values(),key=lambda r:(len(r['pair_blocks'])+len(r.get('pair_paths',[]))+2*len(r.get('even_paths',[])),r['d'],r['name'])),
                mechanism_excluded=forced_good,mechanism_prior_completed='evidence/one_state_structure_20260927/COMPLETED_SUMMARY.txt (C053/C055): equal third type excluded',
                excluded=[dict(group=g,d=d,passport='/'.join(p),evidence=e) for (g,d,p),e in sorted(excluded.items())])


def exact_scan_enumerated(values, scan):
    if int(values['scan next']) >= scan[0]+scan[1]:
        return True
    return (scan[0]==0 and scan[1]>=int(values.get('scan raw class pairs',scan[1]))
            and int(values['scan next'])>=int(values['scan total candidate pairs']))


def exact_scan_complete(values, scan):
    return values.get('scan range resolved')=='true' and exact_scan_enumerated(values,scan)


def empty_exact_coverage():
    return dict(next=0,total=None,enumeration_complete=False,traversal_exhausted=False,
                range_resolved=True,search_complete=False)


def extend_exact_coverage(previous,meta,values,binary_sha,passport):
    """The same strict interval validation for raw replay and RAM snapshots."""
    if meta!=dict(binary_sha256=binary_sha,passport=passport,enumeration='centralizer-v1'):
        raise RuntimeError('Refusing exact cursor reuse with changed binary, passport or enumeration')
    cursor=previous['next'];end=int(values['scan next']);total=int(values['scan total candidate pairs'])
    if (previous['enumeration_complete'] or int(values['scan start'])!=cursor or not cursor<=end<=total
            or (previous['total'] is not None and previous['total']!=total)
            or values.get('scan enumeration')!='centralizer-v1'):
        raise RuntimeError('Exact enumeration coverage has a gap, overlap or changed candidate count')
    exhausted=end==total
    resolved=previous['range_resolved'] and int(values.get('hard unsolved candidates',0))==0
    return dict(next=end,total=total,enumeration_complete=exhausted,traversal_exhausted=exhausted,
                range_resolved=resolved,search_complete=exhausted and resolved)


def exact_coverage(job_dir,binary_sha,passport,prefix=None,memo=None):
    prefixes=coverage_prefix_list(prefix)
    key=(str(job_dir),binary_sha,passport,tuple(map(str,prefixes)))
    if memo is not None and key in memo:return dict(memo[key])
    coverage=exact_coverage(Path(prefixes[-1]),binary_sha,passport,prefixes[:-1],memo) if prefixes else empty_exact_coverage()
    for path in sorted(job_dir.glob('attempt_*/SUMMARY.txt'),key=lambda p:int(p.parent.name.split('_')[-1])):
        if not (path.parent/'VERIFIED').exists():continue
        meta=json.loads((path.parent/'SCAN_INPUT.json').read_text())
        v=dict(l.split(' = ',1) for l in path.read_text().splitlines() if ' = ' in l)
        coverage=extend_exact_coverage(coverage,meta,v,binary_sha,passport)
    if memo is not None:memo[key]=dict(coverage)
    return coverage


def pair_coverage(job_dir, binary_sha, profile=None, prefix=None):
    prefixes=coverage_prefix_list(prefix)
    previous=pair_coverage(Path(prefixes[-1]),binary_sha,profile,prefixes[:-1]) if prefixes else None
    cursor=(previous['next_h'],previous['next_third']) if previous else (0,0)
    resolved=previous['range_resolved'] if previous else True
    exhausted=previous['traversal_exhausted'] if previous else False
    attempts=previous['verified_attempts'] if previous else 0
    for path in sorted(job_dir.glob('attempt_*/SUMMARY.txt'),key=lambda p:int(p.parent.name.split('_')[-1])):
        if not (path.parent/'VERIFIED').exists():continue
        meta=json.loads((path.parent/'PAIR_INPUT.json').read_text())
        if meta['binary_sha256']!=binary_sha:
            raise RuntimeError('Refusing cursor reuse with a changed enumeration binary')
        recorded={k:meta.get(k) for k in ('cycles','paths','even_paths')}
        for key in ('third_pairs','third_order','target_group','quotient_h'):
            if key in meta:recorded[key]=meta[key]
        if profile is None:profile=recorded
        if profile!=recorded:
            raise RuntimeError('Refusing cursor reuse with a changed pair profile')
        values=dict(l.split(' = ',1) for l in path.read_text().splitlines() if ' = ' in l)
        start=(int(values['scan start h']),int(values['scan start third']))
        end=(int(values['scan next h']),int(values['scan next third']))
        if exhausted or start!=cursor or end<start:
            raise RuntimeError('Pair enumeration coverage has a gap, overlap or invalid cursor')
        cursor=end
        resolved &= values['scan range resolved']=='true'
        exhausted=values['traversal exhausted']=='true'
        attempts+=1
    return dict(next_h=cursor[0],next_third=cursor[1],range_resolved=resolved,
                traversal_exhausted=exhausted,search_complete=exhausted and resolved,
                verified_attempts=attempts,binary_sha256=binary_sha)


def completed_pair_profiles():
    """Validate previous prefixes, rather than treating a DONE file as proof."""
    profiles={}
    campaigns=(('campaign_20260926_235453_768706','MECHANISM_JOBS.json'),
               ('campaign_20260927_115848_035745','STRUCTURE_JOBS.json'),
               ('campaign_20260929_160438_677235','STRUCTURED_FOLLOWUP.json'))
    for name,filename in campaigns:
        root=CAMPAIGNS/name
        if not (root/'inputs'/filename).exists():continue
        sha=json.loads((root/'STATE.json').read_text())['binary_sha256']
        rows=json.loads((root/'inputs'/filename).read_text())
        if isinstance(rows,dict):rows=rows['mechanisms']
        pending=candidates_from_files(list(root.glob('jobs/*/attempt_*/hard_candidates.txt'))+list(root.glob('jobs/*/attempt_*/survivors.txt')),weak_only=True)
        for row in rows:
            if row.get('pair_blocks') is None:continue
            folder=root/'jobs'/row['name']
            if not pair_coverage(folder,sha)['traversal_exhausted']:continue
            if pending:
                local=candidates_from_files(list(folder.glob('attempt_*/hard_candidates.txt'))+list(folder.glob('attempt_*/survivors.txt')),weak_only=True)
                if any(seeds & pending.get(group,set()) for group,seeds in local.items()):continue
            # These campaign-level resolved audits are pinned in C053/C055/C069.
            # An unexhausted profile is NEVER imported as an exclusion.
            d=row['d'];m=(d-len(row.get('pair_paths',[]))-2*len(row.get('even_paths',[])))//2
            key=(d,tuple(row['pair_blocks']),tuple(row.get('pair_paths',[])),tuple(row.get('even_paths',[])),row.get('third_pairs',m))
            profiles[key]=str(folder.relative_to(PROJECT))
    return profiles


def max_power_index(degree,order):
    """Exact maximum index of a permutation satisfying c^order=1."""
    dp=[0]*(degree+1)
    for n in range(1,degree+1):
        dp[n]=max(dp[n-l]+l-1 for l in range(1,min(n,order)+1) if order%l==0)
    return dp[degree]


def discovery_plan(completed=None,frontier=False):
    """Frozen disjoint construction strata, with a separate non-census lane."""
    completed=completed_pair_profiles() if completed is None else completed
    lanes={'middle':[],'powers':[],'rigid':[]}
    if frontier:lanes['small']=[]
    excluded=[]
    kinds=sorted([(n,'cycle') for n in range(2,18,2)]+[(n,'odd') for n in range(1,18,2)]+[(2*n,'even') for n in range(2,9,2)])
    def extend(start,d,cycles,paths,even):
        components=len(cycles)+len(paths)+2*len(even)
        m=(d-len(paths)-2*len(even))//2
        if (9 if frontier else 12)<=d<=17 and 1<components<=16 and m>0:
            identical=not even and all(n==2 for n in cycles) and all(n==1 for n in paths)
            if not identical:
                base=dict(d=d,pair_blocks=cycles,pair_paths=paths,even_paths=even)
                for k in range(1,d//2+1):
                    key=(d,tuple(cycles),tuple(paths),tuple(even),k)
                    if k<components-1 or 2*m+k<d+1:continue
                    if key in completed:
                        excluded.append(dict(profile=list(key),evidence=completed[key]));continue
                    if d==12 and not paths and not even and k==m:
                        excluded.append(dict(profile=list(key),evidence='C045: complete A12 fixed-point-free involution symmetry slice'));continue
                    group='A' if m%2==0 and k%2==0 else 'S'
                    lane='small' if d<12 else 'middle'
                    lanes[lane].append(dict(base,lane=lane,kind='pair',group=group,third_pairs=k,third_order=2))
                for q in ((3,4,5,6) if frontier and d<12 else (3,4)):
                    max_index=min(d-3,max_power_index(d,q))
                    if max_index<components-1 or 2*m+max_index<d+1:continue
                    groups=[('A' if m%2==0 else 'S')] if q%2 else (['S','A'] if d<12 else ['auto'])
                    for group in groups:
                        if group=='A' and m%2:continue
                        lane='small' if d<12 else 'powers'
                        lanes[lane].append(dict(base,lane=lane,kind='pair',group=group,third_order=q))
        for i in range(start,len(kinds)):
            weight,kind=kinds[i]
            if d+weight>17:continue
            extend(i,d+weight,cycles+([weight] if kind=='cycle' else []),
                   paths+([weight] if kind=='odd' else []),even+([weight//2] if kind=='even' else []))
    extend(0,0,[],[],[])
    # Base modular representations have arbitrary transitive dessins, unlike
    # the previously exhausted tetrahedral/octahedral/source block library.
    for modulus in range(2,8):
        for base in range(3,17):
            n=modulus*base
            if n>34 or (modulus%2==0 and base%2) or (modulus%3==0 and base%3):continue
            for seed in range(16 if frontier else 8):
                label=f'b{n}:{modulus}:{seed}'
                for blocks in ([label],*[ [label,old] for old in (9,16) if n+old<=34]):
                    for order in ((3,4,5,6,0) if frontier else (3,4,6,0)):
                        lanes['rigid'].append(dict(lane='rigid',kind='rigid',d=sum(map(block_degree,blocks)),
                                                   group='auto',rigid_blocks=blocks,order=order))
    for lane,rows in lanes.items():
        for row in rows:
            tag=hashlib.sha256(json.dumps(row,sort_keys=True).encode()).hexdigest()[:16]
            row['name']=f"discover_{lane}_{row['d']}_{tag}"
        if lane!='rigid':
            rows.sort(key=lambda r:((r['group']!='S',r['d'],r['third_order']==2) if lane=='small' else (),
                                   len(r['pair_blocks'])+len(r['pair_paths'])+2*len(r['even_paths']),r['d'],r.get('third_order'),r['name']))
        else:rows.sort(key=lambda r:(len(r['rigid_blocks']),r['d'],r['order']==0,r['name']))
    return dict(lanes=lanes,excluded=excluded,scope='Selected exact symmetry slices and biased diversity sampling; not whole-degree censuses',
                s11='C069: previous capped candidates resolved; replay only new surviving candidates',
                allocation=['small','powers','small','rigid','middle','diversity'] if frontier else ['powers','middle','powers','rigid','powers','diversity','middle','powers'],
                policy='frontier-v1' if frontier else 'discovery-v1')


def coverage_prefix_list(prefix):
    rows=list(prefix) if isinstance(prefix,list) else [prefix] if prefix else []
    rows=[relocated_path(p) for p in rows]
    if len({str(Path(p).resolve()) for p in rows})!=len(rows):raise RuntimeError('Duplicate/cyclic enumeration prefix')
    return rows


def discovery_prefixes(row):
    return coverage_prefix_list(row.get('coverage_prefixes',row.get('coverage_prefix')))


def rigid_coverage(folder,sha,profile,prefix=None):
    prefixes=coverage_prefix_list(prefix)
    previous=rigid_coverage(Path(prefixes[-1]),sha,profile,prefixes[:-1]) if prefixes else None
    cursor=previous['next_branch'] if previous else 0
    total=previous['total'] if previous else None
    exhausted=previous['traversal_exhausted'] if previous else False
    for path in sorted(folder.glob('attempt_*/SUMMARY.txt'),key=lambda p:int(p.parent.name.split('_')[-1])):
        if not (path.parent/'VERIFIED').exists():continue
        if json.loads((path.parent/'RIGID_INPUT.json').read_text())!=dict(binary_sha256=sha,**profile):
            raise RuntimeError('Changed rigid representation or binary cannot reuse its branch cursor')
        v=dict(l.split(' = ',1) for l in path.read_text().splitlines() if ' = ' in l)
        if 'structural exclusion' in v:
            exhausted=True;break
        if exhausted or int(v['rigid start branch'])!=cursor:raise RuntimeError('Rigid branch coverage gap/overlap')
        n=int(v['rigid total branches']);next_cursor=int(v['rigid next branch'])
        if not cursor<=next_cursor<=n or (total is not None and n!=total):raise RuntimeError('Invalid rigid branch coverage')
        cursor=next_cursor;total=n;exhausted=v['rigid traversal exhausted']=='true'
    return dict(next_branch=cursor,total=total,traversal_exhausted=exhausted)


def discovery_coverage(row,folder,digest,exact_cache=None):
    sha=row.get('enumeration_sha256',digest);prefix=discovery_prefixes(row)
    if 'prior_attempts' in row and sum(len(list(Path(p).glob('attempt_*/VERIFIED'))) for p in prefix)!=row['prior_attempts']:
        raise RuntimeError('Inherited verified attempts are missing or changed; refusing cursor reset')
    if str(folder.resolve()) in {str(Path(p).resolve()) for p in prefix}:raise RuntimeError('Current job cannot be its own prefix')
    if row['kind']=='passport':return exact_coverage(folder,sha,row['passport'],prefix,memo=exact_cache)
    if row['kind']=='pair':return pair_coverage(folder,sha,discovery_pair_profile(row),prefix)
    return rigid_coverage(folder,sha,dict(blocks=row['rigid_blocks'],order=row['order'],enumeration='root-branches-v1'),prefix)


def discovery_chain(parent):
    seen=set()
    while parent:
        parent=Path(relocated_path(parent));identity=str(parent.resolve())
        if identity in seen:raise ValueError('Cycle in discovery predecessors')
        seen.add(identity)
        state=json.loads((parent/'STATE.json').read_text())
        if state['status'] not in ('stopped','complete','failed'):raise ValueError('Stop predecessor before inheriting it')
        if state.get('plan') not in ('discovery','discovery-frontier'):raise ValueError('Predecessor must be a discovery campaign')
        plan=json.loads((parent/'inputs/DISCOVERY_PLAN.json').read_text())
        yield parent,state,plan
        parent=plan.get('predecessor')


def inherit_discovery(schedule,parent):
    """Completed profiles stay excluded; partial cursors keep their OLD executable."""
    chain=list(discovery_chain(parent))
    for root,state,old in reversed(chain):
        prior={r['name']:r for rows in old['lanes'].values() for r in rows}
        for lane,rows in schedule['lanes'].items():
            kept=[]
            for row in rows:
                old_row=prior.get(row['name'])
                if old_row:
                    folder=root/'jobs'/row['name'];coverage=discovery_coverage(old_row,folder,state['binary_sha256'])
                    if coverage['traversal_exhausted']:
                        schedule['excluded'].append(dict(profile=row['name'],reason='verified predecessor enumeration; unresolved leads imported separately',evidence=str(folder)))
                        continue
                    # Retain original binary and the complete ordered prefix.
                    prefixes=discovery_prefixes(old_row)+[str(folder)]
                    coverage_prefix_list(prefixes)
                    sha=old_row.get('enumeration_sha256',state['binary_sha256'])
                    binary=relocated_path(old_row.get('enumeration_binary',str(root/'bin'/sha)))
                    row.update(coverage_prefixes=prefixes,enumeration_sha256=sha,enumeration_binary=str(binary),
                               prior_attempts=old_row.get('prior_attempts',0)+len(list(folder.glob('attempt_*/VERIFIED'))))
                kept.append(row)
            schedule['lanes'][lane]=kept
    schedule['predecessor']=str(chain[0][0])
    return schedule


def candidate_paths(root,verified=True):
    """Read closed attempts by default; live worker files are not snapshots."""
    root=Path(relocated_path(root))
    paths=[p for p in itertools.chain(*(root.glob('jobs/*/attempt_*/'+name) for name in ('hard_candidates.txt','survivors.txt','construction_survivors.txt')))
           if p.stat().st_size and (not verified or (p.parent/'VERIFIED').exists())]
    for name in ('INHERITED_PENDING.txt','RECOVERED_PENDING.txt'):
        inherited=root/'inputs'/name
        if inherited.exists() and inherited.stat().st_size:paths.append(inherited)
    return paths


class VerifiedDiscoveryIndex:
    """In-memory view of immutable verified attempts; no active-worker reads.

    Reconstruct once on resume, then register each attempt exactly once. This
    avoids repeated directory walks on T7 while keeping the raw audit files.
    """
    def __init__(self,root):
        self.root=Path(root);self.attempts=set();self.paths=set();self.metrics={};self.totals=[0,0,0]
        self.records={};self.certificate_paths=set();self.scan_inputs={};self.job_summaries={}
        # Exact strings, never probabilistic hashes. A bounded positive cache
        # may forget an old witness and do extra replay, never certify absence.
        self._resolved=OrderedDict();self._resolved_limit=2_000_000
        self._pending={};self._by_folder={};self._prefix_candidates={};self.candidate_files_read=0
        jobs=self.root/'jobs'
        if jobs.exists():
            # Directory listings amortize ExFAT metadata work; avoid many
            # separate exists/stat/glob calls for each immutable attempt.
            with os.scandir(jobs) as job_entries:
                for job in job_entries:
                    if job.name.startswith('._') or not job.is_dir():continue
                    with os.scandir(job.path) as attempts:
                        for attempt in attempts:
                            if attempt.name.startswith('attempt_') and attempt.is_dir():
                                self.register(Path(attempt.path))
        for name in ('INHERITED_PENDING.txt','RECOVERED_PENDING.txt'):
            path=self.root/'inputs'/name
            if path.exists() and path.stat().st_size:
                self.paths.add(path);self._register_candidates([path])

    def register(self,attempt,row=None,values=None):
        attempt=Path(attempt)
        if attempt in self.attempts:return
        with os.scandir(attempt) as entries:files={e.name:e for e in entries}
        if 'VERIFIED' not in files or 'SUMMARY.txt' not in files:return
        if row is None:row=json.loads((attempt/'JOB_INPUT.json').read_text())
        if values is None:values=dict(l.split(' = ',1) for l in (attempt/'SUMMARY.txt').read_text().splitlines() if ' = ' in l)
        scan_input=json.loads((attempt/'SCAN_INPUT.json').read_text()) if row.get('kind')=='passport' else None
        paths=[]
        for name in ('hard_candidates.txt','survivors.txt','construction_survivors.txt'):
            path=attempt/name
            if name in files and files[name].stat().st_size:paths.append(path)
        # Parse the entire immutable attempt before committing metrics. Replays
        # can resolve earlier leads even when their candidate file is empty.
        self._register_candidates(paths,attempt/'RESULTS.tsv',force=row.get('kind')=='replay',results_exists='RESULTS.tsv' in files)
        metric=self.metrics.setdefault(row['lane'],[0,0,0,0.0]);metric[0]+=1
        generated=int(values.get('accepted generating seeds',0));hard=int(values.get('hard unsolved candidates',0))
        metric[1]+=generated;metric[2]+=hard;metric[3]+=float(values.get('elapsed seconds',0))
        for i,n in enumerate([generated,int(values.get('solved seeds',0)),hard]):self.totals[i]+=n
        self.paths.update(paths)
        self.records[attempt/'SUMMARY.txt']=dict(values)
        if scan_input is not None:self.scan_inputs[attempt]=scan_input
        self.job_summaries.setdefault(attempt.parent,[]).append(attempt/'SUMMARY.txt')
        self.certificate_paths.update(attempt/name for name in files if name.startswith('CERTIFIED_') and name.endswith('.txt'))
        self.attempts.add(attempt)

    def _register_candidates(self,paths,results=None,force=False,results_exists=None):
        resolved=set();local=set();parsed=[]
        if (paths or force) and results is not None and (results_exists if results_exists is not None else results.exists()):
            with results.open() as stream:
                for row in csv.DictReader(stream,delimiter='\t'):
                    status=row.get('status','');seed=row.get('seed')
                    if status.endswith('_found') or status.startswith('certified_'):local.add(seed)
                    if status in ('strong_found','weak_found') or status.startswith('certified_'):resolved.add(seed)
        for path in paths:
            for line in path.read_text().splitlines():
                if not line.strip():continue
                if not line.startswith('group=') or ' seed=' not in line:
                    raise ValueError(f'Malformed candidate in {path}: {line[:100]}')
                gd=line.split()[0].removeprefix('group=');seed=line.split(' seed=',1)[1]
                if seed not in resolved and not (path.name=='survivors.txt' and seed in local):
                    parsed.append((path.parent.parent,gd,seed))
        for seeds in self._pending.values():seeds.difference_update(resolved)
        for seed in resolved:
            self._resolved[seed]=None;self._resolved.move_to_end(seed)
        for folder,gd,seed in parsed:
            if seed in self._resolved:continue
            self._pending.setdefault(gd,set()).add(seed)
            self._by_folder.setdefault(folder,{}).setdefault(gd,set()).add(seed)
        while len(self._resolved)>self._resolved_limit:self._resolved.popitem(last=False)
        self.candidate_files_read+=len(paths)

    def pending(self):
        """An owned snapshot; no filesystem access, no active-attempt reads."""
        return {gd:set(seeds) for gd,seeds in self._pending.items() if seeds}

    def local_pending(self,folders,pending):
        if not pending:return 0
        groups={}
        for folder in folders:
            folder=Path(folder)
            local=self._by_folder.get(folder,{})
            if not folder.is_relative_to(self.root):
                # Inherited pinned prefixes are immutable, but their local
                # lead association is not contained in the flat recovery input.
                # Load it once, so a pending ancestral cap cannot disappear
                # from a row's unresolved count merely because it is old.
                if folder not in self._prefix_candidates:
                    paths=[p for pattern in ('attempt_*/hard_candidates.txt','attempt_*/survivors.txt','attempt_*/construction_survivors.txt')
                           for p in folder.glob(pattern) if (p.parent/'VERIFIED').exists()]
                    self._prefix_candidates[folder]=candidates_from_files(paths,weak_only=True)
                local=self._prefix_candidates[folder]
            for gd,seeds in local.items():
                groups.setdefault(gd,set()).update(seeds & pending.get(gd,set()))
        return sum(map(len,groups.values()))

    def live_report(self):
        known={certificate_identity(p) for p in (self.root/'known').glob('CERTIFIED_*.txt')}
        catalog=known|{certificate_identity(p) for p in (self.root/'resolved').glob('CERTIFIED_*.txt')}
        values=dict(updated=dt.datetime.now().isoformat(),verified_attempts=len(self.attempts),
            accepted_generating_occurrences=self.totals[0],solved_occurrences=self.totals[1],unresolved_historical_events=self.totals[2],
            certified_components=len(catalog),new_components=len(catalog-known),units='occurrences including repeats; not distinct tuples/orbits',
            cached_candidate_files=self.candidate_files_read,cached_positive_seeds=len(self._resolved),
            source='immutable VERIFIED attempts; RAM-indexed pending leads and metrics; full tables at round/stop boundaries')
        atomic_json(self.root/'LIVE_REPORT.json',values)
        pending=self.pending()
        atomic_json(self.root/'DISCOVERY_PENDING.json',{group:len(seeds) for group,seeds in pending.items()})
        with (self.root/'LANE_METRICS.tsv').open('w') as stream:
            writer=csv.writer(stream,delimiter='\t')
            writer.writerow(['lane','verified_attempts','accepted_seed_occurrences','unresolved_events','elapsed_worker_seconds'])
            writer.writerows([lane,*metric] for lane,metric in sorted(self.metrics.items()))
        (self.root/'REPORT.txt').write_text('Hurwitz live campaign report\n'+''.join(f'{k}: {v}\n' for k,v in values.items())+
            'Detailed DISCOVERY_PROGRESS.tsv is rebuilt at round/stop boundaries.\n'
            'No whole-degree census or time-to-discovery claim. Cap/timeout events are not counterexamples.\n')


class VerifiedCoverageCache:
    """Session-local exact coverage, validated once then invalidated on VERIFIED.

    Raw pinned inputs remain the authority. This is deliberately NOT persisted:
    every resume revalidates history, so stale/mutated caches cannot reset a cursor.
    """
    def __init__(self,index=None):
        self.index=index;self._values={};self.hits=0;self.misses=0

    def get(self,row,folder,digest):
        folder=Path(folder)
        contract=(digest,json.dumps(row,sort_keys=True))
        saved=self._values.get(folder)
        if saved is None or saved[0]!=contract:
            self.misses+=1
            if (self.index is not None and row['kind']=='passport' and not discovery_prefixes(row)
                    and not row.get('prior_attempts',0)):
                coverage=empty_exact_coverage()
                for path in sorted(self.index.job_summaries.get(folder,[]),key=lambda p:int(p.parent.name.split('_')[-1])):
                    coverage=extend_exact_coverage(coverage,self.index.scan_inputs[path.parent],self.index.records[path],
                                                  row.get('enumeration_sha256',digest),row['passport'])
            else:coverage=discovery_coverage(row,folder,digest)
            saved=(contract,coverage)
            self._values[folder]=saved
        else:self.hits+=1
        return dict(saved[1])

    def invalidate(self,folder):
        self._values.pop(Path(folder),None)


@contextmanager
def atomic_text(path):
    temporary=path.with_name(path.name+f'.{os.getpid()}.tmp')
    try:
        with temporary.open('w') as stream:yield stream
        temporary.replace(path)
    finally:
        if temporary.exists():temporary.unlink()


def write_discovery_progress(root,schedule,digest,coverage_cache,index,exact_mission=False):
    """One coverage lookup per row, using immutable session caches when present."""
    pending=index.pending() if index is not None else candidates_from_files(candidate_paths(root),weak_only=True)
    with ExitStack() as stack:
        writer=csv.writer(stack.enter_context(atomic_text(root/'DISCOVERY_PROGRESS.tsv')),delimiter='\t')
        writer.writerow(['job','lane','degree','enumeration_complete','unresolved_seeds','search_complete'])
        exact=None
        if exact_mission:
            exact=csv.writer(stack.enter_context(atomic_text(root/'EXACT_PASSPORT_PROGRESS.tsv')),delimiter='\t')
            exact.writerow(['group','degree','genus','passport','next','total','enumeration_complete','unresolved_seeds','search_complete'])
        for row in itertools.chain.from_iterable(schedule['lanes'].values()):
            folder=root/'jobs'/row['name'];c=coverage_cache.get(row,folder,digest)
            folders=[folder,*map(Path,discovery_prefixes(row))]
            if index is not None:n=index.local_pending(folders,pending)
            else:
                paths=[p for f in folders for pattern in ('attempt_*/hard_candidates.txt','attempt_*/survivors.txt','attempt_*/construction_survivors.txt')
                       for p in f.glob(pattern) if (p.parent/'VERIFIED').exists()]
                local=candidates_from_files(paths,weak_only=True) if pending else {}
                n=sum(len(seeds & pending.get(gd,set())) for gd,seeds in local.items())
            writer.writerow([row['name'],row['lane'],row['d'],c['traversal_exhausted'],n,c['traversal_exhausted'] and n==0])
            if exact is not None and row['lane']=='exact':
                exact.writerow([row['group'],row['d'],row['genus'],row['passport'],c['next'],c['total'],c['enumeration_complete'],n,c['enumeration_complete'] and n==0])
    atomic_json(root/'DISCOVERY_PENDING.json',{gd:len(seeds) for gd,seeds in pending.items()})
    return pending


class DiscoveryProofCache:
    """VERIFIED-only exact words and completed representation/order records.

    Children receive immutable hashed snapshots and replay every imported word.
    Misses, caps and partial representations are never negative certificates.
    """
    def __init__(self,root,predecessor=None):
        self.root=Path(root);self.words={};self.representations=set();self.limit=50000;self.snapshots={}
        self.directory=self.root/'proof_cache';self.directory.mkdir(exist_ok=True)
        ancestors=[self.root]+[p for p,_,_ in discovery_chain(predecessor)] if predecessor else [self.root]
        for parent in reversed(ancestors):
            for path in (parent/'proof_cache').glob('WORDS_*.tsv'):self._words(path)
            ledger=parent/'proof_cache/REPRESENTATIONS.tsv'
            if ledger.exists():
                with ledger.open() as stream:
                    for row in csv.DictReader(stream,delimiter='\t'):self.representations.add((row['representation_key'],int(row['entry_order'])))
            plan=parent/'inputs/DISCOVERY_PLAN.json'
            if plan.exists():
                for row in json.loads(plan.read_text()).get('lanes',{}).get('rigid',[]):
                    for summary in (parent/'jobs'/row['name']).glob('attempt_*/SUMMARY.txt'):
                        if (summary.parent/'VERIFIED').exists():self._representation(summary)
        for path in self.root.glob('jobs/*/attempt_*/POSITIVE_WITNESSES.tsv'):
            if (path.parent/'VERIFIED').exists():self._words(path)
        self._save_representations()

    def _words(self,path):
        changed=set()
        with Path(path).open() as stream:
            reader=csv.DictReader(stream,delimiter='\t')
            if reader.fieldnames!=['degree','moves','seed']:raise RuntimeError('Malformed positive-word cache')
            for row in reader:
                d=int(row['degree']);moves=row['moves'].split(',') if row['moves'] else []
                if not 4<=d<=34 or len(moves)>4096 or any(int(m) not in range(6) for m in moves) or not row['seed'].startswith('['):raise RuntimeError('Invalid positive-word row')
                values=self.words.setdefault(d,{})
                if len(values)<self.limit or row['seed'] in values:
                    old=values.get(row['seed'])
                    if old is None or len(moves)<(len(old.split(',')) if old else 0):
                        values[row['seed']]=row['moves'];changed.add(d);self.snapshots.pop(d,None)
        return changed

    def _representation(self,path):
        before=len(self.representations)
        values=dict(l.split(' = ',1) for l in Path(path).read_text().splitlines() if ' = ' in l)
        if values.get('search complete')=='true' and values.get('hard unsolved candidates','0')=='0' and values.get('representation key'):
            self.representations.add((values['representation key'],int(values['entry order divides'])))
        return len(self.representations)!=before

    def _save_representations(self):
        content='representation_key\tentry_order\n'+''.join(f'{key}\t{order}\n' for key,order in sorted(self.representations))
        target=self.directory/'REPRESENTATIONS.tsv';temp=target.with_suffix('.tmp');temp.write_text(content);temp.replace(target)

    def register(self,attempt):
        attempt=Path(attempt)
        if not (attempt/'VERIFIED').exists():return
        path=attempt/'POSITIVE_WITNESSES.tsv'
        if path.exists():
            for d in self._words(path):
                target=self.directory/f'WORDS_{d}.tsv';temp=target.with_suffix('.tmp');temp.write_text(self._content(d));temp.replace(target)
        if (attempt/'SUMMARY.txt').exists() and self._representation(attempt/'SUMMARY.txt'):self._save_representations()

    def _content(self,d):
        return 'degree\tmoves\tseed\n'+''.join(f'{d}\t{word}\t{seed}\n' for seed,word in sorted(self.words.get(d,{}).items()))

    def snapshot(self,d):
        if d in self.snapshots:return self.snapshots[d]
        content=self._content(d);tag=hashlib.sha256(content.encode()).hexdigest()[:20]
        path=self.root/'inputs'/f'POSITIVE_{d}_{tag}.tsv'
        if not path.exists():path.write_text(content)
        self.snapshots[d]=str(path);return str(path)

    def representation_snapshot(self):
        content=(self.directory/'REPRESENTATIONS.tsv').read_bytes();tag=hashlib.sha256(content).hexdigest()[:20]
        path=self.root/'inputs'/f'REPRESENTATIONS_{tag}.tsv'
        if not path.exists():path.write_bytes(content)
        return str(path)


def recover_interrupted_candidates(root):
    """Snapshot stopped, unverified attempts before this worker starts children.

    Original outputs remain untouched. Malformed quiescent files still fail
    validation: ignoring a partial seed would conceal an unresolved lead.
    """
    paths=[p for p in candidate_paths(root,verified=False)
           if not (p.parent/'VERIFIED').exists()]
    pending=candidates_from_files(paths,weak_only=True)
    destination=root/'inputs/RECOVERED_PENDING.txt'
    temporary=destination.with_suffix('.tmp')
    temporary.write_text(''.join(f'group={group} seed={seed}\n'
                                for group,seeds in sorted(pending.items())
                                for seed in sorted(seeds)))
    temporary.replace(destination)
    return pending


def exploration_passport(degree,group,genus,seed):
    """Biased genus-balanced sampling, including four different cycle types."""
    if not 0<=genus<=degree-5:raise ValueError('Weak-counterexample genus must be in 0..d-5')
    def partitions(n,minimum=1):
        if not n:yield ()
        for k in range(minimum,n+1):
            for tail in partitions(n-k,k):yield (k,)+tail
    types=[p for p in partitions(degree) if 3<=len(p)<degree and (group=='S' or (degree-len(p))%2==0)]
    by_index={}
    for p in types:by_index.setdefault(degree-len(p),[]).append(p)
    rng=random.Random(seed)
    @lru_cache(None)
    def possible(slots,total,odd):
        if not slots:return total==0 and (group=='A' or odd)
        return total>0 and any(possible(slots-1,total-i,odd or i%2==1) for i in by_index)
    total=2*degree-2+2*genus;odd=False;four=[]
    if not possible(4,total,odd):return None
    for slots in range(4,0,-1):
        choices=[i for i in by_index if possible(slots-1,total-i,odd or i%2==1)]
        i=rng.choice(choices);four.append(rng.choice(by_index[i]));total-=i;odd|=i%2==1
    return '/'.join('.'.join(str(n) for n in p if n>1) for p in four)


def exploration_slot(lane):
    return 'exact' if lane=='exact' else 'hybrid' if lane=='hybrid' else 'small' if lane in ('small','small-symmetry') else 'exploration' if lane in ('diversity','replay') else 'mechanism'


@lru_cache(None)
def coupled_profiles(d):
    """Finite mixed central-voltage recipes, not a representation census."""
    blocks=[]
    for n in range(1,d+1):
        for m in range(1,min(6,n)+1):
            if n%m or n//m>11:continue
            base=n//m
            if m>1 and base<3:continue
            for z in range(m):
                # Fixed-fiber phases must solve 2s=z or 3s=z respectively.
                possible=any(a+2*b>=base-1 and (base==2*a or any(2*s%m==z for s in range(m))) and
                             (base==3*b or any(3*s%m==z for s in range(m)))
                             for a in range(base//2+1) for b in range(base//3+1))
                if possible:blocks.append((n,m,z))
    profiles=[]
    def visit(start,left,chosen):
        if left==0:
            if any(z for _,m,z in chosen):profiles.append(tuple(f'{n}:{m}:{z}' for n,m,z in chosen))
            return
        if len(chosen)==3:return
        for i in range(start,len(blocks)):
            block=blocks[i]
            if block[0]>left:break
            visit(i,left-block[0],chosen+(block,))
    visit(0,d,())
    return tuple(sorted(profiles,key=lambda p:(len(p)==1,len(p),p)))


def weak_symmetry_rows(epoch,known,resolved=None,include_triangle=True,include_pure=False,include_next=False,include_frontier=False,include_patterns=False,include_fibers=False,include_refinement=False,include_extensions=False,include_stability=False,include_joint=False,include_useful=False):
    """Fresh biased constructions; no old exact cursor or census claim."""
    rows=[]
    if include_frontier:
        for d in (9,10,11):
            for group in ('A','S'):
                rows.append(dict(name=f's3_exact_v1_{group}{d}',kind='sample',lane='small-symmetry',
                    group=group,d=d,family='s3-exact',pure_braid_start=True,finite_construction=True))
    for d in (9,10,11,18,25,32,34,12,13,14,15,16,17,20,21,24,26,27,28,30):
        for group in ('S','A'):
            if include_pure:
                models=('A4','S4','A5')+ (('PSL27',) if include_next and d in range(9,19) else ())
                if include_frontier and d>=12:models=('S3',)+models
                for model in models:
                    for power in ((1,2,3) if include_next else (1,)):
                        name=f'pure_v2_e{epoch:06d}_{model}_r{power}_{group}{d}' if include_next else f'pure_v1_e{epoch:06d}_{model}_{group}{d}'
                        rows.append(dict(name=name,kind='sample',
                            lane='small-symmetry' if d<12 else 'symmetry',group=group,d=d,
                            family='pure-symmetry',pure_braid_start=True,pair_model=model,pure_power=power))
            if include_frontier and d in (*range(9,19),25):
                for genus in (2,4,6):
                    if genus>d-5:continue
                    rows.append(dict(name=f'pair_glue_v1_e{epoch:06d}_g{genus}_{group}{d}',kind='sample',
                        lane='small-symmetry' if d<12 else 'symmetry',group=group,d=d,
                        family='pair-glue',unrestricted_start=True,glue_genus=genus,pair_genus_max=2))
            if include_next and d in range(9,19):
                for genus in (0,1):
                    if group=='A' and d==9 and genus==0:continue # C004/C091 complete genus-zero census
                    rows.append(dict(name=f'glue_v1_e{epoch:06d}_g{genus}_{group}{d}',kind='sample',
                        lane='small-symmetry' if d<12 else 'symmetry',group=group,d=d,
                        family='dessin-glue',unrestricted_start=True,glue_genus=genus))
            if include_triangle:
                rows.append(dict(name=f'triangle_v1_e{epoch:06d}_{group}{d}',kind='sample',
                    lane='small-symmetry' if d<12 else 'symmetry',group=group,d=d,
                    family='triangle-symmetry',braid_fixed_start=True))
            rows.append(dict(name=f'square_v1_e{epoch:06d}_{group}{d}',kind='sample',
                lane='small-symmetry' if d<12 else 'symmetry',group=group,d=d,
                family='square-symmetry',braid_fixed_start=True))
    local=[];sources=[];seen=set();fiber_sources=[]
    for cert in sorted(Path(known).glob('CERTIFIED_*weak.txt')):
        profile=Path(known).parent/'known_profiles'/f'{cert.stem}.tsv'
        fixed=True # curated controls have fixed states; live catalogs supply exact profiles
        if profile.exists():
            with profile.open() as stream:fixed=int(next(csv.DictReader(stream,delimiter='\t'))['any_braid_fixed_states'])>0
        sources.append((cert,fixed))
    if resolved is not None:
        for cert in sorted(Path(resolved).glob('CERTIFIED_*weak.txt')):
            profile=Path(resolved).parent/'classifications'/f"{cert.stem.removeprefix('CERTIFIED_').removesuffix('_weak')}.tsv"
            if profile.exists():
                with profile.open() as stream:row=next(csv.DictReader(stream,delimiter='\t'))
                fixed=int(row['any_braid_fixed_states'])>0
                if fixed or include_pure:sources.append((cert,fixed))
    for cert,fixed in sources:
        identity=certificate_identity(cert)
        if identity in seen:continue
        seen.add(identity)
        h=certificate_header(cert);d=int(h['d'])
        fiber_sources.append((cert,h['group'],d))
        for group in (('S',) if h['group']=='S' else ('A','S')):
            if include_pure:
                for power in ((1,2,3) if include_next else (1,)):
                    name=f'pure_centralizer_v2_e{epoch:06d}_r{power}_{group}_{cert.stem}' if include_next else f'pure_centralizer_v1_e{epoch:06d}_{group}_{cert.stem}'
                    local.append(dict(name=name,kind='sample',
                        lane='symmetry',group=group,d=d,family='pure-centralizer',
                        pure_braid_start=True,source=str(cert),pure_power=power))
            if fixed:
                local.append(dict(name=f'centralizer_v1_e{epoch:06d}_{group}_{cert.stem}',kind='sample',
                    lane='symmetry',group=group,d=d,family='symmetry-centralizer',
                    braid_fixed_start=True,source=str(cert)))
    local.sort(key=lambda r:('genus4' not in r['source'] if include_frontier else False,r['d']!=18,
                            0 if 'XXXX' in r['source'] else 1 if 'XXYY' in r['source'] else 2,r['d']))
    # Known S18 fixed states first; fresh square constructions include degree
    # 32 = 2*4^2, a suggested continuation of 18 = 2*3^2, not an existence theorem.
    patterns=[]
    if include_patterns:
        # New geometries first, not another pass over an exhausted S3 population.
        # Old policies keep their saved schedules and executables unchanged.
        for d in (9,10,11,18,25,12,13,14,15,16,17,20,21,24,28,30):
            for group in ('S','A'):
                lane='small-symmetry' if d<=11 else 'symmetry'
                models=('A4','S4','A5','D5','D7','D9','D11','D13') if d<=18 else (
                    ('D9','D11') if d in (24,25) else ('D13',) if d>=28 else ('D9',))
                for model in models:
                    if model.startswith('D') and int(model[1:])>=d:continue
                    powers=tuple(range(1,(int(model[1:])+1)//2)) if model in ('D9','D11','D13') else (1,2)
                    for power in powers:
                        patterns.append(dict(name=f'marked_v1_e{epoch:06d}_{model}_r{power}_{group}{d}',
                            kind='sample',lane=lane,group=group,d=d,family='marked-symmetry',
                            pure_braid_start=True,pair_model=model,pure_power=power,
                            marked_pairs='all' if model.startswith('D') else 'novel',novel_power=power>3,
                            point_orbit_max=d if d<=11 else 8))
                if d>18 and d!=25:continue
                genera=(0,1,2,3,4,5,6) if d<25 else (0,1,2,3,4,5,6,12)
                for genus in genera:
                    if genus>d-5 or (group=='A' and d==9 and genus==0):continue
                    repair=16 if d==25 and genus in (4,12) else 4 if d>=12 and genus in (3,5) else 0
                    patterns.append(dict(name=f'resistant_glue_v1_e{epoch:06d}_g{genus}_{group}{d}',
                        kind='sample',lane=lane,group=group,d=d,family='resistant-glue',
                        unrestricted_start=True,glue_genus=genus,pair_genus_max=3,
                        expand_pair_actions=True,glue_repair_steps=repair,
                        point_orbit_max=d if d<=11 else 8))
    fibers=[]
    if include_fibers:
        for cert,source_group,d in fiber_sources:
            for group in ('A','S'):
                for radius in (0,2,4):
                    if include_refinement and d==18 and source_group=='S' and group=='A' and radius==0:
                        header=certificate_header(cert)
                        if header['seed'] in ('',None):raise ValueError(f'Missing source seed: {cert}')
                        if all(sum(int(n)-1 for n in part.split('.'))%2 for part in seed_profile(header['seed'],d)[0].split('/')):
                            tag=hashlib.sha256(cert.read_bytes()).hexdigest()[:16]
                            for first in range(0,int(header['states']),8):
                                fibers.append(dict(name=f'fiber_exact_v1_{tag}_A18_s{first:03d}',
                                    kind='sample',lane='symmetry',family='pair-fiber-exact',group='A',d=18,
                                    source=str(cert),fiber_state_start=first,fiber_state_count=8,
                                    finite_construction=True,unrestricted_start=True))
                            continue
                    fibers.append(dict(name=f'fiber_v1_e{epoch:06d}_{group}_{cert.stem}_r{radius}',
                        kind='sample',lane='symmetry',family='pair-fiber',group=group,d=d,
                        source=str(cert),source_group=source_group,fiber_radius=radius,unrestricted_start=True))
        fibers.sort(key=lambda r:(r['d']!=18,r['group']!='A' if r['d']==18 else r['group']!='S',
                                   'genus4' not in r['source'],r.get('fiber_radius',-1),r['source'],r.get('fiber_state_start',0)))
        for d in (*range(9,18),18,25):
            for group in ('S','A'):
                for genus in (0,1,2,3,4,5,6,12):
                    # Only A9 genus zero has a whole-passport census here.
                    if genus>d-5 or (group=='A' and d==9 and genus==0):continue
                    for radius in ((0,2) if d in (9,10,11,18,25) else (0,)):
                        fibers.append(dict(name=f'fiber_glue_v{2 if include_refinement else 1}_e{epoch:06d}_g{genus}_r{radius}_{group}{d}',
                            kind='sample',lane='small-symmetry' if d<=11 else 'symmetry',family='cut-join-glue',
                            group=group,d=d,glue_genus=genus,pair_genus_max=3,expand_pair_actions=True,
                            fiber_radius=radius,glue_repair_steps=0,unrestricted_start=True,post_cut_base_frames=include_refinement,
                            point_orbit_max=d if d<=11 else 8))
        for row in fibers+patterns+local+rows:
            if row.get('family')!='s3-exact':row.update(screen_steps=64,skip_known_starts=True)
    extensions=[]
    if include_extensions:
        # New marked relations and both-pair deformations, not a renaming of
        # the completed common-S3 or one-cut/one-cut bridge populations.
        for d in range(9,18):
            for group in ('A','S'):
                lane='small-symmetry' if d<=11 else 'symmetry'
                bounds=dict(local_block_max=8,local_frame_cap=512,local_combination_cap=20000,
                            point_orbit_max=min(d,8),screen_steps=64,skip_known_starts=True)
                for power in (1,2,4):
                    extensions.append(dict(name=f'independent_v1_e{epoch:06d}_r{power}_{group}{d}',
                        kind='sample',lane=lane,family='independent-symmetry',group=group,d=d,
                        pure_braid_start=True,pure_power=power,**bounds))
                for genus in range(min(d-5,4)+1):
                    if group=='A' and d==9 and genus==0:continue
                    for family in ('independent-glue','double-coset-glue','coordinated-glue'):
                        extensions.append(dict(name=f'{family.replace("-","_")}_v1_e{epoch:06d}_g{genus}_{group}{d}',
                            kind='sample',lane=lane,family=family,group=group,d=d,
                            unrestricted_start=True,glue_genus=genus,pair_genus_max=2,
                            boundary_centralizer_cap=2048,coordinated_rounds=4,repair_beam=6,
                            coordinated_extra_cuts=2,post_cut_base_frames=family=='coordinated-glue',**bounds))
        for cert,_,d in fiber_sources:
            for group in ('A','S'):
                for relaxed in (False,True):
                    extensions.append(dict(name=f'coordinated_fiber_v1_e{epoch:06d}_{group}_{cert.stem}_relaxed{int(relaxed)}',
                        kind='sample',lane='symmetry',family='coordinated-fiber',group=group,d=d,
                        source=str(cert),unrestricted_start=True,coordinated_rounds=6,repair_beam=8,
                        coordinated_extra_cuts=2,allow_partition_change=relaxed,screen_steps=64,skip_known_starts=True))
    stability=[]
    if include_stability:
        for d in (9,10,11,12,13,14,15,16,17,18,25):
            profiles=list(coupled_profiles(d));rng=random.Random(20261008+epoch*1009+d);rng.shuffle(profiles)
            for i,profile in enumerate(profiles[:96 if include_joint else 24]):
                stability.append(dict(name=f'coupled_v1_e{epoch:06d}_d{d}_p{i:03d}',kind='sample',lane='small-symmetry' if d<=11 else 'symmetry',
                    d=d,group='auto',family='coupled-rigid',representation_blocks=profile,order=0,unrestricted_start=True))
        for d in (9,10,11):
            for group in ('A','S'):
                for power in (1,2):
                    block=f'broad_v1_e{epoch:06d}_{group}{d}_r{power}'
                    arms=[(8,512),(d,512),(8,2048),(d,2048)];random.Random(20261008+epoch+d+power).shuffle(arms)
                    for size,cap in arms:
                        stability.append(dict(name=f'{block}_b{size}_f{cap}',seed_key=block,kind='sample',lane='small-symmetry',
                            d=d,group=group,family='independent-symmetry',pure_braid_start=True,pure_power=power,
                            local_block_max=size,local_frame_cap=cap,local_combination_cap=100000,point_orbit_max=d,
                            screen_steps=64,skip_known_starts=True,comparison_block=block))
        for row in extensions+fibers+patterns+local+rows+stability:
            if row['kind']=='sample' and row.get('family')!='coupled-rigid':row['cache_positive_witnesses']=True
            if row.get('family') in ('independent-glue','double-coset-glue'):row['retire_boundary_plans']=True
    joint=[]
    if include_joint:
        for d in (*range(9 if include_useful else 10,18),25):
            for group in ('A','S'):
                for rounds in (0,8):
                    block=f'joint_v1_e{epoch:06d}_{group}{d}'
                    joint.append(dict(name=f'{block}_repair{rounds}',seed_key=block,kind='sample',
                        lane='small-symmetry' if d<=11 else 'symmetry',family='joint-pure',group=group,d=d,
                        pure_braid_start=True,pure_power=1,joint_repair_rounds=rounds,
                        local_block_max=min(d,11),point_orbit_max=d,local_frame_cap=2048,
                        local_combination_cap=100000,screen_steps=64,skip_known_starts=True,
                        cache_positive_witnesses=True,comparison_block=block))
        # Keep degree priority, but alternate control/repair order within blocks.
        for i in range(0,len(joint),2):
            if random.Random(20261010+epoch+i).randrange(2):joint[i],joint[i+1]=joint[i+1],joint[i]
        sources={}
        for cert in sorted(known.glob('CERTIFIED_*weak.txt')):
            h=certificate_header(cert)
            sources.setdefault(certificate_identity(cert),(cert,int(h['d'])))
        for identity,(cert,d) in sources.items():
            tag=hashlib.sha256(str(identity).encode()).hexdigest()[:16]
            for group in ('A','S'):
                types,_=seed_profile(certificate_header(cert)['seed'],d)
                if group=='A' and sum(sum(int(n)-1 for n in t.split('.'))%2==0 for t in types.split('/'))<2:continue
                for rounds in (0,8):
                    block=f'joint_source_v1_e{epoch:06d}_{group}{d}_{tag}'
                    joint.append(dict(name=f'{block}_r{rounds}',seed_key=block,kind='sample',lane='symmetry',
                        family='joint-pure',group=group,d=d,source=str(cert),joint_source_only=True,
                        joint_repair_rounds=rounds,pure_braid_start=True,pure_power=1,screen_steps=64,
                        skip_known_starts=True,cache_positive_witnesses=True,comparison_block=block))
    if include_useful:
        # Preserve policy-eleven rows on old configs. New rows search beyond
        # rejected matches and have protected unconstrained S18-style companions.
        useful=[]
        for old in joint:
            if old['joint_repair_rounds']==0:continue
            r=dict(old,name=old['name'].replace('_v1_','_v2_')+'_admit',
                   joint_admission_aware=True,joint_repair_rounds=12)
            useful.append(r)
            if r['d']<=17 or r['group']=='A':
                useful.append(dict(r,name=r['name']+'_hidden',joint_pure_basis=True,joint_hidden_fiber=True))
        # Broad-only source rows remain in the other construction families.
        joint=useful
        for r in joint+stability+extensions+fibers+patterns+local+rows:
            if r['kind']=='sample' and r.get('family')!='coupled-rigid':
                r['cache_positive_witnesses']=True
                r['construction_witness_limit']=256
    return joint+stability+extensions+fibers+patterns+local+rows


def discovery_worker_cost(row):
    return row.get('workers',2) if row['kind']=='hybrid' else 1


def discovery_round_complete(stopping,deadline_expired,lanes,dead,active,replays,hot):
    """A stop/deadline is not proof the finite round queues were drained."""
    return not stopping and not deadline_expired and set(lanes)<=set(dead) and not (active or replays or hot)


def defect_discovery_row(epoch):
    cases=[(group,d) for d in (*range(9,18),18,25) for group in ('A','S')]
    group,d=cases[epoch%len(cases)]
    return dict(name=f'defect_v1_e{epoch:06d}_{group}{d}',kind='defect',lane='diversity',
                group=group,d=d,draws=8,family='zero-defect',workers=1)


def defect_discovery_summary(attempt,row,elapsed):
    starts=[]
    for i in range(row['draws']):
        folder=attempt/f'draw_{i:06}'
        if not (folder/'DONE').is_file():raise RuntimeError('Incomplete defect draw; not a verified batch')
        path=folder/'SUMMARY.tsv'
        if not path.exists():
            if not (folder/'DONE').read_text().startswith('no_eligible_seed'):raise RuntimeError('Missing eligible defect summary')
            continue
        with path.open() as stream:records=list(csv.DictReader(stream,delimiter='\t'))
        if len(records)!=1:raise RuntimeError('Invalid defect summary cardinality')
        record=records[0]
        if record['group']!=row['group'] or int(record['degree'])!=row['d'] or int(record['genus'])<=0 or record['identity_entry']!='false':
            raise RuntimeError('Wrong defect population')
        if record['status'] not in ('closed','walk_positive_defect','positive_defect_witness','capped','timeout'):
            raise RuntimeError('Unknown defect completion status')
        if record['status']=='closed' and not (folder/'ZERO_DEFECT.txt').is_file():raise RuntimeError('Closed all-zero claim lacks a certificate')
        starts.append(record)
    values={'search mode':'zero-defect','group':f"{row['group']}_{row['d']}",'family':'zero-defect',
            'weak only':'false','accepted generating seeds':0,'solved seeds':0,'hard unsolved candidates':0,
            'defect eligible starts':len(starts),'closed zero-defect starts':sum(r['status']=='closed' for r in starts),
            'unresolved defect starts':sum(r['status'] in ('capped','timeout') for r in starts),
            'elapsed seconds':elapsed,'sampling batch finished':'true','search complete':'false',
            'coverage':'biased zero-start pilot; separate from weak discovery, not a census'}
    (attempt/'SUMMARY.txt').write_text(''.join(f'{k} = {v}\n' for k,v in values.items()))


def discovery_quotas(workers,hybrid,gpu_lanes=1,passport_mission=False,stability_mission=False,useful_mission=False):
    cpu=workers-(2*gpu_lanes if hybrid else 0)
    if passport_mission:
        if hybrid or cpu<3:raise ValueError('Small-passport mission requires at least three CPU-only workers')
        return dict(exact=cpu-1,discovery=1)
    if stability_mission:
        if cpu<5:raise ValueError('Stability mission requires five CPU search slots beside Metal helpers')
        quotas=(dict(exact=1,exploration=1,small=1,mechanism=cpu-3) if useful_mission
                else dict(exact=2,exploration=1,small=1,mechanism=cpu-4))
        if hybrid:quotas['hybrid']=2*gpu_lanes
        return quotas
    exploration=max(1,cpu//5);small=max(1,(cpu-exploration)//2)
    quotas=dict(exploration=exploration,small=small,mechanism=max(1,cpu-exploration-small))
    if hybrid:quotas['hybrid']=2*gpu_lanes
    return quotas


def reservation_allows(group,occupied,available,quotas):
    """Borrow only a fully idle reservation; active jobs retain their slots."""
    borrowed=sum(n for g,n in quotas.items() if g not in available and occupied[g]==0)
    return occupied[group]<quotas[group]+borrowed


def hybrid_discovery_rows(epoch,batch,batches):
    # Fresh named streams, never paired CPU/GPU repetitions. Exact symmetry
    # searches remain in their existing CPU lanes; this is biased sampling.
    cases=[(18,2),(18,3),(25,4)]+[(d,f) for d in range(9,18) for f in (0,1)]
    if epoch%2:cases=cases[3:]+cases[:3]
    return [dict(name=f'hybrid_e{epoch:06d}_{d}_{f}',lane='hybrid',kind='hybrid',d=d,
                 group='A/S',family=f,workers=2,batch=batch,batches=batches) for d,f in cases]


def unfinished_hybrid_rows(root):
    rows=[]
    for folder in sorted((root/'jobs').glob('hybrid_*')):
        if (folder/'DONE').exists():continue
        inputs=sorted(folder.glob('attempt_*/JOB_INPUT.json'))
        if inputs:
            row=json.loads(inputs[-1].read_text())
            if row['kind']!='hybrid' or row['name']!=folder.name:raise RuntimeError('Invalid saved Metal job identity')
            rows.append(row)
    return rows


def hybrid_discovery_summary(attempt):
    data=json.loads((attempt/'SUMMARY.json').read_text())
    if data['backend']!='gpu' or data['screen_mode']!=2 or data['shadow_mismatches']!=0:
        raise RuntimeError('Hybrid discovery must pass exact early-witness audits')
    totals=[sum(row['stats'][i] for row in data['rows']) for i in range(12)]
    values={'search mode':'hybrid discovery','group':f"A/S_{data['degree']}",'family':data['family'],
            'weak only':'true','accepted generating seeds':totals[5],
            'solved seeds':totals[6]+totals[7],'hard unsolved candidates':totals[9],
            'new verified weak components':totals[8],
            'elapsed seconds':data['pipeline_seconds']+data['shadow_audit_seconds'],
            'search complete':'false','sampling batch finished':'true',
            'coverage':'biased proposals, not a degree or passport census'}
    (attempt/'SUMMARY.txt').write_text(''.join(f'{k} = {v}\n' for k,v in values.items()))


def discovery_command(row,binary,attempt,seconds,known,source,cursor=None,seed=1):
    """One CPU worker per exact job; explicit weighted budget for Metal jobs."""
    if row['kind']=='hybrid':
        runtime=Path(row['runtime'])
        cmd=[runtime/'hybrid_pilot',runtime/'libscreen.dylib',runtime/'source/scripts/gpu/gpu_screen.metal',
             attempt,known,'gpu',row['d'],row['family'],row['batch'],row['batches'],row['workers'],seed,50000,250,2]
        if row.get('compact_inputs'):cmd.append(1)
        return list(map(relocated_path,cmd))
    if row['kind']=='defect':
        return list(map(str,[binary,'defect-mission','--mode','walk','--group',row['group'],'--d',row['d'],
            '--threads',1,'--draws',row['draws'],'--walk-steps',4096,'--rng',seed,'--seconds',2,
            '--cap',50000,'--family','low-index','--out',attempt]))
    cmd=[binary]
    if row.get('family')=='coupled-rigid':
        cmd+=['rigid-glue','--representation-blocks',','.join(row['representation_blocks']),
              '--representation-seed',seed,'--entry-order',row['order'],'--start-branch',0,
              '--max-nodes',2_000_000,'--cap',200_000,'--diagnose']
        if row.get('completed_representations'):cmd+=['--exclude-completed-representations',row['completed_representations']]
    elif row['kind']=='passport':
        cmd+=['passport-exact','--group',row['group'],'--d',row['d'],'--threads',1,
              '--family','passport','--passport',row['passport'],'--centralizer-rows',
              '--start',cursor['next'],'--count',18446744073709551615,
              '--weak-only','--all-pairs-intransitive','--steps',64,
              '--cap',200_000,'--candidate-seconds',10,'--seed',seed]
        cmd+=['--probe-witness-limit',256,'--short-weak-witnesses']
        if row.get('cache_positive_witnesses'):cmd+=['--cache-positive-witnesses']
        if row.get('positive_witnesses'):cmd+=['--positive-witnesses',row['positive_witnesses']]
    elif row['kind']=='pair':
        cmd+=['hexagon-search','--pair-blocks',','.join(map(str,row['pair_blocks'])),
              '--pair-paths',','.join(map(str,row['pair_paths'])),
              '--even-paths',','.join(map(str,row['even_paths'])),
              '--third-order',row['third_order'],'--group',row['group'],
              '--start-h',cursor['next_h'],'--start-third',cursor['next_third'],'--cap',200_000,'--quotient-h']
        if row['third_order']==2:cmd+=['--third-pairs',row['third_pairs']]
    elif row['kind']=='rigid':
        cmd+=['rigid-glue','--certificate',source,'--blocks',','.join(map(str,row['rigid_blocks'])),
              '--entry-order',row['order'],'--start-branch',cursor['next_branch'],'--cap',200_000,'--diagnose']
    else:
        cmd+=['smart-hunt','--group',row['group'],'--d',row['d'],'--threads',1,
              '--family',row.get('family','low-index'),'--steps',64,'--cap',row.get('cap',200_000),
              '--candidate-seconds',row.get('candidate_seconds',10),'--seed',seed,
              '--weak-only','--all-pairs-intransitive']
        if row['kind']=='replay':cmd+=['--replay',row['replay'],'--replay-only']
        else:
            if not (row.get('braid_fixed_start') or row.get('pure_braid_start') or row.get('unrestricted_start')):cmd+=['--no-braid-fixed']
            if row.get('pair_model'):cmd+=['--pair-model',row['pair_model']]
            if row.get('pure_power'):cmd+=['--pure-power',row['pure_power']]
            if row.get('marked_pairs'):cmd+=['--marked-pairs',row['marked_pairs']]
            if row.get('novel_power'):cmd+=['--novel-power']
            if row.get('point_orbit_max'):cmd+=['--point-orbit-max',row['point_orbit_max']]
            if row.get('expand_pair_actions'):cmd+=['--expand-pair-actions']
            if row.get('glue_repair_steps') is not None:cmd+=['--glue-repair-steps',row['glue_repair_steps']]
            if row.get('fiber_radius') is not None:cmd+=['--fiber-radius',row['fiber_radius']]
            if row.get('fiber_state_start') is not None:cmd+=['--fiber-state-start',row['fiber_state_start'],'--fiber-state-count',row['fiber_state_count']]
            if row.get('post_cut_base_frames'):cmd+=['--post-cut-base-frames']
            if row.get('screen_steps') is not None:cmd+=['--screen-steps',row['screen_steps']]
            if row.get('skip_known_starts'):cmd+=['--skip-known-starts']
            if row.get('glue_genus') is not None:cmd+=['--glue-genus',row['glue_genus']]
            if row.get('pair_genus_max') is not None:cmd+=['--pair-genus-max',row['pair_genus_max']]
            if row.get('excluded_passports'):cmd+=['--excluded-passports',row['excluded_passports']]
            for field in ('local_block_max','local_frame_cap','local_combination_cap','boundary_centralizer_cap',
                          'coordinated_rounds','repair_beam','coordinated_extra_cuts','joint_repair_rounds','construction_witness_limit'):
                if row.get(field) is not None:cmd+=['--'+field.replace('_','-'),row[field]]
            if row.get('allow_partition_change'):cmd+=['--allow-partition-change']
            if row.get('joint_source_only'):cmd+=['--joint-source-only']
            for flag in ('joint_admission_aware','joint_pure_basis','joint_hidden_fiber'):
                if row.get(flag):cmd+=['--'+flag.replace('_','-')]
            if row.get('source'):
                cmd+=['--certificate',row['source']]
                if row.get('family')=='local-switch':cmd+=['--switches',row['switches'],'--max-genus',row['max_genus']]
            if row.get('preserve'):cmd+=['--preserve-passport']
            if row.get('passport'):cmd+=['--passport',row['passport']]
        if row.get('cache_positive_witnesses'):cmd+=['--cache-positive-witnesses']
        if row.get('positive_witnesses'):cmd+=['--positive-witnesses',row['positive_witnesses']]
        if row.get('retire_boundary_plans'):cmd+=['--retire-boundary-plans']
    cmd+=['--known-certificates',known,'--seconds',seconds,'--out',attempt]
    return list(map(relocated_path,cmd))


def discovery_pair_profile(row):
    profile=dict(cycles=row['pair_blocks'],paths=row['pair_paths'],even_paths=row['even_paths'],
                 third_order=row['third_order'],target_group=row['group'],quotient_h=True)
    if row['third_order']==2:profile['third_pairs']=row['third_pairs']
    return profile


def completed_rigid_profiles(parent=None):
    rows=[]
    for name,filename in (('campaign_20260926_235453_768706','MECHANISM_JOBS.json'),
                          ('campaign_20260927_115848_035745','STRUCTURE_JOBS.json')):
        root=CAMPAIGNS/name
        if not (root/'inputs'/filename).exists():continue
        for row in json.loads((root/'inputs'/filename).read_text()):
            if not row.get('rigid_blocks'):continue
            summaries=[p for p in (root/'jobs'/row['name']).glob('attempt_*/SUMMARY.txt') if (p.parent/'VERIFIED').exists()]
            for path in summaries:
                v=dict(l.split(' = ',1) for l in path.read_text().splitlines() if ' = ' in l)
                if v.get('search complete')=='true' and v.get('hard unsolved candidates','0')=='0':
                    rows.append(dict(blocks=row['rigid_blocks'],order=int(v['entry order divides']),evidence=str(path.relative_to(PROJECT))))
                    break
    if parent:
        for root,state,plan in discovery_chain(parent):
            for row in plan['lanes'].get('rigid',[]):
                folder=root/'jobs'/row['name']
                if discovery_coverage(row,folder,state['binary_sha256'])['traversal_exhausted']:
                    rows.append(dict(blocks=row['rigid_blocks'],order=row['order'],evidence=str(folder)))
    return rows


def deduplicate_representations(schedule,keys,previous):
    old={}
    for row in previous:old.setdefault(keys[','.join(map(str,row['blocks']))],[]).append(row)
    unique=[];seen={}
    for row in schedule['lanes']['rigid']:
        representation=keys[','.join(map(str,row['rigid_blocks']))];order=row['order'];key=(representation,order)
        prior=next((r for r in old.get(representation,[]) if r['order']==0 or (order!=0 and r['order']%order==0)),None)
        if prior or key in seen:
            schedule['excluded'].append(dict(profile=row['rigid_blocks'],order=order,
                evidence=prior['evidence'] if prior else seen[key],reason='exactly conjugate braid representation with covering completed order' if prior else 'same braid representation/order already scheduled'))
            continue
        seen[key]=row['name'];row['representation_key']=representation;unique.append(row)
    schedule['lanes']['rigid']=unique
    return schedule


def rigid_glue_jobs(entry_order=2, tetrahedral=False):
    # Central quotients plus the three-cyclic-fiber B3 family (degree 3m).
    # Four product-one involutions generate a generalized dihedral group,
    # never A_d/S_d for d>=4, so require some center order >2.
    sizes=tuple(sorted({4,8,16,*range(3,35,3)}))
    if tetrahedral:
        # t8:0 and t16:2 are conjugate to existing source-derived blocks.
        sizes+=tuple(f't{4*m}:{phase}' for m in range(2,9) for phase in range(m)
                     if (m,phase) not in ((2,0),(4,2)))
        sizes=tuple(sorted(sizes,key=lambda b:(block_degree(b),str(b))))
    non_involutory={16}|{3*m for m in range(1,12) if m//math.gcd(m,2)>2}
    profiles=set()
    def extend(start, total, blocks):
        new=any(str(b).startswith('t') for b in blocks)
        center_large=non_involutory.intersection(blocks) or any(
            str(b).startswith('t') and (block_degree(b)//4)//math.gcd(block_degree(b)//4,3)>2 for b in blocks)
        if total>=4 and (not tetrahedral or new) and (entry_order!=2 or center_large):
            profiles.add(blocks)
        for i in range(start,len(sizes)):
            size=sizes[i]
            degree=block_degree(size)
            if total+degree<=34:
                extend(i,total+degree,blocks+(size,))
    extend(0,0,())
    priority=[(9,16),(4,9),(3,16),(8,16),(15,16),(15,15),(18,),(12,16),(3,3,16),(9,9,16),(16,16),(9,9,9)]
    if tetrahedral:
        priority=[(9,'t20:2'),(9,'t20:3'),(9,'t16:1'),(3,'t20:0'),(3,'t28:0'),(6,'t28:0')]
    return priority+sorted(profiles-set(priority),key=lambda b:(sum(map(block_degree,b)),tuple(map(str,b))))


def atomic_json(path, value):
    tmp = path.with_name(path.name + f".{os.getpid()}.tmp")
    tmp.write_text(json.dumps(value, indent=2) + "\n")
    tmp.replace(path)


def census_progress(root):
    path = root / "CENSUS_PROGRESS.tsv"
    if not path.exists():
        return None
    with path.open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if not rows:
        return None
    return (sum(row.get("complete") == "True" for row in rows), len(rows),
            sum(int(row.get("verified_scanned_prefix", 0)) for row in rows),
            sum(int(row.get("total_pairs", 0)) for row in rows))


def format_status_bar(state, active=False):
    root = Path(state.get("root", ""))
    progress = census_progress(root) if root.exists() else None
    if progress:
        complete, total, scanned, pairs = progress
        fraction = complete / total if total else 0
        width = 24
        filled = round(width * fraction)
        bar = "#" * filled + "." * (width - filled)
        census = f"A9 exact [{bar}] {complete}/{total} passports ({scanned:,}/{pairs:,} pairs)"
    else:
        census = "A9 exact [........................] n/a"
    jobs = state.get("completed_jobs", "?")
    lock = "running" if active else state.get("status", "unknown")
    campaign = "unbounded campaign" if state.get("forever") else "finite campaign"
    return (f"Hurwitz {census} | phase {state.get('phase', 'initializing')} | {campaign} | "
            f"chunks {jobs} | {lock} | threads {state.get('threads', '?')} | "
            f"job {state.get('job') or 'starting'}")


def stop_signal(signum, frame):
    global STOPPING
    STOPPING = True
    if CHILD is not None and CHILD.poll() is None:
        CHILD.terminate()
    for child in tuple(PAIR_CHILDREN):
        if child.poll() is None:child.terminate()


def command(cmd, logfile, env=None):
    global CHILD
    cmd=list(map(relocated_path,cmd))
    print("RUN", " ".join(map(str, cmd)), flush=True)
    with logfile.open("a") as log:
        CHILD = subprocess.Popen(list(map(str, cmd)), cwd=PROJECT, env=env,
                                 stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
        for line in CHILD.stdout:
            log.write(line)
            log.flush()
            print(line, end="", flush=True)
        rc = CHILD.wait()
    CHILD = None
    if rc:
        print(f"Command exited {rc}; see {logfile}", flush=True)
    return rc


def candidates_from_files(paths,weak_only=False):
    paths = sorted(set(paths))
    # A later strong witness or closed certificate resolves an earlier capped
    # attempt. A weak witness alone must NOT erase a pending strong question.
    globally_resolved = set()
    for folder in {path.parent for path in paths}:
        results = folder / 'RESULTS.tsv'
        if results.exists():
            with results.open() as stream:
                for row in csv.DictReader(stream,delimiter='\t'):
                    if (row.get('status') == 'strong_found' or row.get('status','').startswith('certified_')
                            or (weak_only and row.get('status')=='weak_found')):
                        globally_resolved.add(row.get('seed'))
    groups = {}
    for path in paths:
        # An interrupted smart run leaves survivors before RESULTS is written.
        # Replay those too. Successful/finished rows need no retry.
        resolved = set()
        results = path.parent / "RESULTS.tsv"
        if path.name == "survivors.txt" and results.exists():
            with results.open() as stream:
                for row in csv.DictReader(stream, delimiter="\t"):
                    if row.get("status", "").endswith("_found") or row.get("status", "").startswith("certified_"):
                        resolved.add(row.get("seed"))
        for line in path.read_text().splitlines():
            if not line.strip():
                continue
            if not line.startswith("group=") or " seed=" not in line:
                raise ValueError(f"Malformed candidate in {path}: {line[:100]}")
            group_degree = line.split()[0].removeprefix("group=")
            seed = line.split(" seed=", 1)[1]
            if seed not in resolved and seed not in globally_resolved:
                groups.setdefault(group_degree, set()).add(seed)
    return groups


def report(root,index=None):
    summaries = sorted(index.records) if index is not None else sorted(p for p in root.glob("jobs/*/attempt_*/SUMMARY.txt") if (p.parent/'VERIFIED').exists())
    defect_rows=[]
    for path in summaries:
        if path.parent.parent.name.startswith('defect_v1_'):
            for data in sorted(path.parent.glob('draw_*/SUMMARY.tsv')):
                with data.open() as stream:defect_rows.extend(dict(r,attempt=str(path.parent.relative_to(root))) for r in csv.DictReader(stream,delimiter='\t'))
    with (root/'DEFECT_STARTS.tsv').open('w') as stream:
        fields=['attempt','case_id','group','degree','genus','status','histogram_unit','measured_states','S_quotient_states','seed','witness']
        writer=csv.DictWriter(stream,fields,delimiter='\t',extrasaction='ignore');writer.writeheader();writer.writerows(defect_rows)
    generated = solved = hard = 0
    coverage = []
    scans = []
    for path in summaries:
        values = index.records[path] if index is not None else dict(line.split(" = ", 1) for line in path.read_text().splitlines() if " = " in line)
        generated += int(values.get("accepted generating seeds", 0))
        solved += int(values.get("solved seeds", 0))
        hard += int(values.get("hard unsolved candidates", 0))
        coverage.append([str(path.parent.relative_to(root)), values.get('group'), values.get('family'),values.get('search mode','smart-hunt'),values.get('passport',''),
                         values.get('weak only'), values.get('accepted generating seeds'),
                         values.get('solved seeds'), values.get('hard unsolved candidates'),
                         values.get('elapsed seconds')])
        if 'scan next' in values:
            scans.append([str(path.parent.relative_to(root)),values.get('group'),values.get('passport'),
                          values['scan start'],values['scan next'],values['scan end'],
                          values['scan total candidate pairs'],values['scan range resolved']])
    certificates = {}
    imported = {certificate_identity(p) for p in root.glob('known/CERTIFIED_*.txt')}
    discovered = index.certificate_paths if index is not None else [p for p in root.glob("jobs/*/attempt_*/CERTIFIED_*.txt") if (p.parent/'VERIFIED').exists()]
    for path in sorted(list(root.glob('known/CERTIFIED_*.txt')) + list(discovered)):
        certificates.setdefault(certificate_identity(path), path)
    if (root/'known_profiles').exists():
        known_profiles=[]
        for path in root.glob('known_profiles/*.tsv'):
            if path.name.startswith('._'):continue  # macOS AppleDouble metadata, not a research table
            with path.open() as stream:known_profiles.extend(csv.DictReader(stream,delimiter='\t'))
        def profile_passport(r):return (r['group'],r['degree'],tuple(sorted(r['passport'].split('/'))))
        known_passports={profile_passport(r) for r in known_profiles}
        known_degrees={int(r['degree']) for r in known_profiles}
        discoveries=[]
        for identity,path in certificates.items():
            if identity in imported:continue
            tag=hashlib.sha256(str(identity).encode()).hexdigest();profile=root/'classifications'/f'{tag}.tsv'
            if not profile.exists():continue
            with profile.open() as stream:r=next(csv.DictReader(stream,delimiter='\t'))
            discoveries.append(dict(r,certificate=str(path),new_passport=str(profile_passport(r) not in known_passports),
                                    smaller_than_starting_examples=str(int(r['degree'])<min(known_degrees)) if known_degrees else 'unknown'))
        discoveries.sort(key=lambda r:(r['whole_component_no_braid_fixed']!='true',r['new_passport']!='True',int(r['degree']),int(r['quotient_states'])))
        fields=['group','degree','genus','quotient_states','passport','new_passport','smaller_than_starting_examples','whole_component_no_braid_fixed','any_braid_fixed_states','first_pair_partition_histogram','certificate']
        with (root/'DISCOVERIES.tsv').open('w') as stream:
            writer=csv.DictWriter(stream,fieldnames=fields,delimiter='\t',extrasaction='ignore');writer.writeheader();writer.writerows(discoveries)
    lines = ["Hurwitz campaign report", f"Updated: {dt.datetime.now().isoformat(timespec='seconds')}",
             f"Completed job attempts: {len(summaries)}", f"Accepted generating samples: {generated:,}",
             f"Solved samples (goal depends on job): {solved:,}", f"Capped/timed-out/unprocessed events: {hard:,}",
             f"Certified bad components modulo S_d conjugation: {len(certificates)}", "",
             f"New components beyond starting catalog: {len(set(certificates)-imported)}",
             "Counts are samples, including repeats/replays; not counts of distinct tuples or braid orbits.",
             "The whole degree has NOT been exhausted. A certificate closes one quotient component.",
             "Every newly written certificate is rechecked by the verify command before its job completes.", ""]
    for filename in ('STRUCTURED_PROGRESS.tsv','MECHANISM_PROGRESS.tsv','DISCOVERY_PROGRESS.tsv'):
        path=root/filename
        if path.exists():
            with path.open() as stream:progress=list(csv.DictReader(stream,delimiter='\t'))
            lines += [f"{filename}: {sum(r['search_complete']=='True' for r in progress)}/{len(progress)} planned slices fully resolved",
                      f'Progress table: {path}','']
    slice_file=root/'inputs/EXACT_SLICES.json'
    if slice_file.exists():
        completed=0
        planned=json.loads(slice_file.read_text())
        with (root/'CENSUS_PROGRESS.tsv').open('w') as stream:
            writer=csv.writer(stream,delimiter='\t')
            writer.writerow(['passport','total_pairs','verified_scanned_prefix','verified_resolved_prefix','complete'])
            for passport,total in planned:
                prefixes=[]
                for require_resolved in (False,True):
                    intervals=sorted((int(row[3]),int(row[4])) for row in scans if row[2]==passport
                                     and (root/row[0]/'VERIFIED').exists()
                                     and (not require_resolved or row[7]=='true'))
                    cursor=0
                    for start,end in intervals:
                        if start>cursor: break
                        cursor=max(cursor,end)
                    prefixes.append(cursor)
                complete=prefixes[1]>=total
                completed+=complete
                writer.writerow([passport,total,*prefixes,complete])
        lines += [f'Fully resolved exact passports: {completed}/{len(planned)}',
                  'Census completion requires gap-free verified intervals with no unresolved cases; samples do not establish it.','']
    with (root / 'COVERAGE.tsv').open('w') as stream:
        writer = csv.writer(stream,delimiter='\t')
        writer.writerow(['attempt','group','family','mode','passport','weak_only','samples','solved','unresolved','seconds'])
        writer.writerows(coverage)
    with (root / 'PASSPORT_SCANS.tsv').open('w') as stream:
        writer=csv.writer(stream,delimiter='\t')
        writer.writerow(['attempt','group','passport','start','next','end','total_pairs','range_resolved'])
        writer.writerows(scans)
    with (root / 'CATALOG.tsv').open('w') as stream:
        writer = csv.writer(stream,delimiter='\t')
        writer.writerow(['group','degree','quotient_states','genus','passport','max_pair_order','weak_witness','certificate','seed'])
        for path in certificates.values():
            h = certificate_header(path)
            passport, genus = seed_profile(h['seed'],int(h['d']))
            writer.writerow([h['group'],h['d'],h['states'],genus,passport,h['max_pair_order'],
                             bool(h.get('weak_witness')),path,h['seed']])
            lines += [str(path)] + [f'{k}={v}' for k,v in h.items()] + [f'passport={passport}',f'genus={genus}','']
    tmp = root / "REPORT.txt.tmp"
    tmp.write_text("\n".join(lines) + "\n")
    tmp.replace(root / "REPORT.txt")


def worker(args, root):
    require_results_volume(root)
    previous_state=json.loads((root/'STATE.json').read_text()) if (root/'STATE.json').exists() else {}
    if args.sidecar:
        os.nice(10)  # one low-priority constraint worker; main campaign untouched
    # Lock descriptor is passed by the launcher and kept open for this run.
    lock = os.fdopen(args.lock_fd, "a")
    signal.signal(signal.SIGTERM, stop_signal)
    signal.signal(signal.SIGINT, stop_signal)
    state = {"root": str(root), "pid": os.getpid(), "status": "starting",
             "hours": "unbounded" if args.forever else args.hours, "forever": args.forever,
             "threads": scheduled_threads() if args.core_schedule else args.threads,
             "worker_limit": args.threads,
             "thread_policy": "4 threads 01:00-16:00 America/Chicago; 8 otherwise" if args.core_schedule else "fixed",
             "started": dt.datetime.now().isoformat(), "job": None,
             "phase": "initialization",
             "completed_jobs": len(list(root.glob("jobs/*/DONE"))),
             "discovery_round": previous_state.get('discovery_round',0),
             "diversity_epoch": previous_state.get('diversity_epoch',0),
             "hybrid_epoch": previous_state.get('hybrid_epoch',0),
             "checkpoint_file": str(root / "CHECKPOINT.json"),
             "resume_command": f"./scripts/weak/run_campaign.py {'--sidecar ' if args.sidecar else ''}--resume '{root}' --forever" + (" --core-schedule" if args.core_schedule else f" --threads {args.threads}")}
    revision=subprocess.run(['git','rev-parse','HEAD'],cwd=PROJECT,text=True,capture_output=True)
    if revision.returncode==0:
        state['git_commit']=revision.stdout.strip()
        state['git_dirty']=bool(subprocess.check_output(['git','status','--porcelain'],cwd=PROJECT,text=True).strip())

    def save():
        atomic_json(LATEST, state)
        atomic_json(root / "STATE.json", state)

    def checkpoint(reason):
        state["last_checkpoint"] = dt.datetime.now().isoformat()
        save()
        atomic_json(root / "CHECKPOINT.json", dict(
            state,
            checkpoint_reason=reason,
            recovery="Verified attempts and DONE markers are the recovery boundary; an incomplete attempt is rerun."))

    save()
    caffeinate = None
    if sys.platform == "darwin":
        caffeinate = subprocess.Popen(["caffeinate", "-i", "-w", str(os.getpid())])
    try:
        # A separate target keeps the original executable and saved campaign
        # available. Native flags apply identically to tests and production.
        env = dict(os.environ, CARGO_TARGET_DIR=str(PROJECT / "target/smart"),
                   CARGO_BUILD_JOBS=str(args.threads),
                   RUSTFLAGS=os.environ.get("RUSTFLAGS", "-C target-cpu=native"))
        for test in ('scripts/weak/test_campaign_runner.py','scripts/weak/test_campaign_cache.py'):
            if STOPPING or command([sys.executable,PROJECT / test],root / 'BUILD.log'):
                raise RuntimeError('Runner checks did not finish successfully')
        for cmd in (["cargo", "test", "--release", "--", f"--test-threads={args.threads}"], ["cargo", "build", "--release"]):
            if STOPPING or command(cmd, root / "BUILD.log", env):
                raise RuntimeError("Build/tests did not finish successfully")
        built = PROJECT / "target/smart/release/hurwitz_breaker"
        digest = hashlib.sha256(built.read_bytes()).hexdigest()
        (root / "bin").mkdir(exist_ok=True)
        binary = root / "bin" / digest
        if not binary.exists():
            shutil.copy2(built, binary)
        (root / "source").mkdir(exist_ok=True)
        for path in [PROJECT / "src/main.rs", PROJECT / "src/smart.rs", PROJECT / 'src/smart/tests.rs', PROJECT / 'src/smart/defect_mission.rs', PROJECT / 'src/smart/weak_discovery.rs', PROJECT / 'src/smart/pair_fiber.rs', PROJECT / 'src/smart/braid_frames.rs', PROJECT / 'src/smart/positive_witness.rs', PROJECT / 'src/smart/pure_stability.rs',
                     *sorted((PROJECT/'src/smart/weak_discovery').glob('*.rs')), PROJECT / 'src/hybrid.rs',Path(__file__).resolve(),
                     PROJECT / "Cargo.toml",PROJECT / 'Cargo.lock',PROJECT / 'scripts/weak/test_campaign_runner.py',PROJECT / 'scripts/weak/test_campaign_cache.py',
                     PROJECT / 'README.md',PROJECT / 'docs/research-workflow.md',
                     PROJECT / 'scripts/weak/run_joint_stability_checks.py',
                     PROJECT / 'scripts/weak/check_useful_repair.py',
                     PROJECT / 'scripts/weak/analyze_escape_words.py',
                     PROJECT / 'scripts/weak/verify_braid_symmetry.py']:
            # Source snapshots are named by content, preserving provenance on
            # resumes even if the user has since edited the checkout.
            data = path.read_bytes()
            sha = hashlib.sha256(data).hexdigest()
            destination = root / "source" / f"{sha}_{path.name}"
            if not destination.exists():
                shutil.copy2(path, destination)
            if path.name=='verify_braid_symmetry.py':verifier_script=destination
        state["binary_sha256"] = digest
        if command([binary, "selftest"], root / "BUILD.log"):
            raise RuntimeError("Selftest failed")

        config_path = root / "CONFIG.json"
        if not config_path.exists():
            atomic_json(config_path, {"master_seed": time.time_ns() & ((1 << 63) - 1),
                                      "source": "src/smart.rs", "quotient": "S_d conjugation",
                                      "plan": args.plan or 'discovery-frontier', "degrees": args.degrees,
                                      "hybrid_enabled": not args.cpu_only and sys.platform=='darwin' and args.threads>=5 and args.plan in (None,'discovery','discovery-frontier'),
                                      "gpu_lanes": args.gpu_lanes,
                                      "compact_hybrid_inputs": args.compact_hybrid_inputs,
                                      "entry_order": args.entry_order, "tetrahedral": args.tetrahedral,
                                      "early_witness": args.early_witness,
                                      "small_passport_mission": args.small_passport_mission,
                                      "stability_mission": args.stability_mission,
                                      "useful_repair_mission": args.useful_repair_mission,
                                      "hybrid_batch_sizes": args.hybrid_batch_sizes,
                                      "hybrid_batches": args.hybrid_batches,
                                     "weak_symmetry_version": (12 if args.useful_repair_mission else 11 if args.stability_mission else 10) if args.plan in (None,'discovery-frontier') else 0,
                                      "discovery_rounds": args.discovery_rounds,
                                      "inherit_discovery": str(args.inherit_discovery.resolve()) if args.inherit_discovery else None})
        config = json.loads(config_path.read_text())
        run_event(root,'session_started',pid=os.getpid(),platform=platform.uname()._asdict(),
                  python=sys.version,logical_cpus=os.cpu_count(),worker_limit=args.threads,
                  git_commit=state.get('git_commit'),git_dirty=state.get('git_dirty'),
                  binary_sha256=digest,runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  config=config,previous_status=previous_state.get('status'))
        master = config['master_seed']
        plan = config.get('plan','broad')  # preserve the original campaign schedule on resume
        degrees = config.get('degrees',args.degrees)
        entry_order = config.get('entry_order',2)
        tetrahedral = config.get('tetrahedral',False)
        hybrid_enabled=config.get('hybrid_enabled',False)  # old resumes retain CPU-only semantics
        gpu_lanes=config.get('gpu_lanes',1)
        passport_mission=config.get('small_passport_mission',False)
        stability_mission=config.get('stability_mission',False)
        useful_mission=config.get('useful_repair_mission',False)
        exact_mission=passport_mission or stability_mission
        if stability_mission:discovery_quotas(args.threads,hybrid_enabled,gpu_lanes,stability_mission=True,useful_mission=useful_mission)
        if passport_mission and (plan!='discovery-frontier' or hybrid_enabled or args.threads<3):
            raise ValueError('Small-passport mission requires discovery-frontier and at least three CPU-only workers')
        if hybrid_enabled and (args.core_schedule or not 5<=args.threads<=10):
            raise ValueError('Hybrid discovery requires a fixed shared budget of 5..10 CPU workers')
        if hybrid_enabled and args.threads<2*gpu_lanes+3:
            raise ValueError('Reserve at least three CPU construction/replay workers beside Metal helpers')
        if args.plan and args.plan != plan:
            raise ValueError('Resume preserves its saved plan; use --new to change plans')
        if plan in ('discovery','discovery-frontier') and args.core_schedule:
            raise ValueError('Discovery resumes require a fixed worker pool: --threads or --all-cores')
        state['plan'] = plan
        # Snapshot inputs once. Resume uses the same input files and job seeds.
        input_dir = root / "inputs"
        input_dir.mkdir(exist_ok=True)
        if not (input_dir / "READY").exists():
            paths = set()
            if plan not in ('hybrid-pilot','rigid-glue','s18-family','mechanism','structure','counterexample-local','counterexample-pairs','small-symmetric','structured-followup','discovery','discovery-frontier'):
                paths.update(CAMPAIGNS.glob("campaign_*/hunt/*/*/hard_candidates.txt"))
                paths.update(CAMPAIGNS.glob("smart_*/jobs/*/attempt_*/hard_candidates.txt"))
                paths.update(CAMPAIGNS.glob("smart_*/jobs/*/attempt_*/survivors.txt"))
            for group, seeds in candidates_from_files(paths).items():
                (input_dir / f"{group}.txt").write_text("\n".join(sorted(seeds)) + "\n")
            (input_dir / "READY").touch()

        # Snapshot and verify the known examples, including certificates whose
        # original survivor was resolved and thus is no longer a hard candidate.
        known_dir = root / 'known'
        known_dir.mkdir(exist_ok=True)
        if not (known_dir / 'READY').exists():
            found = {}
            found.update((p.name,p) for p in (PROJECT/'evidence/certificates').glob('CERTIFIED_*weak.txt'))
            if config.get('inherit_discovery'):
                parent=Path(config['inherit_discovery'])
                for p in parent.glob('known/CERTIFIED_*.txt'):found.setdefault(p.name,p)
                for p in parent.glob('jobs/*/attempt_*/CERTIFIED_*.txt'):
                    if (p.parent/'VERIFIED').exists():found.setdefault(p.name,p)
            for path in ([] if plan in ('rigid-glue','s18-family','mechanism','structure','counterexample-local','counterexample-pairs','small-symmetric','structured-followup','discovery','discovery-frontier') else sorted(CAMPAIGNS.glob('smart_*/jobs/*/attempt_*/CERTIFIED_*.txt'))):
                if (path.parent.parent / 'DONE').exists():
                    found.setdefault(path.name,path)
            for name,path in found.items():
                if STOPPING or command([binary,'verify','--certificate',path],root / 'KNOWN_VERIFY.log'):
                    raise RuntimeError(f'Known certificate verification failed: {path}')
                if not (known_dir / name).exists():
                    shutil.copy2(path,known_dir / name)
            (known_dir / 'READY').touch()

        deadline = None if args.forever else time.monotonic() + args.hours * 3600
        state["status"] = "running"
        checkpoint("run-start")
        if plan=='hybrid-pilot':
            state['phase']='paired CPU/Metal pilot: full validation, bounded shared buffers'
            state['job']='hybrid-pilot';checkpoint('hybrid-pilot-start')
            if args.core_schedule or not 6<=args.threads<=10:
                raise ValueError('Hybrid pilot requires a fixed budget of 6..10 CPU workers')
            pilot_dir=root/'hybrid'
            # On resume run the original orchestrator, not edited checkout code.
            script=pilot_dir/'source/scripts/gpu/hybrid_pilot.py'
            if not script.exists():script=PROJECT/'scripts/gpu/hybrid_pilot.py'
            cmd=[sys.executable,script,'--out',pilot_dir,'--max-workers',args.threads,
                 '--min-free-gb',args.min_free_gb]
            if config.get('early_witness'):cmd+=['--early-witness']
            if config.get('compact_hybrid_inputs'):cmd+=['--compact-inputs']
            cmd+=['--batch-sizes',config.get('hybrid_batch_sizes','65536,262144'),
                  '--batches',config.get('hybrid_batches',3)]
            if deadline is not None:cmd+=['--max-seconds',max(1,deadline-time.monotonic())]
            if command(cmd,root/'HYBRID.log') and not STOPPING:
                raise RuntimeError('Hybrid pilot incomplete/failed; saved batches and ledgers retained')
            if (pilot_dir/'RESULTS.json').exists():
                state['hybrid_results']=str(pilot_dir/'RESULTS.json')
                state['completed_jobs']=json.loads((pilot_dir/'RESULTS.json').read_text())['completed_arms']
            checkpoint('hybrid-pilot-end')
            return
        serial = 0

        def job(name, group, d, family="low-index", replay=None, replay_only=False, weak=False, deep=False,
                passport=None, parents=None, deep_level=0, scan=None, orbit_map=False,
                restriction=False, construction=None, rigid_blocks=None, exceptional=None,
                pair_blocks=None,pair_paths=None,position_spectrum=False,order=None,budget=None,
                even_paths=None,resumable=False,diagnose=False,operation=None,
                source_cert=None,switches=None,max_genus=None,pair_sieve=False,centralizer=False):
            nonlocal serial
            serial += 1
            job_dir = root / "jobs" / name
            if (job_dir / "DONE").exists() or (scan and (job_dir/'ENUMERATED').exists()):
                return True
            coverage=exact_coverage(job_dir,digest,passport) if centralizer else None
            if coverage and coverage['enumeration_complete']:return True
            remaining = args.job_seconds if deadline is None else int(deadline - time.monotonic())
            if STOPPING or remaining < 1:
                return False
            if free_result_space(root) < args.min_free_gb * 1024**3:
                raise RuntimeError(f'Less than {args.min_free_gb:g} GiB free; stopping with results retained')
            job_dir.mkdir(parents=True, exist_ok=True)
            previous = set(job_dir.glob("attempt_*/hard_candidates.txt"))
            previous.update(job_dir.glob("attempt_*/survivors.txt"))
            retry = candidates_from_files(previous,weak_only=weak).get(f"{group}_{d}", set())
            attempt = job_dir / f"attempt_{len(list(job_dir.glob('attempt_*'))):03d}"
            attempt.mkdir()
            if retry:
                if replay:
                    retry.update(Path(replay).read_text().splitlines())
                replay = attempt / "REPLAY.txt"
                replay.write_text("\n".join(sorted(retry)) + "\n")
            threads = 1 if orbit_map or restriction or rigid_blocks or pair_blocks is not None else scheduled_threads() if args.core_schedule else args.threads
            seconds = min(remaining, args.job_seconds if not deep else 600)
            if budget is not None:
                seconds=min(seconds,budget)
            if args.core_schedule:
                seconds = min(seconds, seconds_until_thread_change())
            if operation:
                # research's certificate loader requires the CERTIFIED_ prefix.
                source=input_dir/'CERTIFIED_A25_order12_genus12_weak.txt'
                if not source.exists():
                    shutil.copy2(PROJECT/'evidence/certificates/CERTIFIED_A25_order12_genus12_weak.txt',source)
                if operation=='independent':
                    sources=sorted(p.parent/'BRAID_SYMMETRY.json' for p in (root/'jobs/one_state_symmetry').glob('attempt_*/VERIFIED'))
                    if not sources:raise RuntimeError('Independent verification requires a verified symmetry job')
                    cmd=['uv','run','--isolated','--no-project','--with','sympy==1.14.0','python',verifier_script,'--input',sources[-1],'--out',attempt]
                elif operation=='research':
                    cmd=[binary,'research','--certificates',source,'--max-lift',1000,'--null-reps',0,
                         '--subset-max-k',3,'--ordered-max-k',3,'--subset-cap',200_000,
                         '--full-max',128,'--lanczos-steps',192,'--out',attempt]
                elif operation=='rayleigh':
                    cmd=[binary,'rayleigh','--certificate',source,'--ordered-k',3,'--iterations',4000,'--out',attempt]
                else:
                    cmd=[binary,operation,'--certificate',source,'--out',attempt]
            elif position_spectrum:
                source=input_dir/'S18_SOURCE.txt'
                if not source.exists():
                    shutil.copy2(PROJECT/'evidence/certificates/CERTIFIED_S18_genus1_weak.txt',source)
                cmd=[binary,'position-spectrum','--certificate',source,'--out',attempt]
            elif pair_blocks is not None:
                cmd=[binary,'hexagon-search','--pair-blocks',','.join(map(str,pair_blocks)),
                     '--pair-paths',','.join(map(str,pair_paths or [])),
                     '--seconds',seconds,'--cap',200_000,'--out',attempt]
                if even_paths:
                    cmd+=['--even-paths',','.join(map(str,even_paths))]
                if resumable:
                    coverage=pair_coverage(job_dir,digest,dict(cycles=pair_blocks,paths=pair_paths,even_paths=even_paths))
                    cmd+=['--start-h',coverage['next_h'],'--start-third',coverage['next_third']]
                    atomic_json(attempt/'PAIR_INPUT.json',dict(binary_sha256=digest,cycles=pair_blocks,paths=pair_paths,even_paths=even_paths))
            elif restriction or rigid_blocks:
                source = input_dir / 'BRAID_SOURCE.txt'
                if not source.exists():
                    shutil.copy2(PROJECT / 'evidence/certificates/CERTIFIED_A25_weak.txt', source)
                cmd = [binary, 'rigid-glue' if rigid_blocks else 'braid-restrict', '--certificate', source,
                       '--seconds', seconds, '--cap', 64 if rigid_blocks else 200_000, '--out', attempt]
                if rigid_blocks:
                    cmd += ['--blocks', ','.join(map(str, rigid_blocks)), '--entry-order', entry_order if order is None else order]
                    if diagnose:cmd+=['--diagnose']
            elif orbit_map:
                cmd = [binary, 'orbit-map', '--group', group, '--d', d, '--passport', passport,
                       '--seconds', seconds, '--cap', 1024, '--max-states', 250_000,
                       '--seed', (master + serial * 1000003) & ((1 << 64) - 1), '--out', attempt]
            else:
                cmd = [binary, 'passport-exact' if scan else 'smart-hunt', '--group', group, '--d', d,
                       '--threads', threads, '--seconds', seconds, '--family', family,
                       '--steps', 64, '--cap', 1_000_000*(1+min(deep_level,3)) if deep else 200_000,
                       '--candidate-seconds', min(600,60*(1+deep_level)) if deep else 10,
                       '--seed', (master + serial * 1000003) & ((1 << 64) - 1), '--out', attempt]
            if construction:
                cmd += ['--involution-pairs', construction[0], '--genus', construction[1]]
            if exceptional:
                cmd += ['--exceptional-type',exceptional]
            if scan:
                cursor, count = scan
                if centralizer:
                    cursor=coverage['next']
                    count=coverage['total'] if coverage['total'] is not None else count
                    atomic_json(attempt/'SCAN_INPUT.json',dict(binary_sha256=digest,passport=passport,enumeration='centralizer-v1'))
                    cmd+=['--centralizer-rows']
                for old in ([] if centralizer else sorted(job_dir.glob('attempt_*/SUMMARY.txt'))):
                    if (old.parent / 'VERIFIED').exists():
                        values=dict(line.split(' = ',1) for line in old.read_text().splitlines() if ' = ' in line)
                        cursor=max(cursor,int(values.get('scan next',cursor)))
                        if exact_scan_complete(values,scan):
                            (job_dir/'DONE').touch()
                            state["completed_jobs"] += 1
                            checkpoint("already-verified-full-scan")
                            return True
                        if exact_scan_enumerated(values,scan):
                            (job_dir/'ENUMERATED').touch()
                            return True
                if cursor>=scan[0]+count:
                    (job_dir/'ENUMERATED').touch()
                    checkpoint("already-enumerated-range-with-pending-leads")
                    return True
                cmd += ['--start',cursor,'--count',scan[0]+count-cursor]
            if replay and not scan and not rigid_blocks:
                cmd += ["--replay", replay]
            if replay_only:
                cmd += ["--replay-only"]
            if weak:
                cmd += ["--weak-only"]
            if passport and not orbit_map:
                cmd += ['--passport',passport]
            if parents:
                cmd += ['--parents',parents]
            if source_cert:
                cmd += ['--certificate',source_cert,'--switches',switches,'--max-genus',max_genus]
                if pair_sieve:cmd += ['--all-pairs-intransitive']
            state["threads"] = threads
            state["job"] = name
            state["phase"] = (operation if operation else "S18 position spectrum" if position_spectrum else f"degree {d} pair-graph constraints" if pair_blocks is not None else f"degree {d} rigid constraints" if rigid_blocks else "A25 braid symmetry restrictions" if restriction else "A24 orbit map" if orbit_map else f"{group}{d} exact passport scan" if scan else
                               f"{group}{d} {family} {'weak' if weak else 'strong'}")
            state["job_log"] = str(attempt / "RUN.log")
            save()
            atomic_json(attempt / "COMMAND.json", list(map(str, cmd)))
            rc = command(cmd, attempt / "RUN.log")
            if rc==0 and operation=='rayleigh':
                (attempt/'SUMMARY.txt').write_text('experiment = labeled-action exact Rayleigh bound\n'+(attempt/'RAYLEIGH_CERT.txt').read_text())
            if rc == 0 and (attempt / "SUMMARY.txt").exists():
                for certificate in sorted(attempt.glob("CERTIFIED_*.txt")):
                    if STOPPING or command([binary, "verify", "--certificate", certificate], attempt / "VERIFY.log"):
                        raise RuntimeError(f"Certificate verification incomplete/failed: {certificate}")
                (attempt / 'VERIFIED').touch()
                values=dict(line.split(' = ',1) for line in (attempt/'SUMMARY.txt').read_text().splitlines() if ' = ' in line)
                if resumable:
                    coverage=pair_coverage(job_dir,digest,dict(cycles=pair_blocks,paths=pair_paths,even_paths=even_paths))
                    atomic_json(job_dir/'COVERAGE.json',coverage)
                    done=coverage['traversal_exhausted']
                elif centralizer:
                    coverage=exact_coverage(job_dir,digest,passport)
                    atomic_json(job_dir/'COVERAGE.json',coverage)
                    done=coverage['search_complete']
                else:
                    done=not scan or exact_scan_complete(values,scan)
                if done:
                    (job_dir / "DONE").touch()
                    state["completed_jobs"] += 1
                elif scan and not centralizer and exact_scan_enumerated(values,scan):
                    (job_dir/'ENUMERATED').touch()
                checkpoint("verified-attempt")
            elif not STOPPING:
                raise RuntimeError(f"Job failed: {name}; see {attempt / 'RUN.log'}")
            # A full report rereads every completed chunk. Each job has already
            # synced its own summary, so batch report rebuilds without risking
            # resumability or certificate verification.
            if serial % 25 == 0:
                report(root)
                if centralizer:refresh()
            return not STOPPING and (deadline is None or time.monotonic() < deadline)

        if plan in ('discovery','discovery-frontier'):
            frontier=plan=='discovery-frontier'
            runtime=root/'hybrid_runtime'
            if hybrid_enabled:
                preparer=runtime/'source/scripts/gpu/hybrid_pilot.py'
                if not preparer.exists():preparer=PROJECT/'scripts/gpu/hybrid_pilot.py'
                if command([sys.executable,preparer,'--out',runtime,'--prepare-only','--early-witness',
                            '--known',known_dir],root/'HYBRID_BUILD.log'):
                    raise RuntimeError('Hybrid runtime preparation failed; no search launched')
            schedule_file=input_dir/'DISCOVERY_PLAN.json'
            if not schedule_file.exists():
                schedule=discovery_plan(frontier=frontier);parent=None
                if exact_mission:
                    exclusions=inherited_passport_exclusions(config.get('inherit_discovery')) if stability_mission else finite_passport_exclusions()
                    passport_plan=(useful_passport_plan(exclusions) if useful_mission else
                        stability_passport_plan(exclusions) if stability_mission else small_passport_plan(exclusions))
                    if useful_mission:
                        atomic_json(input_dir/'SMALL_DEGREE_SCOPE_AUDIT.json',small_degree_scope_audit(exclusions))
                    atomic_json(input_dir/'SMALL_PASSPORT_PLAN.json',passport_plan)
                    if stability_mission:
                        (input_dir/'FINITE_PASSPORT_EXCLUSIONS.tsv').write_text('# group\tdegree\tpassport\tverified provenance\n'+''.join(
                            f'{g}\t{d}\t{"/".join(p)}\t{source}\n' for (g,d,p),source in sorted(exclusions.items())))
                    schedule['lanes']['exact']=passport_plan['rows']
                    schedule['excluded']+=passport_plan['excluded']
                    schedule['allocation']=(['exact','exact']+schedule['allocation']) if stability_mission else ['exact','diversity','exact','symmetry','exact','rigid']
                    schedule['active_lanes']=list(schedule['lanes'])+['diversity','symmetry']
                    schedule['policy']='useful-repair-v1' if useful_mission else 'joint-stability-v1' if stability_mission else 'small-passports-v1'
                    schedule['scope']=passport_plan['scope']+'; separate biased construction, unrestricted discovery and replay lanes'
                if config.get('inherit_discovery'):
                    parent=Path(relocated_path(config['inherit_discovery']))
                    pending=candidates_from_files(candidate_paths(parent,verified=False),weak_only=True)
                    imported=input_dir/'INHERITED_PENDING.txt'
                    imported.write_text(''.join(f'group={group} seed={seed}\n' for group,seeds in sorted(pending.items()) for seed in sorted(seeds)))
                    schedule=inherit_discovery(schedule,parent)
                previous=completed_rigid_profiles(parent if parent and not pending else None)
                models={','.join(map(str,r['rigid_blocks'])) for r in schedule['lanes']['rigid']}
                models.update(','.join(map(str,r['blocks'])) for r in previous)
                model_file=input_dir/'SYMMETRY_MODELS.txt'
                model_file.write_text('\n'.join(sorted(models))+'\n')
                preflight=root/'symmetry_preflight';preflight.mkdir(exist_ok=True)
                # An interrupted preflight is retained and rerun into a fresh directory.
                if (preflight/'REPRESENTATIONS.tsv').exists():preflight=root/f'symmetry_preflight_{time.time_ns()}'
                if command([binary,'braid-library','--certificate',known_dir/'CERTIFIED_A25_weak.txt',
                            '--profiles',model_file,'--out',preflight],root/'PREFLIGHT.log'):
                    raise RuntimeError('Braid representation preflight failed')
                with (preflight/'REPRESENTATIONS.tsv').open() as stream:
                    keys={r['profile']:r['representation_key'] for r in csv.DictReader(stream,delimiter='\t')}
                schedule=deduplicate_representations(schedule,keys,previous)
                atomic_json(schedule_file,schedule)
            schedule=json.loads(schedule_file.read_text())
            for pinned in {row['enumeration_binary']:row['enumeration_sha256'] for rows in schedule['lanes'].values() for row in rows if row.get('enumeration_binary')}.items():
                if hashlib.sha256(Path(relocated_path(pinned[0])).read_bytes()).hexdigest()!=pinned[1]:raise RuntimeError('Inherited executable does not match its coverage hash')
            # Cache certificates separately: "known" remains the immutable
            # starting catalog, so discoveries never disappear from reports.
            resolved=root/'resolved';resolved.mkdir(exist_ok=True)
            for cert in known_dir.glob('CERTIFIED_*.txt'):
                if not (resolved/cert.name).exists():shutil.copy2(cert,resolved/cert.name)
            if frontier:
                controls=root/'known_profiles';controls.mkdir(exist_ok=True)
                for cert in known_dir.glob('CERTIFIED_*.txt'):
                    target=controls/f'{cert.stem}.tsv'
                    if not target.exists() and command([binary,'certificate-profile','--certificate',cert,'--out',target]+(['--pure-stability'] if useful_mission else []),root/'PROFILE_CONTROLS.log'):
                        raise RuntimeError('Known component profile failed')
            source=known_dir/'CERTIFIED_A25_weak.txt'
            if stability_mission and not useful_mission:
                # Small finite abstract-action pilots; these are NOT tuple
                # realizations and are kept outside the counterexample catalog.
                for cycle_type in ('5.4.1','6.3.1'):
                    action_dir=root/'fiber_action_frontier'/cycle_type
                    if action_dir.exists() and not (action_dir/'SUMMARY.txt').is_file():
                        action_dir=root/'fiber_action_frontier'/f'{cycle_type}_retry_{time.time_ns()}'
                    if not action_dir.exists():
                        if command([binary,'braid-library','--fiber-actions','--d',10,'--fiber-type',cycle_type,
                                    '--seconds',60,'--out',action_dir],root/'FIBER_ACTION.log'):
                            raise RuntimeError('Abstract fiber action pilot failed')
            recover_interrupted_candidates(root)
            state['phase']='loading immutable VERIFIED results into RAM caches';save()
            cache_started=time.monotonic()
            run_event(root,'cache_rebuild_started',scope='one-time session reconstruction; raw files remain authoritative')
            live_index=VerifiedDiscoveryIndex(root) if config.get('weak_symmetry_version',0)>=6 else None
            coverage_cache=VerifiedCoverageCache(live_index)
            proof_cache=DiscoveryProofCache(root,config.get('inherit_discovery')) if config.get('weak_symmetry_version',0)>=10 else None
            run_event(root,'cache_rebuild_finished',seconds=time.monotonic()-cache_started,
                      verified_attempts=len(live_index.attempts) if live_index is not None else None,
                      cached_candidate_files=live_index.candidate_files_read if live_index is not None else None)
            if live_index is not None:live_index.live_report()
            def pending_paths(verified=True):
                return sorted(live_index.paths) if live_index is not None and verified else candidate_paths(root,verified)
            def discovery_refresh():
                started=time.monotonic()
                pending=write_discovery_progress(root,schedule,digest,coverage_cache,live_index,exact_mission)
                run_event(root,'progress_exported',seconds=time.monotonic()-started,
                          coverage_cache_hits=coverage_cache.hits,coverage_cache_misses=coverage_cache.misses,
                          scope='scheduled rows only; detailed export is not a coverage authority')
                return pending

            def diversity_rows(round_no):
                rows=[]
                sources={}
                for cert in sorted(known_dir.glob('CERTIFIED_*.txt'))+sorted(resolved.glob('CERTIFIED_*weak.txt')):
                    sources.setdefault(certificate_identity(cert),cert)
                for cert in sources.values():
                    if hybrid_enabled and cert.name in ('CERTIFIED_S18_genus1_weak.txt','CERTIFIED_S18_XXYY_genus1_weak.txt','CERTIFIED_A25_weak.txt'):
                        continue  # these neighborhoods have a dedicated Metal lane
                    h=certificate_header(cert);group=h['group'];d=int(h['d']);_,genus=seed_profile(h['seed'],d)
                    for radius in (1,2,4,8):
                        for preserve in (False,True):
                            rows.append(dict(name=f'diversity_r{round_no:03d}_{cert.stem}_s{radius}_p{int(preserve)}',
                                             lane='diversity',kind='sample',family='local-switch',group=group,d=d,
                                             source=str(cert),switches=radius,max_genus=min(genus+2,d-5) if round_no%2==0 else d-5,
                                             preserve=preserve))
                for d in range(9 if frontier else 12,18):
                    for group in ('S','A'):
                        family='low-index' if round_no%2==0 else 'passport'
                        row=dict(name=f'diversity_r{round_no:03d}_{group}{d}_{family}',lane='diversity',kind='sample',d=d,group=group,family=family)
                        if family=='passport':
                            if frontier:
                                genus=(round_no//2+master)%(d-4)
                                passport=exploration_passport(d,group,genus,master+round_no*1009+d*7+(group=='S'))
                                if not passport:continue
                                row.update(passport=passport,target_genus=genus)
                                rows.append(row);continue
                            path=input_dir/f'DIVERSITY_PASSPORTS_{group}{d}.json'
                            if not path.exists():
                                candidates=structured_slices(group,(d,),(0,1,2))
                                candidates=[r['passport'] for r in candidates if all(d-sum(int(n)-1 for n in t.split('.'))>=3 for t in r['passport'].split('/'))]
                                atomic_json(path,candidates)
                            passports=json.loads(path.read_text())
                            if not passports:continue
                            row['passport']=passports[(master+round_no)%len(passports)]
                        rows.append(row)
                random.Random(master+round_no).shuffle(rows)
                if frontier:rows.sort(key=lambda r:(r['d'] not in (9,10,11),r['group']!='S',r['d']))
                if passport_mission:rows.sort(key=lambda r:(r['d'] not in (9,10,11),r['family']=='local-switch',r['d'],r['name']))
                return rows

            def exploration_round(round_no):
                pending=discovery_refresh()
                replays=[]
                for group_degree,seeds in sorted(pending.items()):
                    if not seeds:continue
                    group,d=group_degree.split('_');content='\n'.join(sorted(seeds))+'\n'
                    tag=hashlib.sha256(content.encode()).hexdigest()[:16]
                    path=input_dir/f'DISCOVERY_REPLAY_{group_degree}_{tag}.txt'
                    if not path.exists():path.write_text(content)
                    replays.append(dict(name=f'discover_replay_{round_no:03d}_{group_degree}_{tag}',kind='replay',lane='replay',group=group,d=int(d),replay=str(path),cap=1_000_000*min(round_no+1,4),candidate_seconds=min(600,60*(round_no+1))))
                queues={lane:iter(row for row in rows if not (root/'jobs'/row['name']/'ENUMERATED').exists()) for lane,rows in schedule['lanes'].items()}
                if exact_mission:
                    queues['exact']=PassportQueue()
                    for row in schedule['lanes']['exact']:
                        folder=root/'jobs'/row['name'];coverage=coverage_cache.get(row,folder,digest)
                        if not coverage['traversal_exhausted']:
                            queues['exact'].push(row,len(list(folder.glob('attempt_*/VERIFIED')))+row.get('prior_attempts',0))
                diversity_epoch=max(round_no,state.get('diversity_epoch',0))
                diversity=diversity_rows(diversity_epoch)
                if config.get('weak_symmetry_version',0)>=4:diversity=[defect_discovery_row(round_no)]+diversity
                queues['diversity']=iter(row for row in diversity if not (root/'jobs'/row['name']/'DONE').exists())
                symmetry_enabled=1<=config.get('weak_symmetry_version',0)<=12
                if symmetry_enabled:
                    adaptive=config['weak_symmetry_version']>=5
                    fresh=weak_symmetry_rows(round_no,known_dir,resolved,include_triangle=config['weak_symmetry_version']>=2,include_pure=config['weak_symmetry_version']>=3,include_next=config['weak_symmetry_version']>=4,include_frontier=adaptive,include_patterns=config['weak_symmetry_version']>=6,include_fibers=config['weak_symmetry_version']>=7,include_refinement=config['weak_symmetry_version']>=8,include_extensions=config['weak_symmetry_version']>=9,include_stability=config['weak_symmetry_version']>=10,include_joint=config['weak_symmetry_version']>=11,include_useful=config['weak_symmetry_version']>=12)
                    if useful_mission:
                        covered=covered_construction_targets(json.loads((input_dir/'SMALL_DEGREE_SCOPE_AUDIT.json').read_text()))
                        fresh=[r for r in fresh if r['kind']!='sample' or (r['group'],r['d']) not in covered]
                    if adaptive:
                        completed=completed_s3_populations(root,schedule.get('predecessor'))
                        if config['weak_symmetry_version']>=8:completed|=completed_fiber_validations()
                        fresh=[r for r in fresh if r['name'] not in completed]
                        exclusions=input_dir/'FINITE_PASSPORT_EXCLUSIONS.tsv'
                        if not exclusions.exists():
                            exclusions.write_text('# group\tdegree\tpassport\tverified provenance\n'+''.join(
                                f'{g}\t{d}\t{"/".join(p)}\t{source}\n' for (g,d,p),source in sorted(finite_passport_exclusions().items())))
                        for r in fresh:
                            if r['d']<=11:r['excluded_passports']=str(exclusions)
                    if passport_mission:
                        fresh=[r for r in fresh if r.get('family')=='coupled-rigid' or
                               (r.get('family')=='coordinated-fiber' and r.get('allow_partition_change'))]
                    for lane in ('small-symmetry','symmetry'):
                        queues[lane]=iter([r for r in fresh if r['lane']==lane and not (root/'jobs'/r['name']/'DONE').exists()])
                retries=unfinished_hybrid_rows(root) if hybrid_enabled else []
                retry_names={row['name'] for row in retries}
                def hybrid_rows(epoch):
                    size=int(config.get('hybrid_batch_sizes','262144').split(',')[-1])
                    for row in hybrid_discovery_rows(epoch,size,config.get('hybrid_batches',3)):
                        row['runtime']=str(runtime)
                        row['compact_inputs']=config.get('compact_hybrid_inputs',False)
                        if row['name'] not in retry_names and not (root/'jobs'/row['name']/'DONE').exists():yield row
                hybrid_epoch=max(round_no,state.get('hybrid_epoch',0))
                if hybrid_enabled:queues['hybrid']=iter(itertools.chain(retries,hybrid_rows(hybrid_epoch)))
                replay_queue=deque(row for row in replays if not (root/'jobs'/row['name']/'DONE').exists())
                hot_queue=deque()
                issued={row['name'] for row in replays}
                allocation=([] if passport_mission else ['small-symmetry','symmetry','small-symmetry','symmetry'] if symmetry_enabled else [])+schedule['allocation']+(['hybrid']*gpu_lanes if hybrid_enabled else [])
                cycle=itertools.cycle(allocation);dead=set(queues)-set(allocation);active={};launched=0;finished=0
                job_times={};last_health=time.monotonic()-60
                last_replay_check=time.monotonic()
                last_report=time.monotonic()
                quotas=discovery_quotas(args.threads,hybrid_enabled,gpu_lanes,passport_mission,stability_mission,useful_mission)
                def slot(lane):return ('exact' if lane=='exact' else 'discovery') if passport_mission else exploration_slot(lane)
                def permitted(lane):
                    if lane=='hybrid':return sum(v[0]['kind']=='hybrid' for v in active.values())<gpu_lanes and sum(discovery_worker_cost(v[0]) for v in active.values())+2<=args.threads
                    if not frontier and not hybrid_enabled:return True
                    group=slot(lane)
                    occupied={g:sum(discovery_worker_cost(v[0]) for v in active.values() if slot(v[0]['lane'])==g) for g in quotas}
                    available={slot(l) for l in queues if l not in dead}
                    if replay_queue:available.add(slot('replay'))
                    if hot_queue:available.add(slot('symmetry'))
                    return reservation_allows(group,occupied,available,quotas)
                def next_row():
                    nonlocal hybrid_epoch
                    if replay_queue and permitted('replay'):return replay_queue.popleft()
                    if hot_queue and permitted('symmetry'):return hot_queue.popleft()
                    for _ in range(len(allocation)):
                        lane=next(cycle)
                        if lane in dead:continue
                        if not permitted(lane):continue
                        row=queues[lane].pop() if lane=='exact' and exact_mission else next(queues[lane],None)
                        if row is None and lane=='hybrid' and (any(l not in dead for l in schedule['lanes']) or any(v[0]['kind']!='hybrid' for v in active.values())):
                            # Feed Metal immediately; replay's timer must not
                            # turn an exhausted sampling epoch into GPU idleness.
                            hybrid_epoch+=1;state['hybrid_epoch']=hybrid_epoch
                            queues['hybrid']=iter(hybrid_rows(hybrid_epoch));save()
                            row=next(queues['hybrid'],None)
                        if row is not None:
                            if row['kind']!='passport' and row['name'] in issued:continue
                            issued.add(row['name']);return row
                        dead.add(lane)
                    return None
                try:
                    while not STOPPING and (deadline is None or time.monotonic()<deadline):
                        while sum(discovery_worker_cost(v[0]) for v in active.values())<args.threads:
                            row=next_row()
                            if row is None:break
                            if free_result_space(root)<args.min_free_gb*1024**3:raise RuntimeError('Disk free-space limit reached')
                            folder=root/'jobs'/row['name'];folder.mkdir(parents=True,exist_ok=True)
                            attempts=len(list(folder.glob('attempt_*')))
                            budget=min(3600,30*4**min(attempts+row.get('prior_attempts',0),4))
                            if row['kind']=='passport':budget=15 if attempts+row.get('prior_attempts',0)==0 else 120
                            if row['kind']=='sample':budget=30 if round_no==0 else min(args.job_seconds,120)
                            if row.get('finite_construction'):budget=max(120,budget)
                            if row['kind']=='replay':budget=min(3600,max(120,args.job_seconds)*2**min(round_no,3))
                            seconds=budget if deadline is None else min(budget,int(deadline-time.monotonic()))
                            if seconds<1:break
                            attempt=folder/f'attempt_{attempts:03d}';attempt.mkdir();cursor=None
                            if row['kind']=='pair':
                                profile=discovery_pair_profile(row);cursor=coverage_cache.get(row,folder,digest)
                                atomic_json(attempt/'PAIR_INPUT.json',dict(binary_sha256=row.get('enumeration_sha256',digest),**profile))
                            elif row['kind']=='rigid':
                                profile=dict(blocks=row['rigid_blocks'],order=row['order'],enumeration='root-branches-v1')
                                cursor=coverage_cache.get(row,folder,digest)
                                atomic_json(attempt/'RIGID_INPUT.json',dict(binary_sha256=row.get('enumeration_sha256',digest),**profile))
                            elif row['kind']=='passport':
                                cursor=coverage_cache.get(row,folder,digest);profile=None
                                atomic_json(attempt/'SCAN_INPUT.json',dict(binary_sha256=row.get('enumeration_sha256',digest),
                                    passport=row['passport'],enumeration='centralizer-v1'))
                            else:profile=None
                            seed=(master+int(hashlib.sha256(row.get('seed_key',row['name']).encode()).hexdigest()[:16],16))&((1<<64)-1)
                            if proof_cache:
                                if row['kind']=='replay':row['cache_positive_witnesses']=True
                                if row.get('cache_positive_witnesses'):row['positive_witnesses']=proof_cache.snapshot(row['d'])
                                if row.get('family')=='coupled-rigid':row['completed_representations']=proof_cache.representation_snapshot()
                            cmd=discovery_command(row,row.get('enumeration_binary',binary),attempt,seconds,known_dir if row['kind']=='hybrid' else resolved,source,cursor,seed)
                            atomic_json(attempt/'COMMAND.json',cmd);atomic_json(attempt/'JOB_INPUT.json',row)
                            log=(attempt/'RUN.log').open('w')
                            child=subprocess.Popen(cmd,cwd=PROJECT,stdout=log,stderr=subprocess.STDOUT)
                            PAIR_CHILDREN.add(child);active[child]=(row,folder,attempt,profile,log);launched+=1
                            job_times[child.pid]=time.monotonic()
                            run_event(root,'job_started',job=row['name'],kind=row['kind'],degree=row['d'],
                                      attempt=str(attempt),pid=child.pid,seed=str(seed),seconds=seconds,
                                      binary_sha256=row.get('enumeration_sha256',digest),
                                      cursor=cursor,command=cmd,cpu_workers=discovery_worker_cost(row))
                        names=','.join(v[0]['name'] for v in active.values())
                        if state['job']!=names:
                            state.update(phase=f'discovery round {round_no}: '+('hybrid + exact CPU pool' if hybrid_enabled else 'mixed construction pool'),threads=sum(discovery_worker_cost(v[0]) for v in active.values()),job=names,
                                         job_log=str(next(iter(active.values()))[2]/'RUN.log') if active else '')
                            save()
                        if not active:break
                        for child in list(active):
                            rc=child.poll()
                            if rc is None:continue
                            row,folder,attempt,profile,log=active.pop(child);PAIR_CHILDREN.discard(child);log.close()
                            job_elapsed=time.monotonic()-job_times.pop(child.pid)
                            run_event(root,'job_exited',job=row['name'],attempt=str(attempt),returncode=rc,
                                      elapsed_seconds=job_elapsed)
                            if rc==0 and row['kind']=='hybrid':hybrid_discovery_summary(attempt)
                            if rc==0 and row['kind']=='defect':
                                for cert in attempt.glob('draw_*/ZERO_DEFECT.txt'):
                                    if command([binary,'defect-mission','--mode','verify-zero','--certificate',cert],attempt/'VERIFY_ZERO.log'):raise RuntimeError('Zero-defect verification failed')
                                defect_discovery_summary(attempt,row,job_elapsed)
                            if rc!=0 or not (attempt/'SUMMARY.txt').exists():raise RuntimeError(f"Discovery job failed: {row['name']}; see {attempt/'RUN.log'}")
                            for cert in attempt.glob('CERTIFIED_*.txt'):
                                if command([binary,'verify','--certificate',cert],attempt/'VERIFY.log'):raise RuntimeError(f'Certificate verification failed: {cert}')
                                identity=certificate_identity(cert)
                                tag=hashlib.sha256(str(identity).encode()).hexdigest()
                                cached=resolved/f'CERTIFIED_{tag}_weak.txt'
                                if not cached.exists():
                                    if frontier:
                                        profiles=root/'classifications';profiles.mkdir(exist_ok=True)
                                        if command([binary,'certificate-profile','--certificate',cert,'--out',profiles/f'{tag}.tsv'],attempt/'PROFILE.log'):raise RuntimeError('Certificate structural profile failed')
                                    shutil.copy2(cert,cached)
                                    if config.get('weak_symmetry_version',0)>=5:
                                        for followup in weak_symmetry_rows(round_no,known_dir,resolved,include_pure=True,include_next=True,include_frontier=True,include_patterns=config['weak_symmetry_version']>=6,include_fibers=config['weak_symmetry_version']>=7,include_refinement=config['weak_symmetry_version']>=8,include_extensions=config['weak_symmetry_version']>=9,include_stability=config['weak_symmetry_version']>=10,include_joint=config['weak_symmetry_version']>=11,include_useful=config['weak_symmetry_version']>=12):
                                            if followup.get('source')==str(cached) and followup['name'] not in issued:
                                                issued.add(followup['name']);hot_queue.append(followup)
                                        run_event(root,'new_component_prioritized',certificate=str(cached),queued=len(hot_queue))
                            (attempt/'VERIFIED').touch()
                            if live_index is not None:live_index.register(attempt,row)
                            if proof_cache:proof_cache.register(attempt)
                            coverage_cache.invalidate(folder)
                            if row['kind'] in ('pair','rigid','passport'):c=coverage_cache.get(row,folder,digest)
                            else:
                                values=dict(l.split(' = ',1) for l in (attempt/'SUMMARY.txt').read_text().splitlines() if ' = ' in l)
                                c=dict(traversal_exhausted=row.get('finite_construction',False) and values.get('finite construction complete')=='true',
                                       sampling_batch_finished=not row.get('finite_construction',False))
                                if row.get('finite_construction'):c['finite_population_resolved']=values.get('finite population resolved')=='true'
                            atomic_json(folder/'COVERAGE.json',c)
                            if c['traversal_exhausted']:(folder/'ENUMERATED').touch()
                            elif row['kind']=='passport':
                                queues['exact'].push(row,len(list(folder.glob('attempt_*/VERIFIED')))+row.get('prior_attempts',0));dead.discard('exact')
                            if (row['kind'] in ('sample','replay','hybrid','defect') and not row.get('finite_construction')) or c['traversal_exhausted']:
                                if not (folder/'DONE').exists():(folder/'DONE').touch();state['completed_jobs']+=1
                            finished+=1;checkpoint('verified-discovery-attempt')
                            run_event(root,'job_verified',job=row['name'],attempt=str(attempt),coverage=c,
                                      summary=dict(l.split(' = ',1) for l in (attempt/'SUMMARY.txt').read_text().splitlines() if ' = ' in l),
                                      certificates=[str(p) for p in attempt.glob('CERTIFIED_*.txt')])
                        if time.monotonic()-last_health>=60:
                            run_event(root,'resource_snapshot',load_average=os.getloadavg(),
                                      free_results_gib=free_result_space(root)/1024**3,
                                      free_internal_gib=shutil.disk_usage(PROJECT).free/1024**3,
                                      completed_jobs=state['completed_jobs'],
                                      active_jobs=[dict(job=v[0]['name'],pid=p.pid,cpu_workers=discovery_worker_cost(v[0])) for p,v in active.items()])
                            last_health=time.monotonic()
                        # Promote fresh capped leads during a long first pass,
                        # rather than waiting for thousands of pilot profiles.
                        if time.monotonic()-last_replay_check>=60:
                            # Keep sampling while long exact queues are still alive.
                            if frontier and 'diversity' in dead and any(l not in dead for l in schedule['lanes']):
                                diversity_epoch+=1;state['diversity_epoch']=diversity_epoch
                                queues['diversity']=iter(diversity_rows(diversity_epoch));dead.remove('diversity');save()
                            pending=live_index.pending() if live_index is not None else candidates_from_files(pending_paths(verified=True),weak_only=True)
                            for gd,seeds in pending.items():
                                group,d=gd.split('_');content='\n'.join(sorted(seeds))+'\n'
                                tag=hashlib.sha256(content.encode()).hexdigest()[:16]
                                name=f'discover_replay_{round_no:03d}_{gd}_{tag}'
                                if name in issued:continue
                                issued.add(name);path=input_dir/f'DISCOVERY_REPLAY_{gd}_{tag}.txt'
                                if not path.exists():path.write_text(content)
                                replay_queue.append(dict(name=name,kind='replay',lane='replay',group=group,d=int(d),replay=str(path),cap=1_000_000*min(round_no+1,4),candidate_seconds=min(600,60*(round_no+1))))
                            last_replay_check=time.monotonic()
                        if time.monotonic()-last_report>=300:
                            if live_index is not None:live_index.live_report()
                            else:discovery_refresh();report(root)
                            last_report=time.monotonic()
                        time.sleep(0.1)
                finally:
                    for child,(_,_,_,_,log) in active.items():
                        if child.poll() is None:child.terminate()
                        child.wait();PAIR_CHILDREN.discard(child);log.close()
                        run_event(root,'job_interrupted',pid=child.pid,returncode=child.returncode,
                                  attempt=str(active[child][2]),reason='stop requested' if STOPPING else 'pool interrupted')
                complete=discovery_round_complete(STOPPING,deadline is not None and time.monotonic()>=deadline,
                                                  queues,dead,active,replay_queue,hot_queue)
                return launched,finished,complete

            round_no=state['discovery_round']
            round_limit=config.get('discovery_rounds',0)
            while not STOPPING and (deadline is None or time.monotonic()<deadline) and (not round_limit or round_no<round_limit):
                _,_,complete=exploration_round(round_no)
                discovery_refresh();report(root,index=live_index)
                if not complete:
                    checkpoint('discovery-round-interrupted');break
                round_no+=1;state['discovery_round']=round_no;checkpoint('discovery-round-complete')
            discovery_refresh()
            return

        if plan == 'structured-followup':
            schedule_file=input_dir/'STRUCTURED_FOLLOWUP.json'
            if not schedule_file.exists():atomic_json(schedule_file,structured_followup())
            schedule=json.loads(schedule_file.read_text())

            def refresh():
                paths=list(root.glob('jobs/*/attempt_*/hard_candidates.txt'))+list(root.glob('jobs/*/attempt_*/survivors.txt'))
                pending=candidates_from_files(paths,weak_only=True)
                with (root/'STRUCTURED_PROGRESS.tsv').open('w') as stream:
                    writer=csv.writer(stream,delimiter='\t')
                    writer.writerow(['group','degree','genus','passport','next','candidate_total','enumeration_complete','unresolved_seeds','search_complete'])
                    for row in schedule['exact']:
                        folder=root/'jobs'/row['name']
                        c=exact_coverage(folder,digest,row['passport'])
                        local=candidates_from_files(list(folder.glob('attempt_*/hard_candidates.txt'))+list(folder.glob('attempt_*/survivors.txt')),weak_only=True)
                        n=len(local.get(f"{row['group']}_{row['d']}",set()) & pending.get(f"{row['group']}_{row['d']}",set()))
                        writer.writerow([row['group'],row['d'],row['genus'],row['passport'],c['next'],c['total'],c['enumeration_complete'],n,c['enumeration_complete'] and n==0])
                with (root/'MECHANISM_PROGRESS.tsv').open('w') as stream:
                    writer=csv.writer(stream,delimiter='\t')
                    writer.writerow(['job','degree','third_pairs','next_h','next_third','traversal_exhausted','unresolved_seeds','search_complete'])
                    for row in schedule['mechanisms']:
                        folder=root/'jobs'/row['name']
                        profile=dict(cycles=row['pair_blocks'],paths=row.get('pair_paths'),even_paths=row.get('even_paths'),third_pairs=row['third_pairs'])
                        c=pair_coverage(folder,digest,profile)
                        local=candidates_from_files(list(folder.glob('attempt_*/hard_candidates.txt'))+list(folder.glob('attempt_*/survivors.txt')),weak_only=True)
                        n=len(local.get(f"{row['group']}_{row['d']}",set()) & pending.get(f"{row['group']}_{row['d']}",set()))
                        writer.writerow([row['name'],row['d'],row['third_pairs'],c['next_h'],c['next_third'],c['traversal_exhausted'],n,c['traversal_exhausted'] and n==0])
                return pending

            def resolve(round_no):
                for gd,seeds in sorted(refresh().items()):
                    if not seeds:continue
                    group,d=gd.split('_')
                    content='\n'.join(sorted(seeds))+'\n'
                    tag=hashlib.sha256(content.encode()).hexdigest()[:12]
                    replay=input_dir/f'FOLLOWUP_REPLAY_{round_no}_{gd}_{tag}.txt'
                    if not replay.exists():replay.write_text(content)
                    if not job(f'followup_replay_{round_no}_{gd}_{tag}',group,int(d),replay=replay,
                               replay_only=True,weak=True,deep=True,deep_level=min(round_no,3)):return False
                refresh()
                return not STOPPING and (deadline is None or time.monotonic()<deadline)

            refresh()
            for phase in ('alternating','degree11'):
                for row in (r for r in schedule['exact'] if r['phase']==phase):
                    while not exact_coverage(root/'jobs'/row['name'],digest,row['passport'])['enumeration_complete']:
                        if not job(row['name'],row['group'],row['d'],'passport',weak=True,
                                   passport=row['passport'],scan=(0,row['total']),centralizer=True,
                                   budget=30 if phase=='degree11' else None):return
                        if phase=='degree11':break
                if not resolve(0 if phase=='alternating' else 1):return

            def pair_round(budget):
                # Independent symmetry profiles use one process each. The pool
                # shares the same five-worker ceiling as the exact scanner.
                pending=iter(r for r in schedule['mechanisms'] if not (root/'jobs'/r['name']/'DONE').exists())
                active={}
                try:
                    while not STOPPING and (deadline is None or time.monotonic()<deadline):
                        while len(active)<args.threads:
                            row=next(pending,None)
                            if row is None:break
                            if free_result_space(root)<args.min_free_gb*1024**3:
                                raise RuntimeError('Disk free-space limit reached')
                            folder=root/'jobs'/row['name'];folder.mkdir(parents=True,exist_ok=True)
                            profile=dict(cycles=row['pair_blocks'],paths=row.get('pair_paths'),even_paths=row.get('even_paths'),third_pairs=row['third_pairs'])
                            c=pair_coverage(folder,digest,profile)
                            attempt=folder/f"attempt_{len(list(folder.glob('attempt_*'))):03d}";attempt.mkdir()
                            seconds=budget if deadline is None else min(budget,int(deadline-time.monotonic()))
                            if seconds<1:break
                            cmd=[binary,'hexagon-search','--pair-blocks',','.join(map(str,row['pair_blocks'])),
                                 '--pair-paths',','.join(map(str,row.get('pair_paths',[]))),
                                 '--even-paths',','.join(map(str,row.get('even_paths',[]))),
                                 '--third-pairs',row['third_pairs'],
                                 '--seconds',seconds,'--cap',200_000,'--start-h',c['next_h'],'--start-third',c['next_third'],'--out',attempt]
                            atomic_json(attempt/'PAIR_INPUT.json',dict(binary_sha256=digest,**profile))
                            atomic_json(attempt/'COMMAND.json',list(map(str,cmd)))
                            log=(attempt/'RUN.log').open('w')
                            child=subprocess.Popen(list(map(str,cmd)),cwd=PROJECT,stdout=log,stderr=subprocess.STDOUT)
                            PAIR_CHILDREN.add(child);active[child]=(row,folder,attempt,profile,log)
                        names=','.join(v[0]['name'] for v in active.values())
                        if state['job']!=names:
                            state.update(phase='A15–A17 independent pair-graph constraints',threads=len(active),job=names,
                                         job_log=str(next(iter(active.values()))[2]/'RUN.log') if active else '')
                            save()
                        if not active:break
                        for child in list(active):
                            rc=child.poll()
                            if rc is None:continue
                            row,folder,attempt,profile,log=active.pop(child);PAIR_CHILDREN.discard(child);log.close()
                            if rc!=0 or not (attempt/'SUMMARY.txt').exists():
                                raise RuntimeError(f"Pair job failed: {row['name']}; see {attempt/'RUN.log'}")
                            for cert in attempt.glob('CERTIFIED_*.txt'):
                                if command([binary,'verify','--certificate',cert],attempt/'VERIFY.log'):
                                    raise RuntimeError(f'Certificate verification failed: {cert}')
                            (attempt/'VERIFIED').touch()
                            c=pair_coverage(folder,digest,profile);atomic_json(folder/'COVERAGE.json',c)
                            if c['traversal_exhausted']:
                                (folder/'DONE').touch();state['completed_jobs']+=1
                            checkpoint('verified-parallel-pair-attempt')
                        time.sleep(0.1)
                finally:
                    for child,(_,_,_,_,log) in active.items():
                        if child.poll() is None:child.terminate()
                        child.wait();PAIR_CHILDREN.discard(child);log.close()

            pair_round(30)
            round_no=2
            while not STOPPING and (deadline is None or time.monotonic()<deadline):
                if not resolve(round_no):return
                remaining=[r for r in schedule['exact'] if not exact_coverage(root/'jobs'/r['name'],digest,r['passport'])['enumeration_complete']]
                pairs_pending=any(not (root/'jobs'/r['name']/'DONE').exists() for r in schedule['mechanisms'])
                if not remaining and not pairs_pending and not refresh():return
                for row in remaining:
                    if not job(row['name'],row['group'],row['d'],'passport',weak=True,
                               passport=row['passport'],scan=(0,row['total']),centralizer=True):return
                pair_round(args.job_seconds)
                round_no+=1
            refresh()
            return

        if plan == 'counterexample-local':
            sources=sorted(known_dir.glob('CERTIFIED_*weak.txt'))
            if len(sources)<5:raise RuntimeError('Expected all five verified starting certificates')
            radii=(1,2,3,5,8,12)
            round_no=0
            while not STOPPING and (deadline is None or time.monotonic()<deadline):
                conditions=[(source,radius) for source in sources for radius in radii]
                random.Random(master+round_no).shuffle(conditions)
                for source,radius in conditions:
                    h=certificate_header(source)
                    group,d=h['group'],int(h['d'])
                    _,g=seed_profile(h['seed'],d)
                    ceiling=g+2 if round_no%4<2 else d
                    sieve=round_no%2==0
                    name=f'local_r{round_no:03d}_{source.stem}_s{radius}_g{ceiling}_p{int(sieve)}'
                    if not job(name,group,d,'local-switch',weak=True,source_cert=source,
                               switches=radius,max_genus=ceiling,pair_sieve=sieve,
                               budget=30 if round_no==0 else None):return
                # Give capped leads longer exact closure before another round.
                paths=list(root.glob(f'jobs/local_r{round_no:03d}_*/attempt_*/hard_candidates.txt'))
                if paths:
                    for group_degree,seeds in sorted(candidates_from_files(paths).items()):
                        if not seeds:continue
                        group,d=group_degree.split('_')
                        replay=input_dir/f'REPLAY_r{round_no:03d}_{group_degree}.txt'
                        replay.write_text('\n'.join(sorted(seeds))+'\n')
                        if not job(f'replay_r{round_no:03d}_{group_degree}',group,int(d),
                                   replay=replay,replay_only=True,weak=True,deep=True):return
                round_no+=1
            return

        if plan == 'counterexample-pairs':
            schedule_file=input_dir/'COUNTEREXAMPLE_PAIR_JOBS.json'
            if not schedule_file.exists():atomic_json(schedule_file,counterexample_pair_jobs())
            schedule=json.loads(schedule_file.read_text())
            for row in schedule:
                if not job(**row,budget=30):return
            while not STOPPING and (deadline is None or time.monotonic()<deadline):
                pending=[row for row in schedule if not (root/'jobs'/row['name']/'DONE').exists()]
                if not pending:return
                for row in pending:
                    if not job(**row):return
            return

        if plan == 'small-symmetric':
            schedule_file=input_dir/'SMALL_SYMMETRIC_SLICES.json'
            if not schedule_file.exists():atomic_json(schedule_file,small_symmetric_slices())
            schedule=json.loads(schedule_file.read_text())
            for i,row in enumerate(schedule):
                name=f"exact_{i:03d}_S{row['d']}"
                while not any((root/'jobs'/name/marker).exists() for marker in ('DONE','ENUMERATED')):
                    if not job(name,'S',row['d'],'passport',weak=True,passport=row['passport'],
                               scan=(0,row['total'])):return
            # Exact scans promote only capped/time-limited seeds. Resolve those
            # before exploratory controls so a finite closure is not stranded.
            paths=list(root.glob('jobs/exact_*/attempt_*/hard_candidates.txt'))
            for group_degree,seeds in sorted(candidates_from_files(paths).items()):
                if not seeds:continue
                group,d=group_degree.split('_')
                replay=input_dir/f'EXACT_REPLAY_{group_degree}.txt'
                replay.write_text('\n'.join(sorted(seeds))+'\n')
                if not job(f'exact_replay_{group_degree}',group,int(d),replay=replay,
                           replay_only=True,weak=True,deep=True):return
            round_no=0
            passports={d:[r['passport'] for r in schedule if r['d']==d] for d in (9,10)}
            while not STOPPING and (deadline is None or time.monotonic()<deadline):
                controls=[(d,family) for d in (9,10) for family in ('intransitive','low-index','mixed','uniform')]
                random.Random(master+round_no).shuffle(controls)
                for d,family in controls:
                    if not job(f'random_r{round_no:03d}_S{d}_{family}','S',d,family,weak=True):return
                for d in (9,10):
                    passport=passports[d][(master+round_no)%len(passports[d])]
                    if not job(f'passport_r{round_no:03d}_S{d}','S',d,'passport',weak=True,passport=passport):return
                round_no+=1
            return

        if plan == 'structure':
            schedule_file=input_dir/'STRUCTURE_JOBS.json'
            if not schedule_file.exists():atomic_json(schedule_file,structure_jobs())
            schedule=json.loads(schedule_file.read_text())
            # Every investigation gets a first pass before expensive profiles
            # receive more time. Pair jobs continue verified, gap-free prefixes.
            for row in schedule:
                if not job(**row,budget=30):return
            while not STOPPING and (deadline is None or time.monotonic()<deadline):
                pending=[row for row in schedule if row.get('resumable') and not (root/'jobs'/row['name']/'DONE').exists()]
                if not pending:return
                for row in pending:
                    if not job(**row):return
            return

        if plan == 'mechanism':
            schedule_file=input_dir/'MECHANISM_JOBS.json'
            if not schedule_file.exists():
                atomic_json(schedule_file,mechanism_jobs())
            schedule=json.loads(schedule_file.read_text())
            # First give every slice a bounded pass. Longer second-pass jobs
            # revisit only unresolved slices; DONE means budget consumed,
            # whereas SUMMARY's search complete records mathematical coverage.
            for pass_name,budget in [('survey',30),('deep',3600)]:
                for row in schedule:
                    prior=root/'jobs'/('survey_'+row['name'])
                    summaries=sorted(prior.glob('attempt_*/SUMMARY.txt'))
                    if pass_name=='deep' and summaries and 'search complete = true' in summaries[-1].read_text():
                        continue
                    item=dict(row);item['name']=pass_name+'_'+item['name']
                    if not job(**item,budget=budget):return
            return

        if plan == 'rigid-glue':
            schedule_file = input_dir / 'RIGID_GLUE_JOBS.json'
            if not schedule_file.exists():
                atomic_json(schedule_file, rigid_glue_jobs(entry_order,tetrahedral))
            for blocks in json.loads(schedule_file.read_text()):
                if not job('rigid_' + '_'.join(map(str, blocks)), 'A', sum(map(block_degree,blocks)), rigid_blocks=blocks):
                    return
            return

        if plan == 'degree11-involutions':
            schedule_file = input_dir / 'DEGREE11_INVOLUTION_SLICES.json'
            if not schedule_file.exists():
                atomic_json(schedule_file, degree11_involution_slices())
            for i, (group, passport, total) in enumerate(json.loads(schedule_file.read_text())):
                name = f'involution_{group}11_p{i:02d}'
                while not any((root/'jobs'/name/marker).exists() for marker in ('DONE','ENUMERATED')):
                    if not job(name, group, 11, 'passport', weak=True,
                               passport=passport, scan=(0, total)):
                        return
            return

        if plan == 's18-family':
            chunk=0
            while not STOPPING and (deadline is None or time.monotonic()<deadline):
                for group,d,m,exceptional in s18_family_jobs():
                    if not job(f's18_family_{chunk:04d}_{group}{d}',group,d,'braid-fixed',weak=True,
                               construction=(m,1),exceptional=exceptional):
                        return
                chunk+=1
            return

        if plan == 'braid-structure':
            if not job('braid_restrictions', 'A', 25, restriction=True):
                return
            schedule_file = input_dir / 'BRAID_FIXED_JOBS.json'
            if not schedule_file.exists():
                atomic_json(schedule_file, braid_fixed_jobs())
            for i, (group, d, m, genus) in enumerate(json.loads(schedule_file.read_text())):
                if not job(f'braid_{i:03d}_{group}{d}_m{m}_g{genus}', group, d,
                           'braid-fixed', weak=True, construction=(m, genus)):
                    return
            return

        if plan == 'orbit-map':
            passport = '2.2.2.2.2.2.2.2.2.2/2.2.2.2.2.2.2.2.2.2/2.2.2.2.2.2.2.2.2.2/2.2.2.3.3.4.4.4'
            chunk = 0
            while not STOPPING and (deadline is None or time.monotonic() < deadline):
                if not job(f'orbit_A24_{chunk:04d}', 'A', 24, 'passport', passport=passport, orbit_map=True):
                    return
                chunk += 1
            return

        # Resolve the existing leads first, including the old 31 A9 survivors.
        for path in sorted(input_dir.glob("*.txt")):
            if not path.stem.startswith(("A_", "S_")):
                continue
            group, degree = path.stem.split("_")
            if not job(f"replay_{group}{degree}", group, int(degree), replay=path, replay_only=True, deep=True):
                return

        if plan in ('research','genus-zero') and 9 in degrees:
            # Full fixed-first-entry censuses of the manageable A9 passports.
            # Chunks are half-open rank intervals and resume from verified cursors.
            slice_file=input_dir/'EXACT_SLICES.json'
            if not slice_file.exists():
                known=[certificate_header(p) for p in known_dir.glob('CERTIFIED_*.txt')]
                atomic_json(slice_file,exact_slices(known,plan=='genus-zero'))
            for i,(passport,total) in sorted(enumerate(json.loads(slice_file.read_text())),key=exact_priority):
                for start in range(0,total,500_000):
                    name=f'x_A9_p{i:03d}_{start:08d}'
                    while not any((root/'jobs'/name/marker).exists() for marker in ('DONE','ENUMERATED')):
                        if not job(name,'A',9,'passport',passport=passport,
                                   scan=(start,min(500_000,total-start))):
                            return
            # Refresh exact completion before building the post-census schedule.
            report(root)

        round_no = 0
        while not STOPPING and (deadline is None or time.monotonic() < deadline):
            # Resolve saved leads before spending a full sweep on new samples.
            paths = list(root.glob('jobs/*/attempt_*/hard_candidates.txt'))
            paths += list(root.glob('jobs/*/attempt_*/survivors.txt'))
            for group_degree, seeds in candidates_from_files(paths).items():
                group, degree = group_degree.split("_")
                path = input_dir / f"deep_r{round_no:03d}_{group_degree}.txt"
                path.write_text("\n".join(sorted(seeds)) + "\n")
                if not job(f"deep_r{round_no:03d}_{group}{degree}", group, int(degree), replay=path, replay_only=True, deep=True,deep_level=round_no):
                    return
            if plan in ('research','genus-zero'):
                schedule_file = input_dir / f'schedule_r{round_no:03d}.json'
                parent_file = input_dir / f'parents_r{round_no:03d}.txt'
                if not schedule_file.exists():
                    paths = {p.name:p for p in known_dir.glob('CERTIFIED_*.txt')}
                    paths.update({p.name:p for p in root.glob('jobs/*/attempt_*/CERTIFIED_*.txt')
                                  if (p.parent.parent / 'DONE').exists()})
                    known = [certificate_header(p) for p in paths.values()]
                    parent_file.write_text(''.join(f"group={h['group']}_{h['d']} seed={h['seed']}\n" for h in known))
                    resolved_a9 = []
                    progress = root/'CENSUS_PROGRESS.tsv'
                    if progress.exists():
                        with progress.open() as stream:
                            resolved_a9 = [row['passport'] for row in csv.DictReader(stream,delimiter='\t') if row['complete']=='True']
                    atomic_json(schedule_file,research_jobs(degrees,known,round_no,resolved_a9))
                schedule = json.loads(schedule_file.read_text())
                for i,(group,d,family,weak,passport) in enumerate(schedule):
                    name = f'r{round_no:03d}_p{i:03d}_{group}{d}_{family}_{"weak" if weak else "strong"}'
                    if not job(name,group,d,family,weak=weak,passport=passport,
                               parents=parent_file if family == 'deform' else None):
                        return
            else:
                for d in range(8, 17):
                    for group in ('A','S'):
                        if group == 'A' and d == 8:
                            continue
                        for family,weak in [('low-index',False),('intransitive',True),('mixed',False)]:
                            name = f'r{round_no:03d}_{group}{d}_{family}_{"weak" if weak else "strong"}'
                            if not job(name,group,d,family,weak=weak):
                                return
            round_no += 1
    except Exception as error:
        if not STOPPING:
            state["status"] = "failed"
            state["error"] = str(error)
            failure=dict(error_type=type(error).__name__,error=str(error),traceback=traceback.format_exc(),
                         phase=state['phase'],jobs=state['job'],pid=os.getpid())
            atomic_json(root/'FAILURE.json',failure)
            run_event(root,'session_failed',**failure)
            print(f"ERROR: {error}", flush=True)
            raise
    finally:
        if state["status"] != "failed":
            if state.get('plan')=='structured-followup' and 'refresh' in locals():refresh()
            state["status"] = "stopped" if STOPPING else "complete"
        state["ended"] = dt.datetime.now().isoformat()
        save()
        run_event(root,'session_ended',status=state['status'],completed_jobs=state['completed_jobs'])
        report(root,index=locals().get('live_index'))
        checkpoint("run-finished")
        if caffeinate:
            caffeinate.terminate()
        lock.close()
        print(f"{state['status'].upper()}: {root / 'REPORT.txt'}", flush=True)


def main():
    global LATEST
    parser = argparse.ArgumentParser(description=__doc__)
    duration = parser.add_mutually_exclusive_group()
    duration.add_argument("--hours", type=float)
    duration.add_argument("--days", type=float)
    duration.add_argument("--forever", action="store_true", help="run until stopped; completed work is retained")
    parser.add_argument('--plan',choices=['research','genus-zero','broad','orbit-map','degree11-involutions','braid-structure','rigid-glue','s18-family','mechanism','structure','counterexample-local','counterexample-pairs','small-symmetric','structured-followup','discovery','discovery-frontier','hybrid-pilot'],help='discovery-frontier: small degrees first; hybrid-pilot: finite paired CPU/Metal pipeline experiment')
    parser.add_argument('--inherit-discovery',type=Path,help='new discovery-frontier: retain a stopped discovery campaign\'s verified prefixes and unresolved leads')
    parser.add_argument('--discovery-rounds',type=int,default=0,help='new discovery: finite sampling rounds, 0 repeats until stopped; saved on creation')
    parser.add_argument('--preview-frontier',action='store_true',help='read-only successor symmetry job plan; does not acquire a lock, build or launch')
    parser.add_argument('--sidecar',action='store_true',help='one low-priority worker with separate lock/status; does not stop the main campaign')
    parser.add_argument('--entry-order',type=int,choices=[0,2,3,4,5,6],default=2,help='rigid-glue: require a^order=1; 0 removes order restriction; saved on creation, preserved on resume')
    parser.add_argument('--tetrahedral',action='store_true',help='rigid-glue: central-lift representations beyond the original block library; saved on creation')
    parser.add_argument('--early-witness',action='store_true',help='hybrid-pilot: reject valid short braid witnesses on CPU/GPU before full-group checks; saved on creation')
    parser.add_argument('--cpu-only',action='store_true',help='new discovery campaigns: opt out of the default shared CPU/Metal lane; old resumes keep their saved policy')
    parser.add_argument('--small-passport-mission',action='store_true',help='new CPU-only frontier: finite missing A9/S9 passports, reserve all but one workers for exact scans; saved on creation')
    parser.add_argument('--stability-mission',action='store_true',help='new frontier: joint pure constraints and missing S9/A10/S10/A11/S11 passports; two reserved exact slots, compatible with Metal; saved on creation')
    parser.add_argument('--useful-repair-mission',action='store_true',help='new frontier: admission-aware repair and complete-pure hidden outer fibers; one exact reservation, independent CPU/Metal lanes retained; implies stability-mission')
    parser.add_argument('--gpu-lanes',type=int,choices=[1,2],default=1,help='new discovery: concurrent Metal jobs, two CPU helpers each, inside the total worker budget; saved on creation')
    parser.add_argument('--compact-hybrid-inputs',action='store_true',help='new hybrid runs: input seeds + SHA256 and packed flags, with complete regeneration and CPU audit; preserves old raw files')
    parser.add_argument('--hybrid-batch-sizes',default='65536,262144',help='pilot batch sizes; discovery uses the last size, up to 1048576; saved on creation')
    parser.add_argument('--hybrid-batches',type=int,default=3,help='hybrid-pilot batches per arm, 1..128; saved on creation')
    parser.add_argument('--degrees',default='9,10,11,8,12,13,14,15,16',help='research target order, comma-separated degrees 8..16; saved on creation')
    parser.add_argument('--min-free-gb',type=float,default=5,help='stop before a job if disk free space is below this GiB threshold')
    parser.add_argument("--threads", type=int, default=os.cpu_count() or 1)
    parser.add_argument('--all-cores',action='store_true',help='explicitly use every logical CPU for this launch/resume; never changes a running campaign')
    parser.add_argument("--core-schedule", action="store_true",
                        help="use 4 threads 01:00-16:00 America/Chicago and 8 otherwise")
    parser.add_argument("--job-seconds", type=int, default=300)
    parser.add_argument("--resume", type=Path, help="resume a smart campaign directory")
    parser.add_argument("--new", action="store_true", help="start a new campaign instead of resuming a stopped one")
    parser.add_argument("--status", action="store_true")
    parser.add_argument("--status-bar", action="store_true", help="print a compact human-readable progress bar")
    parser.add_argument("--stop", action="store_true", help="stop this runner; unfinished seeds/results are retained")
    parser.add_argument("--worker", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--lock-fd", type=int, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.useful_repair_mission:args.stability_mission=True
    if args.stability_mission and (args.plan!='discovery-frontier' or args.small_passport_mission or args.core_schedule):
        parser.error('stability-mission requires discovery-frontier and cannot combine with small-passport-mission or core-schedule')
    if args.small_passport_mission and (args.plan!='discovery-frontier' or not args.cpu_only or args.threads<3 or args.core_schedule):
        parser.error('small-passport-mission requires discovery-frontier, cpu-only and at least three fixed workers')
    if args.discovery_rounds<0:parser.error('discovery-rounds must be nonnegative')
    if args.discovery_rounds and args.plan not in (None,'discovery','discovery-frontier'):parser.error('discovery-rounds requires a discovery plan')
    if args.all_cores:
        if args.core_schedule or args.sidecar:parser.error('all-cores cannot be combined with core-schedule or sidecar')
        args.threads=os.cpu_count() or 1
    if args.stability_mission:
        try:discovery_quotas(args.threads,not args.cpu_only and sys.platform=='darwin',args.gpu_lanes,stability_mission=True,useful_mission=args.useful_repair_mission)
        except ValueError as e:parser.error(str(e))
    try:hybrid_sizes=list(map(int,args.hybrid_batch_sizes.split(',')))
    except ValueError:parser.error('hybrid-batch-sizes must be comma-separated integers')
    if not hybrid_sizes or any(not 1<=n<=1048576 for n in hybrid_sizes) or not 1<=args.hybrid_batches<=128:
        parser.error('hybrid batch sizes must be 1..1048576 and batches 1..128')
    if (args.new or args.preview_frontier) and args.plan in (None,'discovery','discovery-frontier') and not args.cpu_only and sys.platform=='darwin':
        if args.core_schedule or not 5<=args.threads<=10:parser.error('new hybrid discovery requires a fixed budget of 5..10 workers; use --cpu-only for other budgets')
        if args.threads<2*args.gpu_lanes+3:parser.error('gpu-lanes requires room for two helpers per lane plus three CPU construction/replay workers')
    if args.preview_frontier:
        if not 1<=args.threads<=128:parser.error('threads must be in 1..128')
        rows=weak_symmetry_rows(0,PROJECT/'evidence/certificates',include_pure=True,include_next=True,include_frontier=True,include_patterns=True,include_fibers=True,include_refinement=True,include_extensions=True,include_stability=True,include_joint=args.stability_mission,include_useful=args.useful_repair_mission)
        cpu=max(1,args.threads-(0 if args.cpu_only or sys.platform!='darwin' else 2*args.gpu_lanes))
        print(json.dumps(dict(policy=12 if args.useful_repair_mission else 11 if args.stability_mission else 10,rows=rows,symmetry_jobs=len(rows),cpu_workers=cpu,
            balanced_sampling_minutes=round(sum(not r.get('finite_construction',False) for r in rows)*30/cpu/60,2),
            finite_s3_jobs=sum(r.get('family')=='s3-exact' for r in rows),
            finite_fiber_shards=sum(r.get('family')=='pair-fiber-exact' for r in rows),
            scope='curated-catalog sampling budgets only; excludes finite enumeration, build, diversity, replay, Metal drain and soft overruns; not a discovery-time forecast'),indent=2))
        return
    if args.plan in ('discovery','discovery-frontier') and args.core_schedule:parser.error('discovery uses a fixed shared worker pool; use --threads or --all-cores')
    if args.inherit_discovery:
        if args.plan!='discovery-frontier':parser.error('inherit-discovery requires discovery-frontier')
        for filename in ('CONFIG.json','STATE.json','inputs/DISCOVERY_PLAN.json'):
            if not (args.inherit_discovery/filename).is_file():parser.error(f'Missing predecessor {filename}')
    if args.sidecar:
        if args.plan not in (None, 'rigid-glue','mechanism','structure','counterexample-pairs') or args.core_schedule:
            parser.error('sidecar supports only rigid-glue, mechanism, structure or counterexample-pairs with one worker')
        if not (args.status or args.status_bar or args.stop or args.resume or args.worker) and args.plan not in ('rigid-glue','mechanism','structure','counterexample-pairs'):
            parser.error('sidecar launches require an allowed plan')
        args.threads = 1
        LATEST = CONTROL / 'SIDECAR_LATEST.json'
        existing = args.resume or args.worker
        if existing and (existing / 'CONFIG.json').exists():
            if json.loads((existing / 'CONFIG.json').read_text()).get('plan') not in ('rigid-glue','mechanism','structure','counterexample-pairs'):
                parser.error('sidecar cannot resume or modify a main campaign')
    args.hours = None if args.forever else args.days*24 if args.days is not None else args.hours if args.hours is not None else 10
    try:
        args.degrees = list(dict.fromkeys(int(x) for x in args.degrees.split(',')))
    except ValueError:
        parser.error('degrees must be comma-separated integers')
    if not args.degrees or any(d<8 or d>16 for d in args.degrees):
        parser.error('degrees must be in 8..16 (packed engine limit)')
    if not math.isfinite(args.min_free_gb) or args.min_free_gb<=0:
        parser.error('min-free-gb must be finite and positive')
    if ((not args.forever and (not math.isfinite(args.hours) or args.hours <= 0))
            or not (1 <= args.threads <= 128) or args.job_seconds <= 0):
        parser.error("hours/job-seconds must be positive; threads must be in 1..128")
    if args.worker:
        worker(args, args.worker.resolve())
        return
    CONTROL.mkdir(parents=True, exist_ok=True)
    historical_status=CAMPAIGNS/LATEST.name
    if not LATEST.exists() and historical_status.exists():shutil.copyfile(historical_status,LATEST)
    state = json.loads(LATEST.read_text()) if LATEST.exists() else {}
    lock = (CONTROL / (".rigid.lock" if args.sidecar else ".campaign.lock")).open("a")
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        active = False
    except BlockingIOError:
        active = True
    if args.status_bar:
        print(format_status_bar(state, active))
        return
    if args.status:
        status=dict(state,lock_active=active)
        if state.get('plan')=='hybrid-pilot' and state.get('root'):
            pilot=Path(state['root'])/'hybrid'
            for name,key in [('STATUS.json','hybrid_progress'),('RESULTS.json','hybrid_results')]:
                if (pilot/name).exists():status[key]=json.loads((pilot/name).read_text())
        print(json.dumps(status, indent=2))
        return
    if args.stop:
        if active and state.get("pid"):
            cmd = subprocess.check_output(["ps", "-p", str(state["pid"]), "-o", "command="], text=True).strip()
            if str(Path(__file__).resolve()) not in cmd or "--worker" not in cmd:
                raise RuntimeError("Recorded PID does not identify this runner; refusing to stop it")
            os.kill(state["pid"], signal.SIGTERM)
            print("Stop requested. Completed jobs and partial candidate files are retained.")
        else:
            print("No campaign is active.")
        return
    if active:
        print(f"Already running: {state.get('root', 'starting')}")
        print("Use ./scripts/weak/run_campaign.py --status for progress.")
        return
    # Do not unknowingly compete with the older eight-thread campaign.
    process_text = subprocess.check_output(["ps", "-axo", "command="], text=True)
    if any(line.strip().startswith("./target/release/hurwitz_breaker ") for line in process_text.splitlines()):
        raise RuntimeError("The original hunt is still running. Stop its campaign or wait for it to finish before launching this one.")
    if args.resume:
        root = Path(relocated_path(args.resume)).resolve()
        if not (root / "CONFIG.json").exists():
            parser.error("--resume must name an existing smart campaign")
    elif not args.new and state.get("root") and state.get("status") != "complete":
        root = Path(relocated_path(state["root"]))
    else:
        root = RESULTS / ("campaign_" + dt.datetime.now(CORE_TIMEZONE).strftime("%Y%m%d_%H%M%S_%f"))
    require_results_volume(root)
    root.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, str(Path(__file__).resolve()), "--worker", str(root),
           "--lock-fd", str(lock.fileno()),
           "--threads", str(args.threads), "--job-seconds", str(args.job_seconds),
           '--degrees',','.join(map(str,args.degrees)),'--min-free-gb',str(args.min_free_gb),
           '--entry-order',str(args.entry_order)]
    cmd += ['--forever'] if args.forever else ['--hours',str(args.hours)]
    if args.plan:
        cmd += ['--plan',args.plan]
    if args.inherit_discovery:
        cmd += ['--inherit-discovery',str(args.inherit_discovery.resolve())]
    cmd += ['--discovery-rounds',str(args.discovery_rounds)]
    if args.sidecar:
        cmd += ['--sidecar']
    if args.tetrahedral:
        cmd += ['--tetrahedral']
    if args.early_witness:
        cmd += ['--early-witness']
    if args.cpu_only:
        cmd += ['--cpu-only']
    if args.small_passport_mission:cmd+=['--small-passport-mission']
    if args.stability_mission:cmd+=['--stability-mission']
    if args.useful_repair_mission:cmd+=['--useful-repair-mission']
    cmd += ['--gpu-lanes',str(args.gpu_lanes)]
    if args.compact_hybrid_inputs:cmd+=['--compact-hybrid-inputs']
    cmd += ['--hybrid-batch-sizes',args.hybrid_batch_sizes,'--hybrid-batches',str(args.hybrid_batches)]
    if args.core_schedule:
        cmd += ['--core-schedule']
    with (root / "RUN.log").open("a") as log:
        child = subprocess.Popen(cmd, cwd=PROJECT, stdin=subprocess.DEVNULL, stdout=log,
                                 stderr=subprocess.STDOUT, start_new_session=True, pass_fds=(lock.fileno(),))
    lock.close()
    # The worker writes its own authoritative state. Wait briefly to make the
    # initial status command useful; a second launcher is excluded by the lock.
    for _ in range(20):
        if (root / "STATE.json").exists() or child.poll() is not None:
            break
        time.sleep(0.1)
    print(f"Started/resumed PID {child.pid}: {root}")
    print(f"Follow: tail -f '{root / 'RUN.log'}'")
    selector = ' --sidecar' if args.sidecar else ''
    print(f"Status: ./scripts/weak/run_campaign.py{selector} --status\nStop:   ./scripts/weak/run_campaign.py{selector} --stop")


if __name__ == "__main__":
    main()
