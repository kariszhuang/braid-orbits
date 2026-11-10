#!/usr/bin/env python3
"""Targeted genus-three discovery. Plan by default; one bounded, audited pool.

Closed smoothness is a combinatorial certificate, not an isotriviality verdict.
No publisher, backup deletion, GPU exclusion or old-campaign restart is here.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
import csv
from datetime import datetime
import fcntl
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time

from audit_guided_native import audit_guided_native, audit_component_native, remaining_timeout
from run_forest_fast import compiled_source_files, spawn_launchd, storage_available
from run_smooth_degree9 import now

PROJECT = Path(__file__).resolve().parents[2]
BOUNDARY_TOTAL = 504000


def absolute_deadline(text):
    """Explicit UTC offset required; never reinterpret the host timezone."""
    if text is None:
        return None
    value = datetime.fromisoformat(text.replace('Z', '+00:00'))
    if value.utcoffset() is None:
        raise ValueError('deadline requires an explicit UTC offset')
    return value.timestamp()


def inherited_deadline(requested, predecessor_plan):
    """An overnight successor may shorten, but never extend, its deadline."""
    previous = predecessor_plan.get('deadline_epoch')
    deadlines = [v for v in (requested, previous) if v is not None]
    return min(deadlines) if deadlines else None


def finalized_predecessor(path):
    """Only drain-complete checkpoints may seed a different run."""
    path = Path(path)
    state = json.loads((path / 'STATE.json').read_text())
    final = json.loads((path / 'FINISHED.json').read_text())
    status = state.get('status')
    if status != final.get('status') or status not in ('stopped', 'complete'):
        raise ValueError('predecessor must have finalized its stop and audits')
    if state.get('active_jobs'):
        raise ValueError('predecessor still has active jobs')
    if status == 'complete' and final.get('reason') not in ('time_budget', 'absolute_deadline_drain'):
        raise ValueError('predecessor completion reason is not a safe drain')
    plan=json.loads((path/'PLAN.json').read_text())
    if plan.get('recovery_adapter'):
        provenance=json.loads((path/'ADAPTER_PROVENANCE.json').read_text());record_path=Path(provenance['recovery_receipt'])
        if sha(record_path)!=provenance['recovery_receipt_sha256']:raise ValueError('recovery receipt binding changed')
        record=json.loads(record_path.read_text())
        if not(record.get('passed')is True and record.get('all_durable_leads_resolved')is True and record.get('pending')==0 and record.get('finalized_attention')is True):
            raise ValueError('recovery is not independently complete')
        original=Path(provenance['original_root'])
        for name,digest in provenance['original_receipts'].items():
            if sha(original/name)!=digest:raise ValueError('original recovery attention receipt changed')
        for name,digest in json.loads((path/'BINDINGS.json').read_text()).items():
            if sha(path/name)!=digest:raise ValueError('recovery adapter binding changed')
    return state, final


def successor_deadline(requested, predecessor, clear_reason=None, timestamp=None):
    """Default inheritance is immutable; an expired mission needs NEW authority.

    The explicit command-line reason is persisted in the successor. This never
    edits the old plan/checkpoint or silently extends a living deadline.
    """
    if predecessor is None:
        if clear_reason is not None:
            raise ValueError('clearing an expired deadline requires a predecessor')
        return requested, None
    path = Path(predecessor)
    previous = json.loads((path / 'PLAN.json').read_text())
    if clear_reason is None:
        return inherited_deadline(requested, previous), None
    if not clear_reason.strip():
        raise ValueError('fresh human authorization reason is required')
    state, final = finalized_predecessor(path)
    old_deadline = previous.get('deadline_epoch')
    current = time.time() if timestamp is None else timestamp
    if (old_deadline is None or old_deadline > current or
        final.get('reason') != 'absolute_deadline_drain' or state['status'] != 'complete'):
        raise ValueError('only a completed, expired absolute-deadline mission may be reset')
    checkpoint = json.loads((path / 'CHECKPOINT.json').read_text())
    if (checkpoint.get('deadline_epoch') != old_deadline or
        checkpoint.get('next_epochs') != final.get('epochs') or
        checkpoint.get('completed_jobs') != final.get('completed_jobs')):
        raise ValueError('expired predecessor checkpoint differs from final receipt')
    return requested, dict(reason=clear_reason.strip(), predecessor=str(path.resolve()),
        old_deadline_epoch=old_deadline, new_deadline_epoch=requested,
        authorized_at_epoch=current, finished_sha256=sha(path / 'FINISHED.json'),
        checkpoint_sha256=sha(path / 'CHECKPOINT.json'))


def deadline_dispatch_closed(plan, timestamp=None):
    deadline = plan.get('deadline_epoch')
    if deadline is None:
        return False
    # Reserve a full proposal job, its last closure, and two minutes for the
    # independent audit. Interrupted raw leads are never silently discarded.
    reserve = plan.get('job_seconds', 45) + plan.get('closure_seconds', 60) + 120
    return (time.time() if timestamp is None else timestamp) >= deadline - reserve


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write(path, value):
    """Atomic mutable pointer; immutable receipts always get distinct paths."""
    path = Path(path); tmp = path.with_name(path.name + '.tmp')
    with tmp.open('w') as stream:
        json.dump(value, stream, indent=2); stream.write('\n')
        stream.flush(); os.fsync(stream.fileno())
    tmp.replace(path)


def identity(n, key):
    return hashlib.sha256(json.dumps([n, key], separators=(',', ':')).encode()).hexdigest()


def admit_elementary_reductions(path):
    """Fail closed on missing/changed coverage; geometric filters are unrelated."""
    path=Path(path);r=json.loads(path.read_text())
    if r.get('schema')!=1 or r.get('source_genus')!=3 or r.get('identity_free')is not True or r.get('excluded_degrees')!=[11]:
        raise ValueError('elementary exclusion scope')
    for key in ('checks','mathematics'):
        item=r[key]
        if sha(item['path'])!=item['sha256']:raise ValueError('elementary exclusion binding changed')
    if not json.loads(Path(r['checks']['path']).read_text()).get('passed'):raise ValueError('elementary reduction checks incomplete')
    profiles=set()
    for item in r['finite_degree11_audits']:
        if sha(item['path'])!=item['sha256']:raise ValueError('finite coverage binding changed')
        checked=json.loads(Path(item['path']).read_text());c=checked['coverage']
        if not checked.get('passed')or checked.get('status')!='complete_negative_slice'or c['degree']!=11 or c['pending']!=0:
            raise ValueError('finite exclusion is not complete')
        profile=tuple(c['ordered_indices'])
        if profile in profiles:raise ValueError('duplicate finite scope')
        profiles.add(profile)
    if profiles!={(1,9,8,8),(2,8,8,8)}:raise ValueError('missing finite degree11 scope')
    return dict(path=str(path.resolve()),sha256=sha(path),**r)


def admit_transposition_reduction(path):
    path=Path(path);r=json.loads(path.read_text())
    if not(r.get('schema')==1 and r.get('passed')is True and r.get('theorem')=='no-positive-genus-all-smooth-transposition'
        and r.get('source_genus')=='positive'and r.get('transitive')is True and r.get('product_one')is True
        and r.get('production_target_genus')==3 and r.get('conclusion')=='not_all_smooth'
        and r.get('reviewed_by_parent_and_geometry_against_primary')is True):raise ValueError('transposition reduction scope')
    for key in('mathematics','primary_source'):
        if sha(r[key]['path'])!=r[key]['sha256']:raise ValueError('transposition proof binding changed')
    return dict(path=str(path.resolve()),sha256=sha(path),**r)


def admit_degree_twelve(path):
    path=Path(path);r=json.loads(path.read_text())
    expected='Complete degree12, source genus3 full all-smooth exclusion, using C225 global transposition obstruction'
    if not(r.get('passed')is True and r.get('scope')==expected and set(r.get('claims',[]))=={'C226','C227'}and r.get('target_degree_lower_bound')==13):
        raise ValueError('degree12 exclusion scope')
    if r.get('single_three',{}).get('unresolved')!=0 or r.get('double_transposition',{}).get('unresolved')!=0:
        raise ValueError('degree12 exclusion is incomplete')
    items=r.get('bindings',[])
    required={'DEGREE_TWELVE.md','exact_d12_three.py','d12_three_audit.rs','d12_double.rs','audit_d12_double.py'}
    if len(items)<12 or not required<={Path(v['path']).name for v in items}:raise ValueError('degree12 coverage bindings missing')
    for item in items:
        item_path=Path(item['path']);item_path=item_path if item_path.is_absolute()else PROJECT/item_path
        if sha(item_path)!=item['sha256']:raise ValueError('degree12 proof binding changed')
    return dict(path=str(path.resolve()),sha256=sha(path),**r)


def admit_degree_thirteen(path):
    path=Path(path);r=json.loads(path.read_text())
    if not(r.get('passed')is True and r.get('claim')=='C229'and r.get('source_genus')==3 and r.get('degree')==13
        and r.get('remaining_index_budgets')==[[2,9,9,10]]and set(r.get('remaining_special_types',[]))=={'3','2,2'}):
        raise ValueError('degree13 saturated reduction scope')
    if {v.get('type')for v in r.get('rows',[])}!={'2,2','3,2','2,2,2'}or any(v.get('unresolved')!=0 for v in r['rows']):
        raise ValueError('degree13 saturated coverage incomplete')
    items=r.get('bindings',[]);required={'DEGREE_THIRTEEN_SATURATED.md','d13_saturated.rs','audit_d13_saturated.py'}
    if len(items)<12 or not required<={Path(v['path']).name for v in items}:raise ValueError('degree13 coverage bindings missing')
    for item in items:
        source=Path(item['path']);source=source if source.is_absolute()else PROJECT/source
        if sha(source)!=item['sha256']:raise ValueError('degree13 proof binding changed')
    return dict(path=str(path.resolve()),sha256=sha(path),**r)


def admit_complete_degree_thirteen(path):
    path=Path(path);r=json.loads(path.read_text())
    if not(r.get('passed')is True and r.get('claim')=='C230'and r.get('source_genus')==3 and r.get('degree')==13
        and r.get('target_degree_lower_bound')==14 and r.get('scope')=='Complete identity-free transitive FULL six-move all-smooth exclusion, all groups'):
        raise ValueError('complete degree13 exclusion scope')
    for key,necessary in(('single_three',2808),('double_transposition',6435312)):
        if r.get(key,{}).get('unresolved')!=0 or r[key].get('necessary_e')!=necessary:
            raise ValueError('complete degree13 coverage incomplete')
    items=r.get('bindings',[])
    required={'FULL_D13.md','audit_d13_three.rs','audit_d13_double_middle.py','test_d13_middle_audit.py','ENUMERATION.bin'}
    if len(items)<17 or not required<={Path(v['path']).name for v in items}:raise ValueError('complete degree13 bindings missing')
    for item in items:
        source=Path(item['path']);source=source if source.is_absolute()else PROJECT/source
        if sha(source)!=item['sha256']:raise ValueError('complete degree13 proof binding changed')
    return dict(path=str(path.resolve()),sha256=sha(path),**r)


def admit_saturated_passports(path):
    path=Path(path);r=json.loads(path.read_text())
    if not(r.get('schema')==1 and r.get('claim')=='C237'and r.get('passed')is True and
        r.get('theorem')=='saturated-single-cycle-star-prime-degree-passport-exclusions'and
        r.get('production_target_genus')==3 and r.get('conclusion')=='not_all_smooth'and r.get('reviewed_by_parent_and_geometry')is True):
        raise ValueError('saturated passport theorem scope')
    expected=dict(degree='prime',product_one=True,target_transitive=True,
        special_entry='one nontrivial k-cycle, 2<=k<d, all other letters fixed',
        designated_partner='exactly k point cycles, fixed points included',
        ordinary_pair_index_sums='all strictly greater than d-1',
        full_braid_orbit='all six Hurwitz moves modulo simultaneous S_d conjugacy',
        initial_order_or_full_pair_assumed=False,arbitrary_group_allowed=True)
    if r.get('scope')!=expected:raise ValueError('saturated passport hypothesis scope')
    rules=[
        'prime k: every nondesignated ordinary cycle length is coprime to k',
        'k=3: each nondesignated ordinary entry has one common nonzero cycle-length residue modulo3',
        'k=4: each nondesignated ordinary entry has one common odd cycle-length residue modulo4 in EVERY prime degree, by four-tile periodicity, not universal finite-field invertibility',
        'ell|k prime and <ell> union {0} contains k consecutive residues modulo d: each nondesignated ordinary entry has one common cycle-length residue modulo ell^v_ell(k)',
        'for a congruence modulus q, floor(d/k)<=q implies a genus-zero contradiction']
    if r.get('rules')!=rules:raise ValueError('saturated passport arithmetic rules differ')
    if r.get('checks',{}).get('parent_unresolved')!=0 or r['checks'].get('retained_positive_genus_three_components')!=1793:
        raise ValueError('saturated passport parent controls incomplete')
    required={'MATHEMATICS.md','check_parent.py','PARENT_CHECKS.json','ADVERSARIAL_REVIEW.md','audit_cyclic_reductions.py',
        'CYCLIC_COUNTEREXAMPLE_CHECKS.json','PRIME_POWER_FOUR_CYCLE.md','audit_prime_power_tiles.py','PRIME_POWER_COUNTEREXAMPLE_CHECKS.json'}
    items=r.get('bindings',[])
    if len(items)!=9 or {Path(v['path']).name for v in items}!=required:raise ValueError('saturated passport proof bindings missing')
    for item in items:
        source=Path(item['path']);source=source if source.is_absolute()else PROJECT/source
        if sha(source)!=item['sha256']:raise ValueError('saturated passport proof binding changed')
    return dict(path=str(path.resolve()),sha256=sha(path),**r)


def scheduled_task(slot,epoch,plan):
    n,lane=task(slot,epoch,plan.get('search_policy','default'))
    if lane=='involution-boundary'and plan.get('boundary_cursor',0)==BOUNDARY_TOTAL:n,lane=12,'collision'
    if n==12 and plan.get('degree12_reduction'):n,lane=13,'collision'
    if n in(12,13)and plan.get('complete_degree13_reduction'):n,lane=14,'collision'
    return n,lane


def task(slot, epoch, policy='default'):
    """Target-first strata, not a census or a mathematical degree ceiling.

    63 is this build's checked u64-mask limit, not an existence bound. The
    frozen old controller retains its own old schedule until validated handoff.
    """
    if policy not in ('default','overnight-refactor'):
        raise ValueError('unknown search policy')
    if policy=='overnight-refactor' and slot==2 and epoch%2==0:
        return [18,20,21,24,25][(epoch//2+2)%5],'refactor'
    if policy=='overnight-refactor' and slot==8 and epoch%3!=2:
        # New proposal law, not a scope exclusion: keep the original general
        # fibre action on every third job and all other streams unchanged.
        return [18,20,21,24,25][(epoch//3)%5],'refactor'
    if slot == 0:
        degrees=[21,24,27,28,30,32,33,36,39,40,42,44,45,48,51,52,54,56,57,60,63]
        return degrees[epoch%len(degrees)],'partial-fiber' if (epoch+epoch//len(degrees))%3 else 'full-fiber'
    if slot == 1:
        # Non-full quotient beta, including small inner degree3/4. Avoids the
        # old single-local unramified isotrivial construction at frame creation.
        degrees=[12,15,16,18,20,21,24,27,28,30,32,33,36,40,42,45,48,52,54,56,60,63]
        return degrees[epoch%len(degrees)],'general-fiber'
    if slot == 8:
        degrees=[12,16,18,20,24,27,28,32,36,40,42,48,52,56,60,63]
        return degrees[(epoch+7)%len(degrees)],'general-fiber'
    if slot == 9:
        # C209 excludes the ENTIRE old13 triangle budget by fixed-point twists.
        # This open12 star budget has a two-cycle first product, so that
        # full-cycle obstruction does not decide it. No blocks imposed.
        return (16,'involution-boundary') if epoch%2 else (12,'collision')
    if slot == 5:
        # A small exploration/control lane, not another known-family farm.
        lane=epoch%4
        if lane==0:
            degrees=[15,18,20,21,24,25]
            return degrees[(epoch//4)%len(degrees)],'mutation'
        if lane==1:
            degrees=[18,24,27,30,36,42,48,54,60]
            return degrees[(epoch//4)%len(degrees)],'elliptic'
        if lane==2:
            degrees=[12,15,16,18,20,21,24,28,30,36,42,48,54,60]
            return degrees[(epoch//4)%len(degrees)],'rational'
        return list(range(12,64))[(epoch//4)%52],'legacy-types'
    if slot == 4:
        # Four nonidentity quotient entries; constrained return budgets. The
        # Full decorated quotient audit rules out the marked r=4 pure-A/B
        # construction below52, including disjoint transpositions. Other r
        # stay eligible; budget feasibility is not existence.
        degrees=[33,35,36,40,42,44,45,48,49,50,52,54,55,56,60,63]
        return degrees[epoch%len(degrees)],'four-quotient'
    # Degrees9/10 already have complete positive-genus censuses. C209's scalar
    # proof + exact slices exclude ALL source-genus3 in11; other genera11 are
    # not excluded. These constructors target genus3 only.
    # No blocks imposed. Interleave high and low degrees so each short epoch
    # explores both, rather than waiting through a whole monotone sweep.
    degrees=[x for pair in zip(range(12,38),range(63,37,-1)) for x in pair]
    return degrees[(epoch + (slot - 2) * 17) % len(degrees)], 'collision'


def job_seed(seed, slot, epoch):
    """Old six streams unchanged; extra slots use a disjoint counter namespace.

    Within the admitted epoch range, old counters stay below2^43 and new
    counters above2^62. Odd multiplication modulo2^64 is injective. Never
    change the stride from6 to10: doing that would lose the old streams.
    """
    if not 0 <= slot < 10 or not 0 <= epoch < 2**40:
        raise ValueError('seed slot/epoch outside checked namespace')
    counter = 1+slot+6*epoch if slot<6 else 2**62+1+(slot-6)+4*epoch
    return (seed+0x9E3779B97F4A7C15*counter) % 2**64


def expanded_epochs(epochs, workers):
    if not 1 <= len(epochs) <= workers <= 10 or any(type(e) is not int or e<0 or e>=2**40 for e in epochs):
        raise ValueError('worker-lane cursor mapping cannot shrink or discard lanes')
    return list(epochs)+[0]*(workers-len(epochs))


def available_slots(streams, busy, workers, cursor, weights=None, issued=None):
    """Rotate saved streams fairly, with a hard shared producer/auditor limit.

    Stream identities and seed cursors are independent of CPU concurrency.
    Never queue work beyond the limit: queued jobs would mask later STOPs.
    """
    if not 1 <= workers <= streams <= 10 or not 0 <= cursor < streams or any(s not in range(streams) for s in busy):
        raise ValueError('invalid stream/concurrency allocation')
    if weights is not None:
        if len(weights)!=streams or any(type(w)is not int or not 1<=w<=16 for w in weights):raise ValueError('invalid positive stream weights')
        if issued is None or len(issued)!=streams or any(type(k)is not int or k<0 for k in issued):raise ValueError('invalid issued-job counts')
        ready=[s for s in range(streams) if s not in busy]
        ready.sort(key=lambda s:(Fraction(issued[s]+1,weights[s]),(s-cursor)%streams))
        chosen=ready[:max(0,workers-len(busy))]
        return chosen,(chosen[-1]+1)%streams if chosen else cursor
    room=max(0,workers-len(busy)); chosen=[]
    for step in range(streams):
        if len(chosen)==room:break
        slot=(cursor+step)%streams
        if slot not in busy:chosen.append(slot)
    return chosen, (chosen[-1]+1)%streams if chosen else cursor


def learned_words(rows, limit=64):
    counts = {}
    for word in rows:
        if word and len(word) <= 24 and all(0 <= x < 6 for x in word):
            key = tuple(word); counts[key] = counts.get(key, 0) + 1
    return sorted(counts, key=lambda w: (-counts[w], len(w), w))[:limit]


def next_word_bank(history,inherited):
    result=learned_words(history)
    result.extend(tuple(w) for w in inherited if tuple(w) not in result)
    return result[:64]


def mutation_seed_text(catalog):
    rows=[]
    for item in catalog['components']:
        if item['genus']!=3:continue
        partitions=[]
        for b in item['block_quotients']:
            labels=[-1]*item['degree']
            for i,block in enumerate(b['blocks']):
                for x in block:labels[x-1]=i
            if -1 in labels:raise ValueError('incomplete source block partition')
            partitions.append(','.join(map(str,labels)))
        rows.append(str(item['degree'])+'\t['+', '.join(item['tuple'])+']\t'+';'.join(partitions)+'\n')
    return ''.join(rows)


def profile_component(folder, expected_id):
    # Literal finite arithmetic only: no public relabeling and no automatic
    # promotion of the draft Claim 8 to a verified geometric theorem.
    from export_research_catalog import first_row, block_details, tuple_details
    from classify_smooth_isotrivial import classify
    from verify_smooth_geometry import parse_tuple
    from genus_three_fiber_geometry import full_fiber_witness,partial_fiber_witness
    summary = first_row(folder / 'SUMMARY.tsv'); n = int(summary['degree'])
    t = parse_tuple(summary['seed'], n); lattice = block_details(t)
    witnesses = []
    for system in lattice['proper_systems']:
        eps = system['ramification_excess']; h = system['genus']
        for entry in system['identity_entries']:
            k = entry - 1
            full = full_fiber_witness(t, system['partition'], k)
            if full is not None:
                witnesses.append(full)
            partial=partial_fiber_witness(t,system['partition'],k)
            if partial is not None:
                witnesses.append(partial)
            if eps[k] > 0 and (h >= 2 or h == 1 and sum(eps) > eps[k]):
                witnesses.append(dict(special_entry=entry, quotient_genus=h,
                    relative_excess=eps, partition=system['partition'],
                    criterion='fixed-hyperbolic' if h >= 2 else 'draft-claim8-elliptic-anchor'))
    detail=tuple_details(t)
    original=classify(dict(degree=n,tuple=[p['cycle_notation'] for p in detail['permutations']],
        block_quotients=[dict(block_size=b['block_size'],block_count=b['blocks'],blocks=b['partition'],
            genus=b['genus'],tuple=[q['cycle_notation'] for q in b['quotient_permutations']])
            for b in lattice['proper_systems']]))
    if original['isotrivial'] is True and witnesses:raise ValueError('conflicting geometric criteria: review proofs')
    return dict(id=expected_id, degree=n, genus=int(summary['genus']),
        states=int(summary['states']), representative=detail,original_sufficient_criterion=original,
        block_lattice=lattice, prospective_nonisotrivial_witnesses=witnesses,
        geometry='requires proof review; full smooth closure alone is insufficient')


def bounded_profile_component(folder,expected_id,deadline_epoch=None):
    """Potentially costly exact block/group work must respect the SAME deadline.

    An interrupted profile leaves the already verified closure and raw tables
    intact; it never turns that new component into a negative or an iso label.
    """
    if deadline_epoch is None:return profile_component(folder,expected_id)
    command=[sys.executable,str(Path(__file__).resolve()),'--profile',str(folder),'--id',expected_id]
    value=json.loads(subprocess.check_output(command,text=True,timeout=remaining_timeout(deadline_epoch)))
    if value.get('id')!=expected_id:raise ValueError('profile child identity mismatch')
    return value


def inherited_known(old):
    """Keep identities through zero-new-result successors, not only one edge."""
    old = Path(old)
    if (old/'BINDINGS.json').exists():
        expected=json.loads((old/'BINDINGS.json').read_text()).get('KNOWN_IDS.json')
        if expected is None or sha(old/'KNOWN_IDS.json')!=expected:
            raise ValueError('predecessor known-ID binding changed')
    ids = set(json.loads((old/'KNOWN_IDS.json').read_text()))
    previous = json.loads((old/'NEW_COMPONENTS.json').read_text()) if (old/'NEW_COMPONENTS.json').exists() else []
    for item in previous:
        cert = Path(item['raw_folder'])
        if not (cert.parent/'VERIFIED').exists() or sha(cert/'STATES.tsv') != item['table_sha256'] or sha(cert.parent/'CHECKS.json') != item['checks_sha256']:
            raise ValueError('predecessor component proof changed')
        ids.add(item['id'])
    return sorted(ids)


def run_job(root, plan, slot, epoch, bank):
    n, lane = scheduled_task(slot,epoch,plan); name = f'job_{slot}_{epoch:08}_{n}_{lane}'
    boundary_cursor=plan.get('boundary_cursor',0)
    folder = root / 'jobs' / name
    seed = job_seed(plan['seed'],slot,epoch)
    binary=root/'smooth_genus_three'
    if n<=25 and plan.get('compact_binary_sha256'):binary=root/'smooth_genus_three_compact'
    command = [str(binary), '--d', str(n), '--lane', lane,
        '--seed', str(seed), '--seconds', str(plan['job_seconds']),
        '--closure-seconds', str(plan['closure_seconds']), '--cap', str(plan['cap']),
        '--out', str(folder), '--verifier', plan['verifier'], '--words', str(bank)]
    if lane in ('mutation','refactor'):command.extend(['--mutation-seeds',str(root/'MUTATION_SEEDS.tsv')])
    if lane=='refactor':command.extend(['--refactor-max-index','3'])
    if lane=='refactor' and plan.get('joined_refactor') and (epoch//(2 if slot==2 else 3))%2==1:
        command.append('--refactor-joined')
    if lane=='refactor' and plan.get('adaptive_refactor',0) and (epoch//(2 if slot==2 else 3))%3==0:
        # Preserve two static-control rounds per adaptive round. This trades
        # some depth for new right-pair templates, not a universal speed gain.
        command.extend(['--adaptive-refactor',str(plan['adaptive_refactor'])])
    if n==14 and plan.get('complete_degree13_reduction')and slot==9 and lane=='collision':
        command.extend(['--index-profile',['2,8,11,11','2,9,10,11','2,10,10,10','3,9,10,10'][epoch%4]])
    elif n==13 and plan.get('degree13_reduction')and lane in('collision','legacy-types','joint-genus3'):
        command.extend(['--index-profile','2,9,9,10'])
    elif slot==9 and lane=='collision':
        profile=['2,8,10,10','2,9,9,10','3,9,9,9'][(epoch//2)%3]if plan.get('degree12_reduction')else('2,8,9,9'if plan.get('transposition_reduction')else'1,9,9,9')
        command.extend(['--index-profile',profile])
    if plan.get('transposition_reduction'):command.append('--exclude-transpositions')
    if plan.get('saturated_passport_reduction'):command.extend(['--saturated-passport-cuts','--saturated-stage',plan.get('saturated_passport_stage','late')])
    if plan.get('frame_retirement_reduction'):command.append('--retire-proven-frames')
    if lane=='involution-boundary':command.extend(['--cursor',str(boundary_cursor)])
    if lane!='mutation':command.append('--literal-trie')
    if lane in ['partial-fiber','full-fiber']:command.append('--skip-covered-fiber')
    if plan.get('trace_early',0):command.extend(['--trace-early',str(plan['trace_early'])])
    if lane=='four-quotient' and plan.get('trace_trie',0):command.extend(['--trace-trie',str(plan['trace_trie'])])
    timeout = None
    if plan.get('deadline_epoch') is not None:
        timeout = plan['deadline_epoch'] - time.time()
        if timeout <= 0:raise TimeoutError('absolute research deadline reached before producer launch')
    with (root / 'logs' / (name + '.log')).open('w') as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=timeout)
    checks = audit_guided_native(folder, root / 'smooth_negative_audit', native_components=True,deadline_epoch=plan.get('deadline_epoch'))
    if (folder/'ADAPTIVE.tsv').exists():
        from verify_genus_three_adaptive import verify_adaptive
        checks['adaptive_template_provenance']=verify_adaptive(folder,deadline_epoch=plan.get('deadline_epoch'))
    if plan.get('transposition_reduction'):
        from verify_genus_three_transpositions import verify_theorem_exclusions
        checks['target_theorem_exclusions']=verify_theorem_exclusions(folder)
        checks['target_theorem_exclusions']['proof_receipt_sha256']=sha(root/'TARGET_THEOREM.json')
    if plan.get('saturated_passport_reduction'):
        from verify_genus_three_saturated import verify_saturated_exclusions
        checks['saturated_passport_exclusions']=verify_saturated_exclusions(folder)
        checks['saturated_passport_exclusions']['proof_receipt_sha256']=sha(root/'SATURATED_PASSPORT_THEOREM.json')
    if plan.get('frame_retirement_reduction'):
        from genus_three_frame_retirement import verify as verify_retired_frames
        checks['proposal_frame_retirement']=verify_retired_frames(folder)
        checks['proposal_frame_retirement']['proof_receipt_sha256']=sha(root/'FRAME_RETIREMENT_THEOREM.json')
    boundary=None
    if lane=='involution-boundary':
        from verify_genus_three_boundary import audit_boundary
        boundary=audit_boundary(folder,expected_start=boundary_cursor)
        checks['boundary_coverage']=boundary
    if (folder/'ISOTRIVIAL.tsv').exists():
        from verify_genus_three_exclusions import audit_exclusions
        with (folder/'STAGES.tsv').open() as stream:stages=next(csv.DictReader(stream,delimiter='\t'))
        checks['isotrivial_exclusions']=audit_exclusions(folder/'ISOTRIVIAL.tsv',int(stages['isotrivial_outside_target']))
    checks['word_bank'] = dict(path=str(bank), sha256=sha(bank))
    checks['scope'] = 'finite two-path degree16 interval; not a degree census' if boundary else 'sampled identity-free genus-three proposals; not a degree census'
    with (folder / 'COUNTS.tsv').open() as stream:
        counts = next(csv.DictReader(stream, delimiter='\t'))
    words = []
    with (folder / 'RESOLUTIONS.tsv').open() as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            if row['status'] == 'nonsmooth':
                words.append(list(map(int, row['word'].split(','))))
    write(folder / 'CHECKS.json', checks)
    (folder / 'VERIFIED').write_text('All persisted rejections and exceptional six-move tables audited.\n')
    components = []
    for item in checks['closed_identities']:
        cid = identity(n, item['key']); cert = folder / item['folder']
        components.append(dict(id=cid, degree=n, genus=3, states=int(item['summary']['states']),
            lane=lane, raw_folder=str(cert), component_identity_key=item['key'],
            table_sha256=sha(cert / 'STATES.tsv'), summary_sha256=sha(cert / 'SUMMARY.tsv'),
            checks_sha256=sha(folder / 'CHECKS.json'), representative=item['summary']['seed']))
        if cid not in set(json.loads((root/'KNOWN_IDS.json').read_text())):
            profile=bounded_profile_component(cert,cid,plan.get('deadline_epoch'));write(cert/'GEOMETRY_REVIEW.json',profile)
            components[-1]['original_isotrivial']=profile['original_sufficient_criterion']['isotrivial']
            components[-1]['geometry_review_sha256']=sha(cert/'GEOMETRY_REVIEW.json')
            components[-1]['prospective_witness_count']=len(profile['prospective_nonisotrivial_witnesses'])
    return dict(slot=slot, epoch=epoch, folder=str(folder), counts=counts,
                components=components, words=words, boundary=boundary)


def controller(root):
    plan = json.loads((root / 'PLAN.json').read_text()); began = time.monotonic()
    status = lambda name, **kw: write(root / 'STATE.json', dict(status=name, pid=os.getpid(), updated=now(), **kw))
    status('queued', active_jobs=[])
    lock = Path(plan['campaign_lock']).open('a'); fcntl.flock(lock, fcntl.LOCK_EX)
    def stop(signum, frame):
        (root / 'STOP').write_text(f'Graceful signal {signum}; finish bounded jobs and audits.\n')
    signal.signal(signal.SIGTERM, stop); signal.signal(signal.SIGINT, stop)
    expected = json.loads((root / 'BINDINGS.json').read_text())
    for name, digest in expected.items():
        if sha(root / name) != digest:
            raise ValueError('Frozen input binding changed: ' + name)
    known = set(json.loads((root / 'KNOWN_IDS.json').read_text()))
    initial=json.loads((root/'INHERITED_CHECKPOINT.json').read_text()) if (root/'INHERITED_CHECKPOINT.json').exists() else {}
    streams=plan.get('streams',plan['workers'])
    epochs = initial.get('next_epochs',[0]*streams); active = {}; discoveries = {}; completed = 0
    if len(epochs)!=streams:raise ValueError('saved stream count differs')
    cursor=initial.get('scheduler_cursor',0)
    boundary_cursor=initial.get('boundary_cursor',0)
    if type(boundary_cursor)is not int or not 0<=boundary_cursor<=BOUNDARY_TOTAL:raise ValueError('boundary cursor outside finite scope')
    weights=plan.get('stream_weights')
    origin=initial.get('allocation_origin',epochs.copy()) if initial.get('allocation_weights')==weights else epochs.copy()
    if len(origin)!=streams or any(o>e for o,e in zip(origin,epochs)):raise ValueError('allocation checkpoint differs')
    totals = dict(constructed=0, scope=0, transitive=0, closed_occurrences=0, unresolved=0)
    bank = root / 'banks' / 'words_00000000.txt'
    # A zero-new successor still needs the inherited literal rejection bank;
    # otherwise its first empty journal silently discards every trained word.
    inherited_words=[list(map(int,line.split(','))) for line in bank.read_text().splitlines() if line.strip()]
    word_history = inherited_words.copy()
    stop_reason = None
    with ThreadPoolExecutor(max_workers=plan['workers']) as pool:
        while True:
            if (root / 'STOP').exists(): stop_reason = stop_reason or 'requested'
            if plan['hours'] and time.monotonic() - began >= plan['hours'] * 3600: stop_reason = stop_reason or 'time_budget'
            if deadline_dispatch_closed(plan):stop_reason=stop_reason or 'absolute_deadline_drain'
            if not storage_available(root): stop_reason = stop_reason or 'disk_guard_10GiB'
            if not stop_reason:
                busy = {v[0] for v in active.values()}
                ready,cursor=available_slots(streams,busy,plan['workers'],cursor,weights,[e-o for e,o in zip(epochs,origin)])
                for slot in ready:
                    e = epochs[slot]; epochs[slot] += 1
                    active[pool.submit(run_job, root, dict(plan,boundary_cursor=boundary_cursor), slot, e, bank)] = (slot, e)
            status('stopping' if stop_reason else 'running', active_jobs=list(active.values()),
                completed_jobs=completed, epochs=epochs, counts=totals,
                new_closed_components=len(discoveries), stop_reason=stop_reason,boundary_cursor=boundary_cursor,boundary_total=BOUNDARY_TOTAL)
            if not active: break
            done, _ = wait(active, timeout=2, return_when=FIRST_COMPLETED)
            for future in done:
                slot, e = active.pop(future)
                try: result = future.result()
                except Exception as exc:
                    write(root / f'FAILURE_{slot}_{e}.json', dict(error=repr(exc), updated=now()))
                    stop_reason = 'audit_or_producer_failure'; continue
                completed += 1; c = result['counts']
                if result.get('boundary'):
                    b=result['boundary']
                    if not b.get('passed') or b['start']!=boundary_cursor or not boundary_cursor<=b['next']<=BOUNDARY_TOTAL:
                        stop_reason='boundary_coverage_failure';continue
                    boundary_cursor=b['next']
                    write(root/'BOUNDARY_PROGRESS.json',dict(next=boundary_cursor,total=BOUNDARY_TOTAL,
                        folder=result['folder'],checks_sha256=sha(Path(result['folder'])/'CHECKS.json'),updated=now()))
                for target, key in [('constructed','constructed'),('scope','scope'),('transitive','transitive'),
                                    ('closed_occurrences','closed'),('unresolved','unresolved')]:
                    totals[target] += int(c[key])
                for item in result['components']:
                    if item['id'] not in known and item['id'] not in discoveries:
                        discoveries[item['id']] = item
                        write(root / 'NEW_COMPONENTS.json', list(discoveries.values()))
                        # Continue through established isotrivial constructions.
                        # An unknown is not called non-isotrivial; stop for proof.
                        if item.get('original_isotrivial') is not True:
                            stop_reason = 'new_genus3_component_requires_geometry_review'
                word_history.extend(result['words'])
                if len(word_history) > 10000: word_history = word_history[-10000:]
                bank = root / 'banks' / f'words_{completed:08}.txt'
                prioritized=next_word_bank(word_history,inherited_words)
                bank.write_text(''.join(','.join(map(str,w))+'\n' for w in prioritized[:64]))
                with (root / 'COMPLETED.jsonl').open('a') as stream:
                    stream.write(json.dumps({k:v for k,v in result.items() if k != 'words'})+'\n')
        final = dict(status='stopped' if stop_reason == 'requested' else 'attention' if stop_reason not in ['time_budget','absolute_deadline_drain'] else 'complete',
            reason=stop_reason, completed_jobs=completed, epochs=epochs, counts=totals,
            new_closed_components=len(discoveries), boundary_cursor=boundary_cursor,boundary_total=BOUNDARY_TOTAL,updated=now(), seconds=time.monotonic()-began)
        write(root/'CHECKPOINT.json',dict(next_epochs=epochs,scheduler_cursor=cursor,allocation_origin=origin,allocation_weights=weights,word_bank=str(bank),word_bank_sha256=sha(bank),
            counts=totals,completed_jobs=completed,new_components=list(discoveries),boundary_cursor=boundary_cursor,
            deadline_epoch=plan.get('deadline_epoch'),
            pending='Every capped/timeout seed remains in the immutable job journals; failed jobs require recovery before reuse.'))
        write(root / 'FINISHED.json', final); status(final['status'], active_jobs=[], **{k:v for k,v in final.items() if k not in ['status','updated']})
    lock.close()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run', action='store_true'); p.add_argument('--worker', type=Path)
    p.add_argument('--profile', type=Path); p.add_argument('--id', default='')
    p.add_argument('--binary', type=Path); p.add_argument('--auditor', type=Path)
    p.add_argument('--compact-binary',type=Path,help='validated low-degree specialization; same target/certificates')
    p.add_argument('--source', type=Path); p.add_argument('--checks', type=Path); p.add_argument('--catalog', type=Path)
    p.add_argument('--predecessor',type=Path,help='finalized stopped run: retain epochs/known IDs/literal word library, never restart it')
    p.add_argument('--workers', type=int, default=6); p.add_argument('--hours', type=float, default=0)
    p.add_argument('--deadline',help='absolute ISO8601 deadline with UTC offset; inherited successors cannot extend it')
    p.add_argument('--clear-expired-deadline',metavar='HUMAN_REASON',help='NEW explicit human authorization to resume a completed expired mission; reason is retained, predecessor is immutable')
    p.add_argument('--trace-early',type=int,default=0,help='bounded independently audited early negative witnesses per job, <=10000')
    p.add_argument('--trace-trie',type=int,default=0,help='bounded exact early trie witnesses per four-quotient job')
    p.add_argument('--search-policy',choices=['default','overnight-refactor'],default='default')
    p.add_argument('--adaptive-refactor',type=int,default=0,help='bounded negative-prefix proposal-template cache, never positive/skip admission')
    p.add_argument('--joined-refactor',action='store_true',help='alternate generic cycle-joining and old disjoint-involution refactor rounds')
    p.add_argument('--transposition-proof',type=Path,help='sealed reviewed no-positive-genus-all-smooth-transposition theorem receipt')
    p.add_argument('--degree12-proof',type=Path,help='sealed complete C226/C227 degree-twelve exclusion; requires transposition theorem')
    p.add_argument('--degree13-proof',type=Path,help='sealed C229 saturated reductions; keep index budget2,9,9,10 and both special types')
    p.add_argument('--complete-degree13-proof',type=Path,help='sealed C230 complete all-group degree13 negative; rotate four still-open degree14 star budgets')
    p.add_argument('--saturated-passport-proof',type=Path,help='sealed C237 single-cycle/saturated/STAR prime-degree theorem; separate arithmetic target exclusions')
    p.add_argument('--saturated-passport-stage',choices=['after-short','late'],default='late',help='never enable unmeasured expensive early placement in production')
    p.add_argument('--frame-retirement-proof',type=Path,help='C251 recursively bound NONISOTRIVIAL-goal recipe exclusions BEFORE third entries are drawn')
    p.add_argument('--streams',type=int,help='saved seed streams; defaults to retaining all predecessor streams while limiting concurrent workers')
    p.add_argument('--stream-weights',help='comma-separated positive job shares; heuristics only, never exclusions')
    p.add_argument('--job-seconds', type=int, default=45); p.add_argument('--closure-seconds', type=int, default=60)
    p.add_argument('--cap', type=int, default=1_000_000); p.add_argument('--seed', type=int, default=20261021)
    p.add_argument('--verifier', choices=['macro','ordinary'], default='macro')
    a = p.parse_args()
    if a.profile:
        print(json.dumps(profile_component(a.profile, a.id), indent=2)); return
    if a.worker:
        try: controller(a.worker.resolve())
        except Exception as exc:
            write(a.worker/'CONTROLLER_FAILURE.json',dict(error=repr(exc),updated=now()))
            write(a.worker/'STATE.json',dict(status='failed',pid=os.getpid(),error=repr(exc),updated=now()))
            raise
        return
    if not 1 <= a.workers <= min(10,os.cpu_count() or 1) or a.hours < 0 or min(a.job_seconds,a.closure_seconds,a.cap) <= 0: p.error('invalid bounded pool configuration')
    if not 0<=a.trace_early<=10000 or not 0<=a.trace_trie<=10000:p.error('early trace bounds must be0..10000')
    if not 0<=a.adaptive_refactor<=1024:p.error('adaptive template bound must be0..1024')
    if a.degree12_proof and not a.transposition_proof:p.error('degree12 proof depends on the admitted transposition theorem')
    if a.degree13_proof and not a.degree12_proof:p.error('degree13 schedule requires degree12 and transposition receipts')
    if a.complete_degree13_proof and not a.degree13_proof:p.error('complete degree13 requires all prerequisite receipts')
    try:
        deadline=absolute_deadline(a.deadline)
        deadline,deadline_reset=successor_deadline(deadline,a.predecessor,a.clear_expired_deadline)
    except ValueError as exc:p.error(str(exc))
    prior_streams=len(json.loads((a.predecessor/'CHECKPOINT.json').read_text())['next_epochs']) if a.predecessor else 0
    streams=a.streams if a.streams is not None else max(a.workers,prior_streams)
    if not a.workers<=streams<=10:p.error('streams must retain the worker budget and be <=10')
    weights=list(map(int,a.stream_weights.split(','))) if a.stream_weights else None
    if weights is not None and (len(weights)!=streams or any(not 1<=w<=16 for w in weights)):p.error('one weight1..16 per stream required')
    proposal = dict(objective='one independently closed, provably non-isotrivial source-genus-three all-smooth component',
        workers=a.workers, streams=streams, stream_weights=weights, gpu=False, gpu_reason='no validated end-to-end GPU discovery gain; resident finite screening prototype is separate',
        hours=a.hours, job_seconds=a.job_seconds, closure_seconds=a.closure_seconds, cap=a.cap,
        deadline_epoch=deadline,deadline_reset_authority=deadline_reset,trace_early=a.trace_early,trace_trie=a.trace_trie,search_policy=a.search_policy,adaptive_refactor=a.adaptive_refactor,joined_refactor=a.joined_refactor,
        joined_schedule='half refactor frame rounds; old involution rounds retained'if a.joined_refactor else'old involution refactor only',
        adaptive_schedule='one of three refactor frame rounds; negative-prefix proposal templates only, not known-positive or skip-bank states',
        seed=a.seed, verifier=a.verifier, quotient='simultaneous S_d conjugacy',
        seed_namespace_version=2,
        schedule_revision=8 if a.joined_refactor else(7 if a.search_policy=='overnight-refactor'else 5),focused_profile=dict(slot=9,alternating=[dict(degree=12,indices=[2,8,9,9]if a.transposition_proof else[1,9,9,9]),
            dict(degree=16,indices=[7,7,7,15],constructor='finite two-path transversal matchings',total=BOUNDARY_TOTAL)],
            exhausted_boundary_fallback='degree12 collision; never repeat completed matching intervals',
            excluded='C209 complete genus3 degree11 and entire13 triangle budget6,6,6,12; higher13 star profiles remain eligible'),
        allocation=[dict(slot=i, degree=task(i,0,a.search_policy)[0],lane=task(i,0,a.search_policy)[1]) for i in range(streams)],
        completion='literal full six-move closure AND separate geometric proof; cap/timeout undecided',
        early_stop='new closed genus-three component outside the established isotrivial criterion triggers geometry review', disk_reserve_bytes=10*2**30,
        project=str(PROJECT), campaign_lock=str(PROJECT / '.run-control/.campaign.lock'))
    if a.transposition_proof:proposal['transposition_reduction']=admit_transposition_reduction(a.transposition_proof)
    if a.degree12_proof:
        proposal['degree12_reduction']=admit_degree_twelve(a.degree12_proof)
        proposal['focused_profile']=dict(slot=9,degree=13,rotating_indices=[[2,8,10,10],[2,9,9,10],[3,9,9,9]],
            shape_exclusions='none added; middle-budget single3-cycle retained',old_boundary_cursor=BOUNDARY_TOTAL)
        proposal['schedule_revision']=9
    if a.degree13_proof:
        proposal['degree13_reduction']=admit_degree_thirteen(a.degree13_proof)
        proposal['focused_profile']=dict(slot=9,degree=13,indices=[2,9,9,10],special_types=['3','2,2'],
            excluded_saturated_budgets=[[2,8,10,10],[3,9,9,9]],old_boundary_cursor=BOUNDARY_TOTAL)
        proposal['schedule_revision']=10
    if a.complete_degree13_proof:
        proposal['complete_degree13_reduction']=admit_complete_degree_thirteen(a.complete_degree13_proof)
        proposal['focused_profile']=dict(slot=9,degree=14,rotating_indices=[[2,8,11,11],[2,9,10,11],[2,10,10,10],[3,9,10,10]],
            excluded_triangle=[6,6,7,13],excluded='complete genus3 degrees9--13; no whole degree14 exclusion',old_boundary_cursor=BOUNDARY_TOTAL)
        proposal['schedule_revision']=11
    if a.saturated_passport_proof:
        proposal['saturated_passport_reduction']=admit_saturated_passports(a.saturated_passport_proof)
        proposal['saturated_passport_stage']=a.saturated_passport_stage
    if a.frame_retirement_proof:
        from genus_three_frame_retirement import admit as admit_frame_retirement
        proposal['frame_retirement_reduction']=admit_frame_retirement(a.frame_retirement_proof,PROJECT)
    if not a.run: print(json.dumps(proposal,indent=2)); return
    if deadline_dispatch_closed(proposal):raise ValueError('absolute deadline leaves no safe job/audit interval')
    if not all([a.binary,a.auditor,a.source,a.checks,a.catalog]): p.error('--run needs binary/auditor/source/checks/catalog')
    checks = json.loads(a.checks.read_text())
    if not checks.get('passed') or sha(a.binary)!=checks['binary_sha256'] or sha(a.auditor)!=checks['audit_binary_sha256']:
        raise ValueError('completed validation/binary binding required')
    if a.saturated_passport_proof and not(checks.get('features',{}).get('saturated_passport_filters')is True and
        'src/smooth_saturated_passport.rs'in checks.get('compiled_source_files',{})):
        raise ValueError('saturated passport filters require their completed source-bound validation')
    if a.frame_retirement_proof and not(checks.get('features',{}).get('proposal_frame_retirement')is True and
        'src/genus_three_frame_retirement.rs'in checks.get('compiled_source_files',{})):
        raise ValueError('proposal FRAME retirement requires its completed source-bound validation')
    if a.compact_binary:
        if sha(a.compact_binary)!=checks.get('compact_binary_sha256'):raise ValueError('compact binary validation binding required')
        proposal['compact_binary_sha256']=sha(a.compact_binary)
    proposal['elementary_reductions']=admit_elementary_reductions(PROJECT/'evidence/genus_three_twists_20261022/EXCLUSIONS.json')
    files = compiled_source_files(a.source,checks['compiled_source_files'])
    inherited={}; extra_known=[]; initial_words=''
    if a.predecessor:
        old=a.predecessor.resolve()
        finalized_predecessor(old)
        checkpoint=json.loads((old/'CHECKPOINT.json').read_text())
        next_epochs=expanded_epochs(checkpoint['next_epochs'],streams)
        if sha(checkpoint['word_bank'])!=checkpoint['word_bank_sha256']:raise ValueError('inherited word bank changed')
        initial_words=Path(checkpoint['word_bank']).read_text()
        extra_known=inherited_known(old)
        origin=expanded_epochs(checkpoint.get('allocation_origin',checkpoint['next_epochs']),streams) if checkpoint.get('allocation_weights')==weights else next_epochs
        inherited=dict(next_epochs=next_epochs,scheduler_cursor=checkpoint.get('scheduler_cursor',0),
            boundary_cursor=checkpoint.get('boundary_cursor',0),
            allocation_origin=origin,
            allocation_weights=weights,predecessor=str(old),
            prior_next_epochs=checkpoint['next_epochs'],new_slots=list(range(len(checkpoint['next_epochs']),streams)),
            checkpoint_sha256=sha(old/'CHECKPOINT.json'),finished_sha256=sha(old/'FINISHED.json'),
            known_ids_sha256=sha(old/'KNOWN_IDS.json'),
            inherited_known_ids=len(extra_known),prior_counts=checkpoint['counts'])
    control = PROJECT / '.run-control'; control.mkdir(exist_ok=True)
    lock = (control / '.campaign.lock').open('a'); fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    if not storage_available('data'): raise RuntimeError('T7 mount/10GiB guard')
    root = Path('data') / ('genus3_' + now().replace(':','').replace('-','').replace('.','_'))
    root.mkdir(); (root/'jobs').mkdir(); (root/'logs').mkdir(); (root/'banks').mkdir()
    (root/'runtime/scripts').mkdir(parents=True)
    for source in (PROJECT/'scripts').rglob('*.py'): shutil.copyfile(source,root/'runtime/scripts'/source.name)
    for rel,path in files.items():
        dest=root/'source'/rel; dest.parent.mkdir(parents=True,exist_ok=True); shutil.copyfile(path,dest)
    shutil.copyfile(a.binary,root/'smooth_genus_three'); (root/'smooth_genus_three').chmod(0o700)
    if a.compact_binary:
        shutil.copyfile(a.compact_binary,root/'smooth_genus_three_compact');(root/'smooth_genus_three_compact').chmod(0o700)
    shutil.copyfile(a.auditor,root/'smooth_negative_audit'); (root/'smooth_negative_audit').chmod(0o700)
    shutil.copyfile(a.checks,root/'VALIDATION.json')
    if a.transposition_proof:
        theorem={k:v for k,v in proposal['transposition_reduction'].items()if k not in('path','sha256')}
        (root/'proofs').mkdir()
        for key,name in(('mathematics','TRANSPOSITION_GK_REDUCTION.md'),('primary_source','Gabai_Kazez_1987.pdf')):
            dest=root/'proofs'/name;shutil.copyfile(theorem[key]['path'],dest)
            theorem[key]=dict(theorem[key],path=str(dest))
        write(root/'TARGET_THEOREM.json',theorem)
        proposal['transposition_reduction']=admit_transposition_reduction(root/'TARGET_THEOREM.json')
    if a.degree12_proof:
        theorem={k:v for k,v in proposal['degree12_reduction'].items()if k not in('path','sha256')}
        destdir=root/'proofs/degree12';destdir.mkdir()
        for i,item in enumerate(theorem['bindings']):
            src=Path(item['path']);src=src if src.is_absolute()else PROJECT/src
            dest=destdir/f'{i:03}'/src.name;dest.parent.mkdir();shutil.copyfile(src,dest);item['path']=str(dest)
        write(root/'DEGREE12_EXCLUSION.json',theorem)
        proposal['degree12_reduction']=admit_degree_twelve(root/'DEGREE12_EXCLUSION.json')
    if a.degree13_proof:
        theorem={k:v for k,v in proposal['degree13_reduction'].items()if k not in('path','sha256')}
        destdir=root/'proofs/degree13';destdir.mkdir()
        for i,item in enumerate(theorem['bindings']):
            src=Path(item['path']);src=src if src.is_absolute()else PROJECT/src
            dest=destdir/f'{i:03}'/src.name;dest.parent.mkdir();shutil.copyfile(src,dest);item['path']=str(dest)
        write(root/'DEGREE13_SATURATED.json',theorem)
        proposal['degree13_reduction']=admit_degree_thirteen(root/'DEGREE13_SATURATED.json')
    if a.complete_degree13_proof:
        theorem={k:v for k,v in proposal['complete_degree13_reduction'].items()if k not in('path','sha256')}
        destdir=root/'proofs/degree13_complete';destdir.mkdir()
        for i,item in enumerate(theorem['bindings']):
            src=Path(item['path']);src=src if src.is_absolute()else PROJECT/src
            dest=destdir/f'{i:03}'/src.name;dest.parent.mkdir();shutil.copyfile(src,dest);item['path']=str(dest)
        write(root/'DEGREE13_COMPLETE.json',theorem)
        proposal['complete_degree13_reduction']=admit_complete_degree_thirteen(root/'DEGREE13_COMPLETE.json')
    if a.saturated_passport_proof:
        theorem={k:v for k,v in proposal['saturated_passport_reduction'].items()if k not in('path','sha256')}
        destdir=root/'proofs/saturated_passport';destdir.mkdir(parents=True)
        for i,item in enumerate(theorem['bindings']):
            src=Path(item['path']);src=src if src.is_absolute()else PROJECT/src
            dest=destdir/f'{i:03}'/src.name;dest.parent.mkdir();shutil.copyfile(src,dest);item['path']=str(dest)
        write(root/'SATURATED_PASSPORT_THEOREM.json',theorem)
        proposal['saturated_passport_reduction']=admit_saturated_passports(root/'SATURATED_PASSPORT_THEOREM.json')
    if a.frame_retirement_proof:
        from genus_three_frame_retirement import admit as admit_frame_retirement
        theorem={k:v for k,v in proposal['frame_retirement_reduction'].items()if k not in('path','sha256','recursively_checked_files')}
        destdir=root/'proofs/frame_retirement';destdir.mkdir(parents=True)
        for item in theorem['bindings']:
            src=Path(item['path']);src=src if src.is_absolute()else PROJECT/src
            dest=destdir/src.name;shutil.copyfile(src,dest);item['path']=str(dest)
        write(root/'FRAME_RETIREMENT_THEOREM.json',theorem)
        proposal['frame_retirement_reduction']=admit_frame_retirement(root/'FRAME_RETIREMENT_THEOREM.json',PROJECT)
    catalog = json.loads(a.catalog.read_text())
    write(root/'KNOWN_IDS.json',sorted(set([v['id'] for v in catalog['components'] if v['genus']==3]+extra_known)))
    (root/'MUTATION_SEEDS.tsv').write_text(mutation_seed_text(catalog))
    proposal.update(catalog=dict(path=str(a.catalog.resolve()),sha256=sha(a.catalog),
        genus3_components=sum(v['genus']==3 for v in catalog['components'])),
        known_component_ids_total=len(json.loads((root/'KNOWN_IDS.json').read_text())))
    if inherited:write(root/'INHERITED_CHECKPOINT.json',inherited)
    proposal['inheritance']=inherited or None
    allocation_plan=dict(proposal,boundary_cursor=inherited.get('boundary_cursor',0))
    proposal['allocation']=[dict(slot=i,degree=scheduled_task(i,inherited.get('next_epochs',[0]*streams)[i],allocation_plan)[0],
        lane=scheduled_task(i,inherited.get('next_epochs',[0]*streams)[i],allocation_plan)[1]) for i in range(streams)]
    write(root/'PLAN.json',proposal); (root/'banks/words_00000000.txt').write_text(initial_words)
    bindings={str(path.relative_to(root)):sha(path) for parent in [root/'source',root/'runtime',root/'proofs'] for path in parent.rglob('*') if path.is_file()}
    for name in ['PLAN.json','KNOWN_IDS.json','MUTATION_SEEDS.tsv','VALIDATION.json','smooth_genus_three','smooth_negative_audit','banks/words_00000000.txt']:
        bindings[name]=sha(root/name)
    if a.compact_binary:bindings['smooth_genus_three_compact']=sha(root/'smooth_genus_three_compact')
    if a.transposition_proof:bindings['TARGET_THEOREM.json']=sha(root/'TARGET_THEOREM.json')
    if a.degree12_proof:bindings['DEGREE12_EXCLUSION.json']=sha(root/'DEGREE12_EXCLUSION.json')
    if a.degree13_proof:bindings['DEGREE13_SATURATED.json']=sha(root/'DEGREE13_SATURATED.json')
    if a.complete_degree13_proof:bindings['DEGREE13_COMPLETE.json']=sha(root/'DEGREE13_COMPLETE.json')
    if a.saturated_passport_proof:bindings['SATURATED_PASSPORT_THEOREM.json']=sha(root/'SATURATED_PASSPORT_THEOREM.json')
    if a.frame_retirement_proof:bindings['FRAME_RETIREMENT_THEOREM.json']=sha(root/'FRAME_RETIREMENT_THEOREM.json')
    if inherited:bindings['INHERITED_CHECKPOINT.json']=sha(root/'INHERITED_CHECKPOINT.json')
    write(root/'BINDINGS.json',bindings)
    command=[sys.executable,str(root/'runtime/scripts/run_genus_three.py'),'--worker',str(root)]
    label='local.hurwitz.'+root.name
    pid=spawn_launchd(root,command,label)
    write(root/'LAUNCH.json',dict(pid=pid,label=label,updated=now()))
    write(control/'genus-three-current.json',dict(root=str(root),pid=pid,label=label))
    lock.close(); print(json.dumps(dict(root=str(root),pid=pid,workers=a.workers),indent=2))


if __name__ == '__main__': main()
