"""Independent streaming rejection audit; exceptional components use Python.

The native binary imports no production kernels. Every persisted seed/word is
checked, not sampled. This is not an audit of all unpromoted random proposals.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import json
import hashlib
from pathlib import Path
import subprocess
import threading
import tempfile
import time

from run_smooth_degree9 import rows, sha
from run_smooth_mission import audit_component,values
from odd_order_spin import hurwitz,parse_seed
from run_s18_linked_roots import state_key

_KNOWN = {}
_KNOWN_LOCK=threading.Lock()

def remaining_timeout(deadline_epoch):
    """Optional absolute wall clock budget shared by producer and auditor."""
    if deadline_epoch is None:return None
    remaining=deadline_epoch-time.time()
    if remaining<=0:raise TimeoutError('absolute research deadline reached before independent audit')
    return remaining

def closed_record_digest(records):
    return hashlib.sha256(json.dumps(records,separators=(',',':')).encode()).hexdigest()

def cache_independent_known(path,digest,keys):
    """Derived cache fed only by full independent catalog/closure audits."""
    with _KNOWN_LOCK:
        _KNOWN[(str(path),digest)]=frozenset(keys)
        while len(_KNOWN)>8:del _KNOWN[next(iter(_KNOWN))]

def known_keys(path):
    path=Path(path);digest=sha(path)
    with _KNOWN_LOCK:
        if (str(path),digest) not in _KNOWN:
            keys=set()
            for line in path.read_text().splitlines():
                degree,seed=line.split('\t',1);n=int(degree)
                keys.add((n,state_key(parse_seed(seed,n))))
            _KNOWN[(str(path),digest)]=keys
            # Bounded derived cache; eviction only causes independent parsing
            # again. Keep the freshly built trusted snapshot.
            while len(_KNOWN)>8:del _KNOWN[next(iter(_KNOWN))]
    return _KNOWN[(str(path),digest)]


def prime_known_cache(path,binary):
    """Exact native-derived keys for an already proof-bound positive snapshot.

    Literal checks alone do not establish closure; caller must bind the prior
    complete-closure provenance. No sampled or capped tables are admitted.
    """
    path=Path(path);digest=sha(path)
    with tempfile.TemporaryDirectory(prefix='hurwitz-known-keys-') as tmp:
        export=Path(tmp)/'keys.tsv'
        report=json.loads(subprocess.check_output([str(binary),'--known-keys',str(path),str(export)],text=True))
        keys=set();rows=0
        with export.open() as stream:
            for line in stream:
                degree,key=line.split('\t',1);keys.add((int(degree),tuple(tuple(p) for p in json.loads(key))));rows+=1
        if rows!=report['records'] or rows!=len(keys) or sha(path)!=digest:raise ValueError('native known-key mass/binding')
        cache_independent_known(path,digest,keys)
        return report


def audit_component_native(folder,binary,collect_records=False,deadline_epoch=None):
    summary=rows(Path(folder)/'SUMMARY.tsv')[0];n=int(summary['degree'])
    with tempfile.TemporaryDirectory(prefix='hurwitz-independent-keys-') as tmp:
        export=Path(tmp)/'keys.tsv'
        report=json.loads(subprocess.check_output([str(binary),'--component',str(folder),str(n),str(export)],text=True,timeout=remaining_timeout(deadline_epoch)))
        keys=[];records=[]
        with export.open() as stream:
            for line in stream:
                value,seed=line.rstrip('\n').split('\t',1)
                key=tuple(tuple(p) for p in json.loads(value));keys.append(key)
                if collect_records and report['closed']:records.append((key,seed))
        if len(keys)!=report['states'] or len(set(keys))!=len(keys):raise ValueError('native exported state mass')
        result=dict(summary=summary,independent_six_move_closure=report['closed'],
            state_keys=sorted(keys),table_sha256=sha(Path(folder)/'STATES.tsv'))
        if collect_records and report['closed']:result['_closed_records']=records
        return result


def audit_guided_native(folder, binary, collect_records=False,native_components=False,deadline_epoch=None):
    folder, binary = Path(folder), Path(binary)
    counts = rows(folder / 'COUNTS.tsv')
    if len(counts) != 1:
        raise ValueError('guided counts must have one row')
    n = int(counts[0]['degree'])
    result = subprocess.run([str(binary), '--guided', str(folder), str(n)],
                            check=True, capture_output=True, text=True,timeout=remaining_timeout(deadline_epoch))
    native = json.loads(result.stdout)
    if native['records'] != sum(native[k] for k in ['nonsmooth', 'closed', 'unresolved']):
        raise ValueError('native resolution mass')
    if (native['closed'], native['unresolved']) != (int(counts[0]['closed']), int(counts[0]['unresolved'])):
        raise ValueError('native status counts')
    components = [];identities=[];records=[]
    for row in native.pop('exceptions'):
        check = (audit_component_native(folder / row['folder'],binary,collect_records=collect_records,deadline_epoch=deadline_epoch)
                 if native_components else audit_component(folder / row['folder'],collect_records=collect_records))
        records.extend((n,key,seed) for key,seed in check.pop('_closed_records',[]))
        keys=check.pop('state_keys')
        if check['summary']['status'] != row['status']:
            raise ValueError('resolution/certificate mismatch')
        if parse_seed(check['summary']['seed'], n) != parse_seed(row['seed'], n):
            raise ValueError('resolution/component seed mismatch')
        components.append(check)
        if row['status']=='closed_smooth':
            identities.append(dict(degree=n,key=min(keys),summary=check['summary'],folder=row['folder']))
    banks={}
    if (folder/'INPUT_BANKS.tsv').exists():
        for line in (folder/'INPUT_BANKS.tsv').read_text().splitlines():
            flag,path=line.split('\t',1)
            if flag in banks or flag not in ['--known-bank','--negative-bank']:raise ValueError('input bank shape')
            banks[flag]=dict(path=path,sha256=sha(path))
            if flag=='--negative-bank':
                banks[flag]['audit']=json.loads(subprocess.check_output(
                    [str(binary),'--negative-bank',path,str(n)],text=True,timeout=remaining_timeout(deadline_epoch)))
    aliases=0;alias_ids=set()
    if (folder/'KNOWN.tsv').exists():
        previous=0
        for line in (folder/'KNOWN.tsv').read_text().splitlines():
            candidate,seed=line.split('\t',1);candidate=int(candidate)
            if '--known-bank' not in banks or not previous<candidate<=int(counts[0]['constructed']):raise ValueError('known alias input/order')
            if (n,state_key(parse_seed(seed,n))) not in known_keys(banks['--known-bank']['path']):raise ValueError('false known alias')
            previous=candidate;aliases+=1;alias_ids.add(candidate)
        if aliases>int(counts[0]['duplicate']):raise ValueError('alias/duplicate count')
        if alias_ids:
            import csv
            with (folder/'RESOLUTIONS.tsv').open() as stream:
                if any(int(row['candidate']) in alias_ids for row in csv.DictReader(stream,delimiter='\t')):raise ValueError('alias also resolved')
    exported=None
    if (folder/'NEGATIVE_CACHE.tsv').exists():
        exported=json.loads(subprocess.check_output(
            [str(binary),'--negative-bank',str(folder/'NEGATIVE_CACHE.tsv'),str(n)],text=True,timeout=remaining_timeout(deadline_epoch)))
        exported['sha256']=sha(folder/'NEGATIVE_CACHE.tsv')
    if (folder/'DESIGN.tsv').exists() and int(rows(folder/'DESIGN.tsv')[0]['known_hits'])!=aliases:raise ValueError('known alias count')
    source_names=['COUNTS.tsv', 'LEADS.tsv', 'RESOLUTIONS.tsv']
    triangle_filters=folder/'TRIANGLE_FILTERS.tsv'
    if triangle_filters.is_file():
        source_names.append('TRIANGLE_FILTERS.tsv')
        rows_t=rows(triangle_filters)
        if len(rows_t)!=1 or int(rows_t[0]['degree'])!=n:raise ValueError('triangle filter record degree/mass')
        if rows_t[0].get('third_method') not in ['random','joint','strong','mixed','centralizer','pairwise_low','middle_repair']:
            raise ValueError('unknown triangle third-entry method')
        exact_passport=rows_t[0].get('exact_passport','')
        if exact_passport:
            profile=tuple(map(int,exact_passport.split(',')))
            if len(profile)!=4 or tuple(sorted(profile))!=profile or profile[3]>=n:
                raise ValueError('invalid exact TRIANGLE passport binding')
        for key in ['profile_reject','equal_low_involution_reject','equal_low_two_special_reject']:
            if int(rows_t[0][key])<0:raise ValueError('negative exact-profile theorem counter')
        if int(rows_t[0]['passport_profiles'])<0:raise ValueError('negative triangle passport-profile count')
        if any(int(rows_t[0][key])<0 for key in ['budget_reject','shard_reject','genus_reject','high_fixed_point_reject']):
            raise ValueError('negative triangle filter count')
        short_fields=[f'short_{j}' for j in range(5)]
        if any(int(rows_t[0][key])<0 for key in short_fields):raise ValueError('negative short-screen count')
        short_total=sum(int(rows_t[0][key]) for key in short_fields)
        if short_total>int(counts[0]['short_reject']):raise ValueError('short-screen counts exceed total screen rejects')
        examples=folder/'TRIANGLE_SHORT_EXAMPLES.tsv'
        if not examples.is_file():raise ValueError('triangle short-witness examples missing')
        source_names.append('TRIANGLE_SHORT_EXAMPLES.tsv')
        expected_words=['1','4','0,1','0,4','3,1'];seen_factors=set()
        for example in rows(examples):
            j=int(example['factor']);word=example['word']
            if j not in range(5) or word!=expected_words[j] or j in seen_factors:
                raise ValueError('triangle short-witness label/order')
            t=parse_seed(example['seed'],n)
            if not values(t)['smooth']:raise ValueError('short witness seed was not initially smooth')
            moved=t
            for move in map(int,word.split(',')):moved=hurwitz(moved,move)
            if values(moved)['smooth']:raise ValueError('persisted short witness does not exhibit nonsmoothness')
            if int(rows_t[0][f'short_{j}'])==0:raise ValueError('short witness has zero recorded count')
            seen_factors.add(j)
        if seen_factors!={j for j in range(5) if int(rows_t[0][f'short_{j}'])>0}:
            raise ValueError('short-witness examples do not match the nonzero categories')
        deep_count=int(rows_t[0]['deep_reject'])
        if deep_count<0 or deep_count>int(counts[0]['short_reject']):raise ValueError('invalid deep-screen reject count')
        deep_examples=folder/'TRIANGLE_DEEP_EXAMPLES.tsv'
        if not deep_examples.is_file():raise ValueError('triangle deep-witness examples missing')
        source_names.append('TRIANGLE_DEEP_EXAMPLES.tsv')
        deep_rows=rows(deep_examples)
        if len(deep_rows)>min(32,deep_count):raise ValueError('triangle deep-witness example cap/mass')
        for example in deep_rows:
            if not example['failure'] or not example['word']:raise ValueError('empty deep failure/witness')
            t=parse_seed(example['seed'],n)
            if not values(t)['smooth']:raise ValueError('deep witness seed was not initially smooth')
            moved=t
            for move in map(int,example['word'].split(',')):moved=hurwitz(moved,move)
            if values(moved)['smooth']:raise ValueError('persisted deep witness does not exhibit nonsmoothness')
    result=dict(sampled=True, construction=counts[0], components=components,
                closed_identities=identities,known_aliases=aliases,input_banks=banks,exported_negative_bank=exported,
                independent_native=native, verifier_sha256=sha(binary),
                sources={name: sha(folder / name) for name in source_names})
    if triangle_filters.is_file():result['triangle_filters']=rows_t[0]
    if collect_records:
        result['closed_record_binding']=closed_record_digest(records)
        result['_closed_records']=records
    return result
