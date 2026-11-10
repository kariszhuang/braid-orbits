"""Immutable per-degree positive snapshots; publish only fully audited closures.

Thread-safe across all worker lanes. Every old input remains immutable while a
producer/auditor uses it. Hitting the RAM/row bound causes redundant work only;
no partial component becomes a positive or negative membership assertion.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import json
from pathlib import Path
import threading

from run_smooth_degree9 import sha
from audit_guided_native import cache_independent_known,known_keys,closed_record_digest


class VerifiedSmoothBank:
    def __init__(self, folder, initial, expected_sha, limit=250000):
        self.folder=Path(folder);self.folder.mkdir()
        if sha(initial)!=expected_sha:raise ValueError('initial verified bank changed')
        self.limit=limit;self.lock=threading.Lock();self.rows={};self.version={};self.snapshots={}
        self.initial_sha=expected_sha
        self.keys={}
        for n,key in known_keys(initial):self.keys.setdefault(n,set()).add(key)
        for line in Path(initial).read_text().splitlines():
            degree,seed=line.split('\t',1);n=int(degree)
            if not 2<=n<=63:raise ValueError('bank degree')
            self.rows.setdefault(n,set()).add(seed)
        if any(n<2 or n>63 for n in self.rows):raise ValueError('initial bank degree outside supported wide range')
        if any(len(v)>limit for v in self.rows.values()):raise ValueError('initial bank exceeds RAM bound')

    def _snapshot(self,n,proof):
        version=self.version.get(n,0);path=self.folder/f'degree_{n}_v{version:07}.tsv'
        path.write_text(''.join(f'{n}\t{s}\n' for s in sorted(self.rows.get(n,set()))))
        digest=sha(path)
        path.with_suffix('.json').write_text(json.dumps(dict(degree=n,version=version,
            rows=len(self.rows.get(n,set())),sha256=digest,initial_sha256=self.initial_sha,
            provenance=proof),indent=2)+'\n')
        self.snapshots[n]=(path,digest)
        cache_independent_known(path,digest,{(n,key) for key in self.keys.get(n,set())})
        return path,digest

    def snapshot(self,n):
        with self.lock:
            return self.snapshots.get(n) or self._snapshot(n,dict(kind='bound initial catalog'))

    def publish(self,job,checks,records):
        job=Path(job)
        # Require the just-written independent report and immutable verified
        # marker, not a producer's claims or an incomplete state table.
        if not (job/'VERIFIED').exists():raise ValueError('unverified positive bank publication')
        on_disk=json.loads((job/'CHECKS.json').read_text())
        if on_disk!=json.loads(json.dumps(checks)):raise ValueError('positive publication report differs')
        if checks.get('closed_record_binding')!=closed_record_digest(records):raise ValueError('verified records binding differs')
        closed=[c for c in checks['components'] if c['summary']['status']=='closed_smooth']
        expected=sum(int(c['summary']['states']) for c in closed)
        if len(records)!=expected:raise ValueError('only full closure records may enter bank')
        if any(not c['independent_six_move_closure'] for c in closed):raise ValueError('partial closure publication')
        identities=checks['closed_identities']
        if len(identities)!=len(closed):raise ValueError('closure identity mass')
        for c,item in zip(closed,identities):
            path=job/item['folder']/'STATES.tsv'
            if c['summary']!=item['summary'] or sha(path)!=c['table_sha256']:raise ValueError('published closure table changed/missing')
        n=int(checks['construction']['degree'])
        if any(degree!=n for degree,_,_ in records):raise ValueError('positive publication degree')
        with self.lock:
            old=self.rows.setdefault(n,set());add={seed for _,_,seed in records}-old
            if len(old)+len(add)>self.limit:
                return dict(degree=n,added=0,at_bound=True,rows=len(old))
            if not add:return dict(degree=n,added=0,rows=len(old))
            old.update(add);self.keys.setdefault(n,set()).update(key for _,key,_ in records)
            self.version[n]=self.version.get(n,0)+1
            path,digest=self._snapshot(n,dict(kind='independently verified six-move closures',
                job=str(job),checks_sha256=sha(job/'CHECKS.json'),added=len(add),
                previous_version=self.version[n]-1,table_sha256=[c['table_sha256'] for c in closed]))
            return dict(degree=n,added=len(add),rows=len(old),path=str(path),sha256=digest)
