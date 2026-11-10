#!/usr/bin/env python3
"""Read-only genus-three dashboard. Refresh every five seconds; Ctrl-C exits.

No search imports, network access, locks, signals to workers or orbit scans.
Counters cover completed audited jobs, so advance at job boundaries (~45s).
"""
import argparse
from collections import Counter, deque
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys
import time

PROJECT = Path(__file__).resolve().parents[2]


def stamp(text):
    value = datetime.fromisoformat(text)
    if value.tzinfo is None:
        raise ValueError('status timestamp needs a timezone')
    return value.timestamp()


class JsonCache:
    def __init__(self):
        self.values = {}
        self.reads = 0

    def read(self, path, optional=False):
        path = Path(path)
        try:
            info = path.stat()
        except FileNotFoundError:
            self.values.pop(path, None)
            if optional:
                return None
            raise
        key = (info.st_ino, info.st_mtime_ns, info.st_size)
        old = self.values.get(path)
        if old is not None and old[0] == key:
            return old[1]
        # Only small control/index files enter this cache, never state tables.
        value = json.loads(path.read_text())
        self.reads += 1
        self.values[path] = (key, value)
        return value


def process_snapshot(root, controller):
    """Inspect only the controller and its children, not every command line."""
    children = subprocess.run(['pgrep','-P',str(controller)],
        capture_output=True,text=True,timeout=3)
    if children.returncode not in [0,1]:
        raise OSError('cannot enumerate controller children')
    pids = [str(controller)] + [p for p in children.stdout.split() if p.isdigit()]
    result = subprocess.run(['ps','-p',','.join(pids),'-o','pid=,ppid=,%cpu=,command='],
        capture_output=True,text=True,timeout=3)
    if result.returncode not in [0,1]:
        raise OSError('cannot inspect selected processes')
    alive = False; search = 0; audits = 0; cpu = 0.
    for line in result.stdout.splitlines():
        fields = line.strip().split(None, 3)
        if len(fields) != 4:
            continue
        pid, ppid = int(fields[0]), int(fields[1]); command = fields[3]
        if pid == controller and f'--worker {root}' in command and 'run_genus_three.py' in command:
            alive = True
        if ppid != controller:
            continue
        if any(command.startswith(str(root/name)+' ')for name in ['smooth_genus_three','smooth_genus_three_compact']):
            search += 1; cpu += float(fields[2])
        elif command.startswith(str(root/'smooth_negative_audit')+' '):
            audits += 1; cpu += float(fields[2])
    return dict(controller_alive=alive, search_workers=search, audit_workers=audits,
                worker_cpu_percent=cpu)


class Monitor:
    """One sample interface hides caching, ancestry and reporting conventions."""
    def __init__(self, project=PROJECT, processes=process_snapshot, clock=time.time):
        self.project = Path(project)
        self.processes = processes; self.clock = clock
        self.files = JsonCache(); self.catalog = None; self.history = deque()
        self.admissions = {}

    def baseline(self, path):
        path = Path(path); info = path.stat()
        key = (str(path), info.st_ino, info.st_mtime_ns, info.st_size)
        if self.catalog is None or self.catalog[0] != key:
            # Read the 45MB publication once on startup/change, retain only IDs.
            # Do not retain/reparse the entire scientific catalogue each tick.
            with path.open() as stream:
                data = json.load(stream)
            self.catalog = (key, dict(version=data['dataset_version'],
                ids={e['id'] for e in data['components']},
                genus3=sum(e['genus'] == 3 for e in data['components'])))
        return self.catalog[1]

    def discoveries(self, root, published):
        found = {}; roots = set(); warnings = []; pending = 0
        while root not in roots:
            roots.add(root)
            raw = self.files.read(root/'NEW_COMPONENTS.json', optional=True)
            cached = self.admissions.get(root)
            if cached is not None and cached[0] is raw:
                entries = cached[1]
            else:
                entries = []; markers = {}; missing = 0
                for item in raw or []:
                    parent = Path(item['raw_folder']).parent
                    if parent not in markers:
                        markers[parent] = (parent/'VERIFIED').is_file()
                    if markers[parent]:
                        entries.append(item)
                    else:
                        missing += 1
                # Completed receipts are immutable. Retry missing markers, but
                # otherwise never restat hundreds of old certificates per tick.
                if not missing:
                    self.admissions[root] = (raw, entries)
                pending += missing
            for item in entries:
                if item.get('genus') != 3 or item['id'] in published:
                    continue
                # NEW_COMPONENTS is controller-admitted metadata, not a raw hit.
                # Do not perform heavyweight scientific re-verification here.
                cid = item['id']; status = item.get('original_isotrivial')
                label = 'isotrivial' if status is True else 'non_isotrivial' if status is False else 'unknown'
                if cid in found and found[cid]['classification'] != label:
                    warnings.append('Conflicting classification: '+cid)
                    label = 'unknown'
                found[cid] = dict(classification=label, degree=item['degree'])
            inherited = self.files.read(root/'INHERITED_CHECKPOINT.json', optional=True)
            if inherited is None:
                # A recovered attention root is reached through a separate,
                # stopped adapter. Its small proof-bound provenance is an
                # ancestry link, not a positive record or new classification.
                plan=self.files.read(root/'PLAN.json',optional=True)
                if not(plan or {}).get('recovery_adapter'):
                    break
                provenance=self.files.read(root/'ADAPTER_PROVENANCE.json',optional=True)
                record=(self.files.read(Path(provenance['recovery_receipt']),optional=True)if provenance else None)
                if not(record and record.get('passed')is True and record.get('finalized_attention')is True and
                    record.get('all_durable_leads_resolved')is True and record.get('pending')==0 and
                    record.get('admissible_predecessor')==str(root)and record.get('root')==provenance.get('original_root')):
                    warnings.append('Incomplete recovery ancestry: '+str(root));break
                inherited=dict(predecessor=provenance['original_root'])
            root = Path(inherited['predecessor'])
            if not root.is_dir():
                warnings.append('Missing predecessor: '+str(root)); break
        else:
            warnings.append('Predecessor cycle; totals may be incomplete')
        if pending:
            warnings.append(f'{pending} index records lack VERIFIED markers and are not counted')
        counts = Counter(e['classification'] for e in found.values())
        return dict(total=len(found), isotrivial=counts['isotrivial'],
            non_isotrivial=counts['non_isotrivial'], unknown=counts['unknown'],
            by_degree=dict(sorted(Counter(e['degree'] for e in found.values()).items())),
            warnings=warnings)

    def sample(self):
        now = self.clock()
        pointer = self.files.read(self.project/'.run-control/genus-three-current.json')
        root = Path(pointer['root'])
        plan = self.files.read(root/'PLAN.json'); state = self.files.read(root/'STATE.json')
        launch = self.files.read(root/'LAUNCH.json')
        baseline = self.baseline(plan['catalog']['path'])
        new = self.discoveries(root, baseline['ids'])
        warnings = list(new.pop('warnings'))
        try:
            processes = self.processes(root, int(state['pid']))
        except (OSError, subprocess.SubprocessError, ValueError) as exc:
            processes = dict(controller_alive=None, search_workers=None, audit_workers=None,
                             worker_cpu_percent=None)
            warnings.append('Process inspection unavailable: '+str(exc))
        recorded = state['status']
        if recorded in ['running','queued','stopping'] and processes['controller_alive'] is False:
            status = 'NOT RUNNING (saved status is stale)'
        elif recorded in ['running','queued','stopping'] and processes['controller_alive'] is None:
            status = recorded.upper()+' (process unconfirmed)'
        else:
            status = recorded.upper()
        counts = state.get('counts', {})
        completed_draws = int(counts.get('constructed', 0))
        elapsed = max(0., now-stamp(launch['updated']))
        active = recorded in ['running','queued','stopping'] and processes['controller_alive'] is True
        if not active:
            final = self.files.read(root/'FINISHED.json',optional=True)
            elapsed = (float(final['seconds']) if final is not None and 'seconds' in final
                       else max(0.,stamp(state['updated'])-stamp(launch['updated'])))
        # Mutable counters advance only after auditing. Do not claim per-tuple
        # instantaneous telemetry from those stepwise values.
        avg = completed_draws/elapsed if elapsed > 0 else None
        if self.history and self.history[-1][0] != root:
            self.history.clear()
        self.history.append((root, now, completed_draws))
        while len(self.history)>2 and self.history[1][1] <= now-120:
            self.history.popleft()
        _, began, first = self.history[0]
        recent = ((completed_draws-first)/(now-began)
                  if active and now-began >= 60 and completed_draws >= first else None)
        return dict(status=status, saved_status=recorded, root=str(root),
            controller_pid=state['pid'], processes=processes,
            saved_counters_age_seconds=max(0.,now-stamp(state['updated'])),
            saved_counters_updated=state['updated'], completed_jobs=state.get('completed_jobs',0),
            proposal_draws=completed_draws, average_draws_per_second=avg,
            recent_draws_per_second=recent,
            completed_qualified=int(counts.get('scope',0)),
            boundary=(dict(next=int(state['boundary_cursor']),total=int(state['boundary_total']))
                      if 'boundary_cursor'in state and 'boundary_total'in state else None),
            unresolved=int(counts.get('unresolved',0)),
            published_baseline_version=baseline['version'],
            published_genus3=baseline['genus3'], new_genus3=new,
            known_genus3_total=baseline['genus3']+new['total'], warnings=warnings)


def rate(value):
    return '—' if value is None else f'{value/1_000_000:,.2f}M draws/s'


def render(value, interval, color=False):
    def ink(text, code):
        return f'\033[{code}m{text}\033[0m' if color else str(text)
    p = value['processes']; new = value['new_genus3']
    workers = (f"{p['search_workers']} workers · {p['audit_workers']} auditing"
               if p['search_workers'] is not None else 'unconfirmed')
    cpu = f" · CPU {p['worker_cpu_percent']:.0f}%" if p['worker_cpu_percent'] is not None else ''
    running = value['status'] == 'RUNNING'
    average = rate(value['average_draws_per_second'])
    mode = 'avg' if running else 'historical avg'
    updated = datetime.fromisoformat(value['saved_counters_updated']).astimezone().strftime('%H:%M:%S')
    age = value['saved_counters_age_seconds']
    lines = [ink('GENUS 3 · ALL-SMOOTH', '1'),
        ink(value['status'], '1;32' if running else '1;31')+f'  {workers}{cpu}',
        '',
        f"Speed              {ink(average, '1;36')} {mode}",
        f"Recent             {ink(rate(value['recent_draws_per_second']), '36')}",
        '',
        f"New (vs v{value['published_baseline_version']})       {ink(format(new['total'], ','), '1;36')}",
        f"Non-isotrivial     {ink(format(new['non_isotrivial'], ','), '1;32' if new['non_isotrivial'] else '1;33')}",
        f"Isotrivial         {format(new['isotrivial'], ',')}",
        f"Unknown            {ink(format(new['unknown'], ','), '1;33' if new['unknown'] else '0')}",
        f"Total genus 3      {format(value['known_genus3_total'], ',')}",
        '',
        f"Jobs {value['completed_jobs']:,} · Draws {value['proposal_draws']:,} · Caps "+ink(value['unresolved'], '1;33' if value['unresolved'] else '0'),
        ink(f"Updated {updated} · Age {age:.0f}s · Refresh {interval:g}s", '1;33' if age>30 else '0')]
    if value.get('boundary'):
        b=value['boundary']
        lines.insert(-1,f"Degree16 matching  {ink(format(b['next'],','), '1;36')} / {b['total']:,} constructed")
    lines.extend(ink('WARNING: '+warning, '1;31') for warning in value['warnings'])
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--once',action='store_true',help='print one snapshot and exit')
    parser.add_argument('--json',action='store_true',help='one machine-readable snapshot')
    parser.add_argument('--interval',type=float,default=5,help='refresh seconds, minimum 2 (default 5)')
    args = parser.parse_args()
    if args.interval < 2 or not args.interval < float('inf'):
        parser.error('interval must be finite and at least two seconds')
    try:
        os.nice(10)  # Only this viewer; never change worker scheduling.
    except OSError:
        pass
    monitor = Monitor(); once = args.once or args.json or not sys.stdout.isatty()
    colored = sys.stdout.isatty() and 'NO_COLOR' not in os.environ and os.environ.get('TERM') != 'dumb'
    try:
        while True:
            try:
                value = monitor.sample()
                content = json.dumps(value,indent=2) if args.json else render(value,args.interval,color=colored)
            except (OSError, ValueError, KeyError, TypeError) as exc:
                content = json.dumps(dict(status='UNAVAILABLE',error=str(exc))) if args.json else 'Status unavailable: '+str(exc)
                if once:
                    print(content); return 1
            if not once:
                print('\033[H\033[2J',end='')
            print(content,flush=True)
            if once:
                return 0
            time.sleep(args.interval)
    except KeyboardInterrupt:
        print(); return 0


if __name__ == '__main__':
    sys.exit(main())
