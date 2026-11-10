#!/usr/bin/env python3
"""Finalize an inactive genus-three controller without changing raw journals.

Audited completed counters remain separate from interrupted journal adapters.
Every durable unfinished seed is resolved before a successor is admissible.
This is recovery, not a population-completeness or discovery claim.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import csv
import fcntl
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess

from audit_guided_native import audit_guided_native
from run_genus_three import PROJECT, identity, inherited_known, next_word_bank, now, sha, write

HEADER = 'candidate\tstatus\tword\tfolder\tseed\n'
COUNTS = 'degree\tlane\tseed\tconstructed\tscope\ttransitive\tpassport\tshort_reject\tduplicate\tclosed\tunresolved\tseconds\n'


def inactive(root, pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        pass
    else:
        raise RuntimeError('controller PID still exists; no recovery or signal sent')
    for line in subprocess.check_output(['ps', '-axo', 'pid,command'], text=True).splitlines()[1:]:
        fields = line.strip().split(None, 1)
        if len(fields) != 2:
            continue
        try:
            argv = shlex.split(fields[1])
        except ValueError:
            continue
        if not argv:
            continue
        if argv[0] in [str(root/name) for name in ('smooth_genus_three', 'smooth_genus_three_compact', 'smooth_negative_audit')]:
            raise RuntimeError('bound producer/auditor still exists')
        if '--worker' in argv and str(root) in argv:
            raise RuntimeError('bound controller still exists')


def durable_journal(folder):
    """Accept complete durable lines only; malformed complete data fail closed."""
    leads = []
    if (folder/'LEADS.tsv').exists():
        for line in (folder/'LEADS.tsv').read_text().splitlines(keepends=True):
            if not line.endswith('\n'):
                break
            candidate, seed = line.rstrip('\n').split('\t', 1)
            candidate = int(candidate)
            if candidate <= (leads[-1][0] if leads else 0):
                raise ValueError('durable lead ordering')
            leads.append((candidate, seed))
    resolutions = []
    if (folder/'RESOLUTIONS.tsv').exists():
        text = (folder/'RESOLUTIONS.tsv').read_text()
        if text and not text.startswith(HEADER):
            if '\n' in text or not HEADER.startswith(text):
                raise ValueError('resolution header')
        for line in text.splitlines(keepends=True)[1:]:
            if not line.endswith('\n'):
                break
            fields = line.rstrip('\n').split('\t')
            at = len(resolutions)
            if len(fields) != 5 or at >= len(leads) or (int(fields[0]), fields[4]) != leads[at]:
                raise ValueError('resolution and durable lead differ')
            resolutions.append(fields)
    return leads, resolutions


def recover(root,finalized_attention=False):
    root = root.resolve()
    launch = json.loads((root/'LAUNCH.json').read_text())
    inactive(root, int(launch['pid']))
    if (root/'RECOVERY').exists():
        raise ValueError('already finalized or recovery exists; never overwrite a receipt')
    plan = json.loads((root/'PLAN.json').read_text())
    original = (root/'STATE.json').read_bytes()
    state = json.loads(original)
    original_final=(root/'FINISHED.json').read_bytes()if(root/'FINISHED.json').exists()else None
    if finalized_attention:
        final_state=json.loads(original_final)if original_final is not None else{}
        if not(state.get('status')==final_state.get('status')=='attention'and
            state.get('reason')==final_state.get('reason')=='audit_or_producer_failure'and
            state.get('epochs')==final_state.get('epochs')and not state.get('active_jobs')):
            raise ValueError('only finalized inactive audit/producer attention may use explicit recovery')
    elif original_final is not None:
        raise ValueError('already finalized; explicit attention recovery is required')
    for name, digest in json.loads((root/'BINDINGS.json').read_text()).items():
        if sha(root/name) != digest:
            raise ValueError('frozen binding changed: '+name)
    known = inherited_known(root)
    complete = [json.loads(line) for line in (root/'COMPLETED.jsonl').read_text().splitlines()]
    rows = {}
    for item in complete:
        folder = Path(item['folder'])
        if folder.parent != root/'jobs' or folder.name in rows or not (folder/'VERIFIED').exists():
            raise ValueError('completed job ownership or duplicate')
        check = json.loads((folder/'CHECKS.json').read_text())
        if check['verifier_sha256'] != sha(root/'smooth_negative_audit'):
            raise ValueError('completed verifier differs')
        for name, digest in check['sources'].items():
            if sha(folder/name) != digest:
                raise ValueError('completed journal changed')
        if check['construction'] != item['counts']:
            raise ValueError('completed counters differ')
        rows[folder.name] = item
    recovery = root/'RECOVERY'
    recovery.mkdir()
    (recovery/'ORIGINAL_STATE.json').write_bytes(original)
    if finalized_attention:
        (recovery/'ORIGINAL_FINISHED.json').write_bytes(original_final)
        shutil.copyfile(root/'CHECKPOINT.json',recovery/'ORIGINAL_CHECKPOINT.json')
    epochs = list(state['epochs'])
    if len(epochs) != plan.get('streams', plan['workers']):
        raise ValueError('stream shape differs')
    old = json.loads((root/'INHERITED_CHECKPOINT.json').read_text())
    latest = root/'banks'/f'words_{len(complete):08}.txt'
    words = [list(map(int, line.split(','))) for line in latest.read_text().splitlines() if line]
    recovered = []
    for folder in sorted((root/'jobs').iterdir()):
        if not folder.is_dir() or folder.name in rows:
            continue
        fields = folder.name.split('_')
        if len(fields) != 5 or fields[0] != 'job':
            raise ValueError('unexpected job folder')
        slot, epoch, degree = map(int, fields[1:4])
        if slot >= len(epochs) or epochs[slot] != epoch+1:
            raise ValueError('interrupted stream cursor differs')
        leads, resolutions = durable_journal(folder)
        adapter = recovery/folder.name
        adapter.mkdir()
        source_hashes = {p.name: sha(p) for p in folder.iterdir() if p.is_file()}
        for candidate, seed in leads[len(resolutions):]:
            replay = adapter/f'pending_{candidate}'
            binary = root/('smooth_genus_three_compact' if degree <= 25 else 'smooth_genus_three')
            subprocess.run([str(binary), '--d', str(degree), '--tuple', seed,
                '--out', str(replay), '--cap', str(plan['cap']),
                '--closure-seconds', str(plan['closure_seconds']), '--verifier', plan['verifier']],
                check=True, capture_output=True, text=True)
            with (replay/'SUMMARY.tsv').open() as stream:
                result = next(csv.DictReader(stream, delimiter='\t'))
            if result['status'] != 'nonsmooth':
                raise RuntimeError('pending seed remains exceptional; retained replay requires review: '+str(replay))
            resolutions.append([str(candidate), 'nonsmooth', result['word'], '', seed])
        if any(row[1] != 'nonsmooth' for row in resolutions):
            raise RuntimeError('interrupted exceptional journal requires component review')
        (adapter/'LEADS.tsv').write_text(''.join(f'{candidate}\t{seed}\n' for candidate, seed in leads))
        (adapter/'RESOLUTIONS.tsv').write_text(HEADER+''.join('\t'.join(row)+'\n' for row in resolutions))
        count = len(leads)
        (adapter/'COUNTS.tsv').write_text(COUNTS+f'{degree}\treplay\t0\t{leads[-1][0] if leads else 0}\t{count}\t{count}\t{count}\t0\t0\t0\t0\t0\n')
        check = audit_guided_native(adapter, root/'smooth_negative_audit', native_components=True)
        if check['independent_native']['nonsmooth'] != count:
            raise ValueError('recovery proof mass differs')
        write(adapter/'CHECKS.json', check)
        words.extend(list(map(int, row[2].split(','))) for row in resolutions)
        if any(sha(folder/name) != digest for name, digest in source_hashes.items()):
            raise ValueError('raw interrupted source changed during recovery')
        recovered.append(dict(folder=str(folder), raw_files=source_hashes,
            durable_leads=count, recovered_pending=count-len(durable_journal(folder)[1]),
            negative_audit=str(adapter/'CHECKS.json'), negative_audit_sha256=sha(adapter/'CHECKS.json')))
    bank = recovery/'words.txt'
    bank.write_text(''.join(','.join(map(str, word))+'\n' for word in next_word_bank(words, words)))
    totals = {target: sum(int(row['counts'][key]) for row in complete) for target, key in
        [('constructed','constructed'), ('scope','scope'), ('transitive','transitive'),
         ('closed_occurrences','closed'), ('unresolved','unresolved')]}
    components = [item for row in complete for item in row['components']]
    if any(item['id'] not in known for item in components):
        raise RuntimeError('unreconciled discovery requires review')
    destination=recovery/'ADMISSIBLE'if finalized_attention else root
    if finalized_attention:destination.mkdir()
    record = dict(passed=True, root=str(root), updated=now(), original_state_sha256=sha(recovery/'ORIGINAL_STATE.json'),
        completed_jobs=len(complete), interrupted=recovered, next_epochs=epochs,
        all_durable_leads_resolved=True, pending=0,
        finalized_attention=finalized_attention,admissible_predecessor=str(destination),
        qualification='completed counters only; interrupted adapters are not original proposal totals; no finite population coverage inferred')
    if finalized_attention:
        record['original_receipts']={name:sha(recovery/('ORIGINAL_'+name))for name in('STATE.json','FINISHED.json','CHECKPOINT.json')}
    write(recovery/'CHECKS.json', record)
    checkpoint = dict(next_epochs=epochs, scheduler_cursor=0,
        allocation_origin=old['allocation_origin'], allocation_weights=old['allocation_weights'],
        word_bank=str(bank), word_bank_sha256=sha(bank), counts=totals,
        completed_jobs=len(complete), new_components=[], boundary_cursor=state['boundary_cursor'],
        deadline_epoch=plan.get('deadline_epoch'),recovery_sha256=sha(recovery/'CHECKS.json'), pending='all durable interrupted starts independently negative; random partial proposals are not a census')
    write(destination/'CHECKPOINT.json', checkpoint)
    final = dict(status='stopped', reason='inactive_controller_recovered', updated=now(),
        completed_jobs=len(complete), counts=totals, epochs=epochs, new_closed_components=0,
        boundary_cursor=state['boundary_cursor'], boundary_total=504000,
        recovery_sha256=sha(recovery/'CHECKS.json'))
    write(destination/'FINISHED.json', final)
    write(destination/'STATE.json', dict(final, active_jobs=[], workers=0, pid=launch['pid']))
    if finalized_attention:
        adapter_plan=dict(plan,recovery_adapter=True,inheritance=dict(predecessor=str(root),kind='verified attention recovery adapter, not a scientific run'))
        write(destination/'PLAN.json',adapter_plan);shutil.copyfile(root/'KNOWN_IDS.json',destination/'KNOWN_IDS.json')
        write(destination/'ADAPTER_PROVENANCE.json',dict(original_root=str(root),recovery_receipt=str(recovery/'CHECKS.json'),
            recovery_receipt_sha256=sha(recovery/'CHECKS.json'),original_receipts=record['original_receipts']))
        write(destination/'BINDINGS.json',{name:sha(destination/name)for name in('PLAN.json','KNOWN_IDS.json','CHECKPOINT.json','FINISHED.json','STATE.json','ADAPTER_PROVENANCE.json')})
        for name,digest in record['original_receipts'].items():
            if sha(root/name)!=digest:raise ValueError('original attention receipt changed during recovery')
    return record


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('root', type=Path)
    p.add_argument('--finalize', action='store_true')
    p.add_argument('--finalize-attention',action='store_true',help='explicitly audit a finalized attention root into a separate stopped checkpoint adapter; original receipts stay unchanged')
    a = p.parse_args()
    if not(a.finalize or a.finalize_attention):
        inactive(a.root.resolve(), int(json.loads((a.root/'LAUNCH.json').read_text())['pid']))
        print('Inactive; use --finalize to audit durable journals and create a stopped checkpoint.')
        return
    with (PROJECT/'.run-control/.campaign.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        print(json.dumps(recover(a.root,finalized_attention=a.finalize_attention), indent=2))


if __name__ == '__main__':
    main()
