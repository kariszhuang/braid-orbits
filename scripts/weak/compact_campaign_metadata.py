"""Preserve then pack AppleDouble metadata of finished historical weak jobs.

Scientific outputs remain at their original paths. Real resource forks/extended
attributes are not discarded: every sidecar is archived and individually hashed
before removal, using the existing non-overwriting restoration implementation.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import fcntl
import json
import os
from pathlib import Path
import shutil
import struct

from compact_smooth_outputs import compact,restore,sha

PROJECT=Path(__file__).resolve().parents[2]


def metadata_plan(root,protected=(),progress=None):
    root=root.resolve();protected={Path(p).resolve() for p in protected};files=[];skipped=0;attempts=0;allocated=logical=0
    jobs=root/'jobs'
    if jobs.is_symlink() or not jobs.is_dir():raise ValueError('real campaign jobs directory required')
    def add(path):
        nonlocal skipped,allocated,logical
        sibling=path.with_name(path.name[2:])
        # root/jobs/job/attempt and both siblings have been checked against
        # symlinks. Re-resolving every ancestor per sidecar is unnecessary.
        if path.is_symlink() or not path.is_file() or not sibling.exists() or sibling.is_symlink() or sibling in protected or path in protected:
            skipped+=1;return
        with path.open('rb') as stream:header=stream.read(26)
        if len(header)!=26 or struct.unpack('>II',header[:8])!=(0x00051607,0x00020000):
            skipped+=1;return
        stat=path.stat();allocated+=stat.st_blocks*512;logical+=stat.st_size;files.append(path)
    for number,job in enumerate(jobs.iterdir(),1):
        if progress and number%1000==0:progress(number,len(files))
        if job.is_symlink() or not job.is_dir() or not (job/'DONE').is_file():continue
        for attempt in job.iterdir():
            if not attempt.name.startswith('attempt_') or attempt.is_symlink() or not attempt.is_dir() or not (attempt/'VERIFIED').is_file():continue
            if any((attempt/name).exists() and (attempt/name).stat().st_size for name in ['hard_candidates.txt','survivors.txt']):
                skipped+=1;continue
            attempts+=1
            for path in attempt.iterdir():
                if path.name.startswith('._'):add(path)
            sidecar=attempt.with_name('._'+attempt.name)
            if sidecar.is_file():add(sidecar)
    if len(files)!=len(set(files)) or any(not f.is_relative_to(root/'jobs') for f in files):raise ValueError('metadata target escape/duplicate')
    return sorted(files),dict(verified_finished_attempts=attempts,metadata_files=len(files),skipped=skipped,
        allocated_bytes=allocated,logical_bytes=logical)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path)
    p.add_argument('--execute',action='store_true');p.add_argument('--restore',type=Path)
    args=p.parse_args();root=args.root.resolve()
    if root.parent!=Path('data') or not root.name.startswith('campaign_'):raise ValueError('exact historical campaign only')
    lock=(PROJECT/'.run-control/.campaign.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    state=json.loads((root/'STATE.json').read_text())
    if state['status'] not in ['stopped','complete','completed']:raise ValueError('campaign must be inactive')
    try:os.kill(state['pid'],0)
    except ProcessLookupError:pass
    else:raise ValueError('campaign PID alive/reused')
    if args.restore:restore(args.restore.resolve(),root);print('RESTORED',flush=True);return
    manifest=json.loads((PROJECT/'evidence/manifest.json').read_text())
    files,summary=metadata_plan(root,[(PROJECT/e['path']).resolve() for e in manifest['sources']],
        progress=lambda jobs,files:print('SCANNED_JOBS',jobs,'METADATA',files,flush=True))
    print(json.dumps(summary,indent=2),flush=True)
    if not args.execute:return
    if args.out is None or args.out.exists() or args.out.resolve().parent!=root.parent:raise ValueError('new explicit sibling archive folder required')
    args.out.mkdir();free_before=shutil.disk_usage(root).free
    (args.out/'PLAN.json').write_text(json.dumps(dict(root=str(root),selection=summary,
        ordinary_scientific_files_untouched=True,mode='recoverable AppleDouble metadata compaction'),indent=2)+'\n')
    batches=[]
    for offset in range(0,len(files),25000):
        result=compact(root,args.out/f'batch_{offset//25000:04}',files[offset:offset+25000],[])
        batches.append(dict(folder=str(args.out/f'batch_{offset//25000:04}'),archive_sha256=result['archive_sha256'],
            removed_files=result['removed_files'],original_allocated_bytes=result['original_allocated_bytes'],archive_bytes=result['archive_bytes']))
        report=dict(passed=True,root=str(root),selection=summary,batches=batches,
            removed_metadata_files=sum(b['removed_files'] for b in batches),ordinary_scientific_files_untouched=True,
            free_bytes_before=free_before,free_bytes_after=shutil.disk_usage(root).free,
            recovery='compact_campaign_metadata.py --root EXACT_ROOT --restore BATCH/REPEATED_OUTPUTS.tar.gz; refuses conflicting existing files')
        temp=args.out/'RESULTS.tmp';temp.write_text(json.dumps(report,indent=2)+'\n');temp.replace(args.out/'RESULTS.json')
        print('BATCH_VERIFIED',offset//25000,'removed',result['removed_files'],'free_GiB',round(report['free_bytes_after']/2**30,2),flush=True)


if __name__=='__main__':main()
