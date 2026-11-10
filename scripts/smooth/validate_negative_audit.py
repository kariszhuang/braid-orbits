#!/usr/bin/env python3
"""Full independent small audit, full degree10 frames and measured paired checks."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse, json, random, statistics, subprocess, time
from pathlib import Path
from run_negative_audit import Native, checked_frame, hashes, audit, save, sha
from run_smooth_degree9 import COLS, rows
from odd_order_spin import parse_seed

PROJECT=Path(__file__).resolve().parents[2]
OLD=Path('data/smooth_frontier_20261012_123427_295724')

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);args=p.parse_args();args.out.mkdir(parents=True,exist_ok=False)
    binary=PROJECT/'target/audit-fast/release/smooth_negative_audit';start=time.monotonic()
    subprocess.run(['cargo','test','--bin','smooth_negative_audit'],cwd=PROJECT,check=True,env=__import__('os').environ|{'CARGO_TARGET_DIR':str(PROJECT/'target/audit-fast'),'CARGO_BUILD_JOBS':'1'})
    subprocess.run(['python3','scripts/smooth/test_negative_audit.py'],cwd=PROJECT,check=True)
    root=args.out/'small';root.mkdir();small=audit(root,OLD/'validation/focused',binary,1);save(args.out/'SMALL_CHECKS.json',small)
    folder=OLD/'degree_10';frames=rows(folder/'FRAMES.tsv');counts=rows(folder/'FRAME_COUNTS.tsv');server=Native(binary,10);full=[]
    try:
        for fi in [0,1000,2000,4000]:
            row=frames[fi];a,b,_,_=parse_seed(f"[{row['a']}, {row['b']}, (), ()]",10);local=folder/'frames'/f'frame_{fi:05}'
            result=server.frame(local,a,b,True);checked_frame(result,local,fi,[int(counts[fi][c]) for c in COLS],True,10)
            result.update(frame=fi,sources=hashes(local,True));old=folder/'audit_frames'/f'frame_{fi:05}.json'
            if old.exists():
                check=json.loads(old.read_text());assert check['count']==result['count'] and check['sources']==result['sources']
            result['python_full_frame_compared']=old.exists();full.append(result)
    finally:server.close()
    save(args.out/'FULL_FRAMES.json',full)
    bench=Path('data/degree10_audit_benchmark_long_20261012.json');paired=json.loads(bench.read_text());rng=random.Random(101);intervals=[]
    for row in paired:
        ratios=[a/b for a,b in zip(row['timings']['python'],row['timings']['native'])]
        bootstrap=sorted(statistics.median(rng.choices(ratios,k=len(ratios))) for _ in range(10000))
        intervals.append(dict(frame=row['frame'],paired_median_speedup=statistics.median(ratios),descriptive95_interval=[bootstrap[249],bootstrap[9749]]))
    save(args.out/'CHECKS.json',dict(passed=True,binary_sha256=sha(binary),source_sha256=sha(PROJECT/'src/bin/smooth_negative_audit.rs'),
        regressions=dict(rust=2,python=3,literal_differential=16824),small_sha256=sha(args.out/'SMALL_CHECKS.json'),full_frames_sha256=sha(args.out/'FULL_FRAMES.json'),
        full_degree10_frames=4,full_degree10_rank_occurrences=sum(r['count'][0] for r in full),
        benchmark_sha256=sha(bench),benchmark_path=str(bench),paired=intervals,elapsed_seconds=time.monotonic()-start,
        limitation='Negative-only frames; every word replayed, no production kernels reused. Prefix timing excludes initialization, input hashing and complete population check; no discovery speed claim.'))
    print('CHECKS='+str(args.out/'CHECKS.json'))

if __name__=='__main__':main()
