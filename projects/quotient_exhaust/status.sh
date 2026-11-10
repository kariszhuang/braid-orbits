#!/bin/bash
ROOT=data/quotient_exhaust_20261030
pgrep -f "$ROOT/bin/quotient_exhaust" >/dev/null && echo RUNNING || echo NOT-RUNNING
tail -4 $ROOT/campaign.log
python3 - "$ROOT" <<'PY'
import json,sys,os,glob,collections
root=sys.argv[1]
for d in sorted(glob.glob(root+'/h*')):
    if not os.path.isdir(d): continue
    done=0; agg=collections.defaultdict(lambda:[0]*7)
    seen=set()
    for l in open(d+'/done.jsonl') if os.path.exists(d+'/done.jsonl') else []:
        j=json.loads(l)
        if j['job'] in seen: continue
        seen.add(j['job']); done+=1
        for c,v in j.get('counts_built_walkpass_dup_nonsmooth_capped_closed_closedpos',{}).items():
            for i in range(7): agg[c][i]+=v[i]
    pend=[json.loads(l).get('jobs_pending') for l in open(d+'/CONFIG.jsonl')][0] if os.path.exists(d+'/CONFIG.jsonl') else '?'
    tot=[sum(v[i] for v in agg.values()) for i in range(7)]
    pos=len(set(l for l in open(d+'/POSITIVE_CANDIDATES.jsonl'))) if os.path.exists(d+'/POSITIVE_CANDIDATES.jsonl') else 0
    print(f"{os.path.basename(d):7s} {'COMPLETE' if os.path.exists(d+'/COMPLETE') else 'partial ':8s} jobs {done}/{pend} built {tot[0]} walkpass {tot[1]} nonsmooth {tot[3]} capped {tot[4]} closed {tot[5]} closed_lambda>0 {tot[6]}  POSITIVES {pos}")
PY
