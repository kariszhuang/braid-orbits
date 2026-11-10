#!/bin/bash
# one-core queue: wait for r1 k10, then exhaustive r2 for k=7,8,9 (low priority)
cd "$(dirname "$0")"
while pgrep -f "elevel --exhaust-r1" >/dev/null; do sleep 15; done
for k in 7 8 9; do bash r2_sweep.sh $k 300000 > elevel_out/r2_k$k.txt 2>&1; done
echo done > elevel_out/QUEUE_DONE
