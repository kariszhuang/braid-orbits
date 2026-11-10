#!/bin/bash
# after r3 queue: exhaustive r=1 at k=11 (one core, nice; K=16 binary to save memory)
cd "$(dirname "$0")"
while [ ! -e elevel_out/QUEUE_R3_DONE ]; do sleep 30; done
for p in "2 4" "4 2" "2 3,2" "3,2 2" "2 2,2,2" "2,2,2 2" "3 3" "3 2,2" "2,2 3" "2,2 2,2"; do
  set -- $p
  nice -n 15 ./target2/release/elevel --exhaust-r1 --k 11 --sig $1 --tau $2 --cap 300000 >> elevel_out/r1_k11.txt 2>&1
done
echo done > elevel_out/QUEUE_K11_DONE
