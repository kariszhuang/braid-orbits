#!/bin/bash
# after the r2 queue: exhaustive r=3 (four transpositions) for k=4..9, low priority, one core
cd "$(dirname "$0")"
while [ ! -e elevel_out/QUEUE_DONE ]; do sleep 20; done
for k in 4 5 6 7 8 9; do nice -n 15 ./target2/release/elevel --exhaust-gen --k $k --types "2;2;2;2" --cap 300000 > elevel_out/r3_k$k.txt 2>&1; done
echo done > elevel_out/QUEUE_R3_DONE
