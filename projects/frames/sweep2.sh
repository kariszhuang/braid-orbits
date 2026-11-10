#!/bin/bash
# frames sweep for larger k (one core). Waits for the k=5,6 sweep to finish first.
cd "$(dirname "$0")"
while [ "$(wc -l < out/k6_r3.txt 2>/dev/null || echo 0)" -lt 25 ]; do sleep 20; done
for k in 7 8 9; do for r in 1 2; do
  ./target/release/frames --k $k --trials 10 --r $r --seed $((k*100+r)) --cap 3000000 > out/k${k}_r${r}.txt 2>&1
done; done
echo finished > out/SWEEP2_DONE
