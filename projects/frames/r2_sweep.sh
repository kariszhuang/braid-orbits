#!/bin/bash
# exhaustive single-moving r=2 check (sigma moving), all class triples with total index 4
K="$1"; CAP="${2:-300000}"
for p in "3 2 2" "2,2 2 2" "2 3 2" "2 2,2 2" "2 2 3" "2 2 2,2"; do
  set -- $p
  nice -n 15 ./target/release/elevel --exhaust-r2 --k $K --sig $1 --tau1 $2 --tau2 $3 --cap $CAP
done
