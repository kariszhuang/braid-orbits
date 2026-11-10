#!/bin/bash
# exhaustive single-moving r=1 check, all class pairs with total index 4
K="$1"; CAP="${2:-300000}"
for p in "2 4" "4 2" "2 3,2" "3,2 2" "2 2,2,2" "2,2,2 2" "3 3" "3 2,2" "2,2 3" "2,2 2,2"; do
  set -- $p
  nice -n 15 ./target/release/elevel --exhaust-r1 --k $K --sig $1 --tau $2 --cap $CAP
done
