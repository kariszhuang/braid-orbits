#!/bin/bash
# genus-0 quotient lane, single moving point, random samples (one core, nice)
cd "$(dirname "$0")"
for kr in "6 4" "6 5" "7 4" "7 5" "8 4" "8 5" "9 3" "9 4" "10 3" "10 4" "12 3" "12 4"; do
  set -- $kr
  nice -n 15 ./target3/release/elevel --p1 --k $1 --r $2 --trials 200 --seed $(( $1 * 10 + $2 )) --cap 300000
done
