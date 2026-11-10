#!/bin/bash
# Exhaustive fixed-quotient genus-3 campaign (detached). Stages run in order, each resumable.
# Stop: touch data/quotient_exhaust_20261030/STOP   (current job chunk finishes)
# Status: bash "$(dirname "$0")/status.sh"
ROOT=data/quotient_exhaust_20261030
BIN=$ROOT/bin/quotient_exhaust
THREADS=${THREADS:-10}
STAGES=(
  "h1k3a 1 3 3,4,5,6,7"
  "h0k2  0 2 3,4,5,6,7,8"
  "h0k3a 0 3 3,4,5"
  "h1k4a 1 4 3,4,5,6"
  "h1k5a 1 5 3,4"
  "h1k6a 1 6 3"
  "h1k7a 1 7 3"
  "h1k8a 1 8 3"
  "h1k6b 1 6 4"
  "h1k5b 1 5 5"
  "h1k7b 1 7 4"
  "h1k6c 1 6 5"
  "h1k5c 1 5 6"
  "h1k8b 1 8 4"
  "h1k7c 1 7 5"
  "h1k6d 1 6 6"
  "h1k3b 1 3 8"
  "h1k4b 1 4 7"
  "h0k3b 0 3 6"
  "h0k4a 0 4 3,4"
)
for st in "${STAGES[@]}"; do
  set -- $st; name=$1; g=$2; k=$3; ms=$4
  [ -e "$ROOT/STOP" ] && { echo "$(date) STOP seen before $name" >> $ROOT/campaign.log; exit 0; }
  if [ -e "$ROOT/$name/COMPLETE" ]; then continue; fi
  echo "$(date) start $name genus=$g k=$k m=$ms" >> $ROOT/campaign.log
  "$BIN" --out "$ROOT/$name" --genus $g --k $k --m $ms --threads $THREADS --walk 300 --cap 400000 --closure-seconds 120 2>> "$ROOT/$name.stderr"
  # complete only if every job is done and no STOP interrupted
  if [ ! -e "$ROOT/STOP" ] && grep -q '"stopped":false' "$ROOT/$name/RECEIPT.jsonl" 2>/dev/null; then
    touch "$ROOT/$name/COMPLETE"; echo "$(date) complete $name" >> $ROOT/campaign.log
  else echo "$(date) interrupted $name" >> $ROOT/campaign.log; exit 0; fi
done
echo "$(date) all stages complete" >> $ROOT/campaign.log
