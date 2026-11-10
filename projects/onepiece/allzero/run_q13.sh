#!/bin/zsh
cd "$(dirname "$0")"
python3 make_collineations.py 13 > q13_gen.log 2>&1
../target/release/collin_big --coll coll_q13.txt --orders 2,2,2,0 --threads 4 --out big_q13_2_2_2_any.jsonl 2> big_q13.log
echo done >> big_q13.log
