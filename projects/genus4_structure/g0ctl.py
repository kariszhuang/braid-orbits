import json,sys,collections
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../common'))
from braid import parse
from blocks import analyse
d=json.load(open('../../exports/latest-discoveries-20261023-v12/compact-v12/all-smooth.json'))['components']
cnt=collections.Counter(); ex={}
for c in d:
    if c['genus'] not in (2,3) or c['degree']>24: continue
    t=parse(c['tuple'],c['degree'])
    for a in analyse(t):
        if a['trivial']:
            key=(c['genus'],a['qgenus'],a['k'],a['m'])
            cnt[key]+=1; ex.setdefault(key,(c['degree'],c['states'],c['tuple'],a['rel']))
for k,v in sorted(cnt.items()): print(v,k, ex[k][0], ex[k][1])
