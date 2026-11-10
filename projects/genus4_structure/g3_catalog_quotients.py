"""For every catalogued source-genus-3 all-smooth component: genus-1 block quotients (with a block-trivial
entry) and the minimal genus-1 block size; components with no genus-1 quotient are listed separately."""
import json,sys,collections
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../common'))
from braid import parse
from blocks import analyse
d=json.load(open('../../exports/latest-discoveries-20261023-v12/compact-v12/all-smooth.json'))['components']
g3=[c for c in d if c['genus']==3]
stat=collections.Counter(); none=[]
for c in g3:
    t=parse(c['tuple'],c['degree'])
    an=analyse(t)
    g1=[a for a in an if a['qgenus']==1]
    kmin=min([a['k'] for a in g1],default=None)
    triv=[a for a in g1 if a['trivial']]
    stat[(c['degree'],kmin,bool(triv),c['isotriviality']['status'])]+=1
    if kmin is None: none.append((c['degree'],c['states'],c['quotient_genera']))
print('genus-3 components:',len(g3))
for k,v in sorted(stat.items(),key=lambda x:(x[0][0],x[0][1] or 0)): print(' degree %s  min genus-1 block size %s  block-trivial entry %s  iso %s : %d'%(k[0],k[1],k[2],k[3],v))
print('components with NO genus-1 quotient:',len(none),none[:10])
