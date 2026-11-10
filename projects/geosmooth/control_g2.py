"""Control for test_g2.py: random transitive positive-genus tuples WITHOUT the block
hypothesis do have single cusps that are geometrically but not combinatorially smooth."""
import sys, random
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../common'))
from braid import compose, inverse, index, smooth, move
from geo_check import geometric_cusp, orbit_partition
random.seed(7)
def rperm(n): p = list(range(n)); random.shuffle(p); return tuple(p)
tested = geo_only = 0
for trial in range(6000):
    n = random.choice([6, 8, 9, 10, 12])
    g1, g2, g4 = rperm(n), rperm(n), rperm(n)
    g3 = compose(inverse(compose(g1, g2)), inverse(g4))
    t = (g1, g2, g3, g4)
    if orbit_partition(list(t), n)[1] != 1: continue
    if sum(index(p) for p in t) < 2 * n: continue
    for step in range(5):
        tested += 1
        if geometric_cusp(t) and not smooth(t): geo_only += 1
        t = move(t, random.randrange(6))
print('cusps tested', tested, 'geometric-only cusps', geo_only)
