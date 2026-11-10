"""Randomized check of Theorem G2 at single cusps: if a tuple has a block system
whose quotient has positive genus and a block-trivial entry, then at every cusp,
geometric smoothness (Lemma G1) implies combinatorial smoothness."""
import sys, random
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../common'))
from braid import compose, inverse, index, smooth, move
from geo_check import geometric_cusp, orbit_partition
random.seed(20261030)
def rperm(n): p = list(range(n)); random.shuffle(p); return tuple(p)
def wreath(q, k):
    # lift block permutation q (on m blocks of size k) with random fibre maps
    m = len(q); out = [0] * (m * k)
    for b in range(m):
        s = rperm(k)
        for i in range(k): out[b * k + i] = q[b] * k + s[i]
    return tuple(out)
tested = geo_only = both = 0
for trial in range(4000):
    m = random.choice([3, 4, 5]); k = random.choice([2, 3, 4])
    q1, q2 = rperm(m), rperm(m); q3 = inverse(compose(q1, q2))
    hq = 1 - m + (index(q1) + index(q2) + index(q3)) // 2
    if hq < 1: continue
    g1 = wreath(q1, k); g2 = wreath(q2, k); g4 = wreath(tuple(range(m)), k)
    g3 = compose(inverse(compose(g1, g2)), inverse(g4))
    t = (g1, g2, g3, g4)
    if orbit_partition(list(t), m * k)[1] != 1: continue
    for step in range(30):
        tested += 1
        g = geometric_cusp(t); c = smooth(t)
        if g and not c: geo_only += 1
        if g and c: both += 1
        t = move(t, random.randrange(6))
print('cusps tested', tested, 'geometric and combinatorial', both, 'geometric only (should be 0)', geo_only)
