try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
from itertools import combinations, permutations
import json
from math import lcm
from pathlib import Path
import random
import re
import unittest

import genus_three_map_core as maps
import genus_three_double_marks as double


def perm(n,edges):
    a=list(range(n))
    for x,y in edges:a[x],a[y]=y,x
    return tuple(a)


class DoubleTests(unittest.TestCase):
    def test_entire_forest_can_share_a_colour_and_matches_all_S7(self):
        a=perm(7,((0,1),(2,3)));shared=perm(7,((1,2),))
        self.assertTrue(double.forest(((0,1),(2,3)),maps.labels(shared)))
        self.assertEqual(maps.index(maps.compose(a,shared)),3)
        for v in permutations(range(7)):
            literal=maps.index(maps.compose(a,v))==2+maps.index(v)
            self.assertEqual(double.forest(((0,1),(2,3)),maps.labels(v)),literal)

    def test_one_face_exact_noncrossing_population_up_to_actual_rotations(self):
        for n in range(4,9):
            p=tuple(range(1,n))+(0,);expected=set()
            for points in combinations(range(n),4):
                x,y,z,w=points
                for edges in(((x,y),(z,w)),((x,z),(y,w)),((x,w),(y,z))):
                    a=perm(n,edges);c=maps.compose(a,p)
                    if maps.index(a)+maps.index(c)!=maps.index(p):continue
                    key=min(tuple(sorted(tuple(sorted(((u+t)%n,(v+t)%n)))for u,v in edges))for t in range(n))
                    expected.add(key)
            actual={tuple(sorted(tuple(sorted(edge))for edge in edges))for edges in double._one_face(tuple(range(n)))}
            self.assertEqual(actual,expected)

    def test_two_faces_never_assume_independent_rotations(self):
        for m,n in((2,3),(3,4),(4,4),(4,6),(5,5)):
            expected=set()
            for a in combinations(range(m),2):
                for b in combinations(range(n),2):
                    key=min((tuple(sorted(((a[0]+t)%m,(a[1]+t)%m))),
                             tuple(sorted(((b[0]+t)%n,(b[1]+t)%n))))for t in range(lcm(m,n)))
                    expected.add(key)
            actual=set(double._two_faces(tuple(range(m)),tuple(range(n))))
            self.assertEqual(actual,expected)
        # In equal-size faces, relative phase is genuinely retained.
        self.assertGreater(len(list(double._two_faces(tuple(range(4)),tuple(range(4))))),3)

    def test_phase_palettes_equal_literal_synchronous_powers(self):
        rng=random.Random(202610250249)
        for sizes in((7,),(3,4),(4,4),(2,6)):
            n=sum(sizes);p=list(range(n));at=0
            for size in sizes:
                for x in range(size):p[at+x]=at+(x+1)%size
                at+=size
            edges=list(double.chord_pairs(p));period=lcm(*sizes)
            for _ in range(30):
                b=list(range(n));e=list(range(n));rng.shuffle(b);rng.shuffle(e)
                for pair in edges:
                    literal=True;moving=pair
                    for _ in range(period):
                        literal&=double.forest(moving,maps.labels(b))and double.forest(moving,maps.labels(e))
                        moving=tuple((p[x],p[y])for x,y in moving)
                    self.assertEqual(double.first_twists(pair,p,b,e),literal)

    def test_middle_exact_whole_forest_matches_literal_square_moves(self):
        rng=random.Random(202610250249);checked=0
        for n in(14,19,23):
            for _ in range(12):
                b,e=maps.grow(n,rng);p=maps.inverse(maps.compose(b,e))
                pairs=list(double.chord_pairs(p));rng.shuffle(pairs)
                for edges in pairs[:8]:
                    a=perm(n,edges);c=maps.compose(a,p);t=(a,c,b,e);colour=maps.labels(c)
                    positive=double.forest(tuple((e[x],e[y])for x,y in edges),colour)
                    ai,ei=maps.inverse(a),maps.inverse(e)
                    negative=double.forest(tuple((ai[ei[x]],ai[ei[y]])for x,y in edges),colour)
                    self.assertEqual(positive,maps.smooth(maps.braid(maps.braid(t,2),2)))
                    self.assertEqual(negative,maps.smooth(maps.braid(maps.braid(t,3),3)))
                    checked+=1
        self.assertGreater(checked,100)

    def test_known_positive_double_transposition_representatives_retained(self):
        root=Path(__file__).resolve().parents[2]
        catalog=json.loads((root/'exports/latest-discoveries-20261023-v12/release-v12/all-smooth.json').read_text())
        checked=0
        for record in catalog['components']:
            if record['genus']<2 or sorted(x for x in record['cycle_types'][0]if x>1)!=[2,2]:continue
            n=record['degree'];t=[]
            for text in record['tuple']:
                a=list(range(n))
                for raw in re.findall(r'\(([^()]*)\)',text):
                    row=[int(x)-1 for x in raw.split()]
                    for x,y in zip(row,row[1:]+row[:1]):a[x]=y
                t.append(tuple(a))
            t=tuple(t);edges=[row for row in maps.cycles(t[0])if len(row)>1]
            self.assertTrue(double.first_twists(edges,maps.compose(t[0],t[1]),t[2],t[3]))
            self.assertIsNone(double.middle_word(t));checked+=1
        self.assertGreater(checked,0)


if __name__=='__main__':unittest.main()
