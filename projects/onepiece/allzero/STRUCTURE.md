# Structure of the primitive one-piece counterexamples (November 8, 2026)

## Projective planes: tuples of central collineations of PG(2,q)

Verified primitive orbits, all one-piece at every state (hence all-zero), all isotrivial,
none all-smooth. Files: COUNTEREXAMPLES.jsonl, GENUS3_COUNTEREXAMPLES.jsonl,
PSL37_COUNTEREXAMPLES.jsonl, COLLINEATION_COUNTEREXAMPLES_q457.jsonl (collin search, independently
re-verified with verify_cex.py).

| q | points | genus | states | generated group | entry orders | count |
|---|---|---|---|---|---|---|
| 4 | 21 | 1 | 28 | PSL(3,4), 20160 | 2,2,2 (elations) + 7 (fixed-point-free) | 5 |
| 5 | 31 | 1 | 12 | PSL(3,5), 372000 | 5,2,2,4 / 5,4,2,2 | 2 |
| 5 | 31 | 3 | 4 | PSL(3,5) | 2,4,4,4 (homologies) | 2 |
| 7 | 57 | 1 | 4 | PSL(3,7), 1876896 | 2,2,2 (homologies) + 7 (elation) | 2 |
| 7 | 57 | 4 | 12 / 4 | PGL(3,7), 5630688 | 3,2,2,6 / 3,2,3,3 | 4 |

Geometry (deep_psl.py):
* Each entry is a central collineation (homology: fixes an axis pointwise plus a centre;
  elation: fixes an axis pointwise, centre on the axis), except the order-7 element for q = 4.
* Centres and axes are in general position (no three centres collinear, no three axes
  concurrent, no centre on another entry's axis; one exception for q = 4).
* Mechanism at a cusp (checked at every state for q = 5, 7): the two entries of a side both fix
  the point P = (axis1 meet axis2) and the line L through the two centres. The orbits of the side
  group are {P}, the points of L (or of L minus P when P lies on L), and the rest of the plane
  (for q = 5: 1 + 6 + 24 or 1 + 5 + 25; for two involutions the group is dihedral and the orbits
  are finer). The genus sits in the big orbit; the point and line pieces are rational.
* Duality: the same matrices acting on LINES give a cover of the same genus whose orbit is also
  one-piece at every state.

Proof outline to pursue: for central collineations in general position, the pair groups preserve
a point-line pair, the pieces are {P}, the line, and the complement, and cycle counting of
p = g1 g2 (which fixes P, L and the corresponding pair of the other side) gives defect 0 and a
single positive-genus piece at every cusp. Not yet a proof.

## Product action on 5 x 5

Two entries swap the coordinates, so the cover factors through a double cover Y -> P1 branched at
those two points (Y of genus 0). Over Y the row and column actions are Galois-conjugate degree-5
covers R, C of Y: genus 0 with group S5 (60-state orbits) or genus 1 with group A5 (12-state
orbits). R x_Y C has genus 6, and X = (R x_Y C)/iota (iota swaps the factors) has genus 1.

## Galois closures (Riemann-Hurwitz: 2g-2 = |G|(-2 + sum(1 - 1/ord)))

| q | |G| | entry orders | genus of Galois closure |
|---|---|---|---|
| 4 | 20160 | (2, 2, 2, 7) | 3601 |
| 5 | 372000 | (2, 4, 4, 4) | 139501 |
| 5 | 372000 | (5, 2, 2, 4) | 102301 |
| 7 | 1876896 | (2, 2, 2, 7) | 335161 |
| 7 | 5630688 | (3, 2, 2, 6) | 1407673 |
| 7 | 5630688 | (3, 2, 3, 3) | 1407673 |

## Nov 8 (late): more planes

* q = 8 (collin, full central-collineation search, 73 points): exactly 2 orbits, genus 3, 4 states,
  PSL(3,8) (order 16482816), three elations of order 2 + one homology of order 7, isotrivial
  (COLLINEATION_COUNTEREXAMPLES_q8.jsonl, independently verified).
* q = 11 (collin_big, passport filter): (2,2,2,11) gives NO all-zero orbit (2.6M tuples, none with
  defect 0 at the first cusp); (2,2,2,any) gives 4 orbits of genus 6, 4 states, PSL(3,11)
  (order 212427600), three involution homologies + one homology of order 5, isotrivial, one-piece
  everywhere (COLLINEATION_COUNTEREXAMPLES_q11.jsonl, independently verified).
* Genus for three involution homologies + a homology of order k (odd q):
  2g - 2 = (q^2 - 1)(5/2 - 1/k) - 2(q^2 + q + 1); q = 11, k = 5 gives g = 6.
* Genus grows with q: 1 (q=4), 3 (q=5), 4 (q=7), 3 (q=8), 6 (q=11). All isotrivial, none all-smooth.
* q = 13 (collin_big, (2,2,2,any), 183 points): 4 orbits of genus 7, 4 states, PSL(3,13) (order
  270178272), three involution homologies + one homology of order 4, isotrivial, one-piece
  everywhere (COLLINEATION_COUNTEREXAMPLES_q13.jsonl, independently verified). Formula check:
  (169-1)(5/2-1/4) - 2*183 = 12 = 2g-2.
* "Three involutions + one central collineation" so far: q=7 elation of order 7 (g 1); q=8 homology
  of order 7 = q-1 (g 3); q=11 homology of order 5 = (q-1)/2 (g 6); q=13 homology of order 4
  = (q-1)/3 (g 7). Which fourth entry occurs is not yet explained.
* q = 9 (collin, full search, 91 points): 2 orbits of genus 5, 4 states, PSL(3,9) (order 42456960),
  three involution homologies + a homology of order 8 = q-1, isotrivial
  (COLLINEATION_COUNTEREXAMPLES_q9.jsonl, independently verified). (80)(5/2-1/8) - 182 = 8.
* Determinant constraint: the four matrices multiply to a scalar c I, so for three involution
  homologies (det -1 each) and a homology with multiplier a: a = -c^3, i.e. a is minus a
  cube. q = 7: cubes {+-1}, so a = 1 (an elation) or -1 -- found the elation. q = 13: cubes =
  order-4 subgroup, so a has order dividing 4 -- found order 4. q = 8, 9, 11: every element is a
  cube, no constraint; found orders q-1, q-1, 5. A further geometric condition must choose these.

## Nov 9 correction

q = 4 has 2 distinct orbits (not 5) up to S_21 relabelling; the plane family totals 24 orbits
(2, 4, 6, 2, 2, 4, 4 for q = 4, 5, 7, 8, 9, 11, 13), and with the 5x5 grid 29 primitive one-piece
orbits. Checked by evidence/onepiece_publication_20261108/verify.py.
