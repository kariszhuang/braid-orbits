# Full-cycle twists: stronger elementary obstructions

October22,2026. Source genus is denoted by g, degree by d. Quadruples have
product1, no identity entry, and transitive monodromy. States are modulo
simultaneous S_d conjugation. Every assertion about an all-smooth component
means the full six-move braid orbit, not a sample or a pure-braid suborbit.
The proof arguments below are elementary drafts for mathematical review;
the finite certificates have separate literal independent checks. No
literature-priority or global nonexistence claim is made. C209/CHECKS.json
binds the computations and their precise populations.

Write Q1(a,b,c,e)=(aba^-1,a,c,e). Its square simultaneously conjugates a,b
by p=ab, leaving c,e unchanged. Thus Q1^(2k) replaces a by p^k a p^-k.
The inverse middle move Q2^-1 then puts (p^k a p^-k,c) first. Entry indices
are only permuted by braids. This simple action is the common mechanism
behind the following results.

## 1. The fixed-point obstruction is not restricted to involutions

**Lemma.** Suppose p=ab is a full d-cycle,

    ind(a)+ind(c)=d-1,       ind(b)+ind(e)>d-1,

and a and c each fix a point. Then (a,b,c,e) has a nonsmooth braid translate.
There is a literal witness of length at most d for odd d, or d+1 in general.

**Proof.** Choose fixed points x of a and y of c. Some k modulo d satisfies
p^k(x)=y. After Q1^(2k), a_k=p^k a p^-k and c both fix y. Apply Q2^-1.
The first product a_k c fixes y, so its index is at most d-2. This cannot
equal the first pair sum d-1 or the complementary sum greater than d-1.
Choose the shorter positive/negative outer twist, giving length
2 min(k,d-k)+1. No assumption of involutions, unique fixed points, full S_d
generation, blocks, or geometry is needed. QED.

**Corollary (an infinite index family).** Let d be odd and g>0. A quadruple
whose indices, in some order, are

    ((d-1)/2, (d-1)/2, (d-1)/2, (d+4g-1)/2)

cannot be all-smooth, whenever these are valid indices.

**Proof.** Put the last index fourth. The complementary pair sum is
d+2g-1>d-1, so smoothness forces ab to be a full cycle. A permutation of
index r moves at most2r points: each nontrivial cycle of length l contributes
l-1 and satisfies l<=2(l-1). Each of the first three entries therefore
fixes a point. The lemma applies. QED.

For source genus3 this excludes the ENTIRE degree13 triangle budget
(6,6,6,12), not merely the passport (2^6 1)^3(13). It also excludes
(7,7,7,13) in degree15 and the analogous valid odd-degree budgets. The
whole degree13 is NOT excluded: its star budgets remain open. The1000
sampled degree13 seeds replayed by the new independent routine include nine
different non-involution first-entry types. Those are regressions of the
formula, not exhaustive type coverage or the proof of the corollary.

This improves the general triangle-regime lower bound to d>=4g+2. Indeed,
the triangle inequality b+c<=d-1 and Riemann--Hurwitz force a>=2g and
d>=4g+1. Equality would force indices(2g,2g,2g,4g), precisely the excluded
half-rank family. In source genus3, a triangle target must therefore have
degree at least14, not13.

## 2. A single cycle forces a balanced partition

**Theorem.** Suppose g>0 and an all-smooth quadruple has an entry a consisting
of one nontrivial k-cycle, and two other entries b,c with exactly k cycles
each, fixed points included. Then k divides d, and b and c both have uniform
cycle type (d/k)^k.

**Proof.** Arrange these entries as (a,b,c,e) using braids. Their indices are

    k-1, d-k, d-k, 2g+k-1.

The last index follows from Riemann--Hurwitz. Hence the complementary pair
sum ind(c)+ind(e)=d+2g-1 exceeds d-1. Smoothness forces ab to be a full
d-cycle p.

For every k' modulo d, twist the first pair and apply Q2^-1. The new first
pair is (p^k' a p^-k',c). Its complementary index sum again exceeds d-1,
so its product must be a full cycle.

Multiplication by a k-cycle increases the index by k-1 exactly when its k
support points belong to k DISTINCT cycles of the partner permutation.
Indeed, expand the k-cycle into k-1 star transpositions. Every step must
join, not split, two cycles; repeated initial cycle labels force a split.
Thus the shifted support S_k'=p^k'(supp a) meets each cycle C of c once.

Sum |S_k' intersect C| over all d shifts. Since p is a full cycle, each
point of C occurs in exactly k shifted supports. Consequently

    d = k |C|.

Every cycle of c has length d/k. In particular k divides d. Interchange b,c
by braids, retaining their cycle counts and the single-cycle type of a, and
repeat to obtain the same conclusion for b. QED.

For k=2 this is stronger than an odd-cycle-colouring argument: a transposition
and two two-cycle entries force both latter types to be (d/2)^2, in every
positive genus. In odd degree the configuration is impossible outright.

## 3. More generally: cycle-support bounds without Hodge theory

Let a have index s, active support size S, and q nontrivial cycles. Thus
S=s+q. Suppose an all-smooth positive-genus quadruple has two other entries
b,c each with exactly s+1 cycles. Then each cycle C of b or c satisfies

    ceil(d/S) <= |C| <= floor(q d/S).                 (3.1)

To prove this, order the tuple (a,b,c,e). The fourth index is2g+s, so the
same complementary-sum argument makes ab and every shifted a_k c full cycles.
Contract cycles of c in a minimal transposition expansion of a_k. Additivity
means the resulting s-edge graph on s+1 vertices is a tree. Equivalently,
the incidence graph of q active a_k-cycles and s+1 c-cycles is a bipartite
tree with S edges. Every c-cycle is hit at least once. A c-cycle is hit at
most once by each active a_k-cycle; two hits would create a cycle in the
incidence graph. Therefore

    1 <= |supp(a_k) intersect C| <= q

at every shift. Summation gives d<=S|C|<=qd, proving(3.1).
The single-cycle theorem is exactly q=1, not an unrelated special trick.

This replaces a weak forest-start condition by a whole cyclic family of
forest constraints. It is still only necessary: it does not certify the
other braid moves or a full all-smooth orbit.

## 4. A complete negative source-genus-three result in degree11

For sorted indices a<=b<=c<=e, positive-genus-three Riemann--Hurwitz gives
a+b+c+e=2d+4. Each of the three class pairings occurs in the full braid orbit.
At least one pair sum in each pairing must be <=d-1. Equivalently,

    a+c<=d-1,       min(a+e,b+c)<=d-1.

In degree11 the only possible budgets are

    (1,7,9,9), (1,8,8,9), (2,8,8,8).

The first is excluded by the single-cycle theorem with k=2: two entries have
two cycles, but2 does not divide11. In the third budget the index2 entry is
either a3-cycle or a double transposition. The3-cycle case is excluded by
the same theorem with k=3. The other two cases have exact finite certificates:

| Budget / special type | First-pair frames | Constructed thirds | Qualifying states | Pending |
|---|---:|---:|---:|---:|
| (1,8,8,9), transposition | 5 | 190080 | 108339 | 0 |
| (2,8,8,8), double transposition | 60 | 13200 | 7051 | 0 |

Here these are necessary S11 states in a chosen ordered full-product slice,
after the cyclic forest condition, not braid components or all labelled
quadruples. All115390 qualifying states have literal negative words of length
at most12. Independent native replay and independent Python reconstruction
of EVERY qualifying state both pass.

### Why this is complete

In the transposition case, arrange indices(1,9,8,8). Smoothness forces the
first product p to be full. Normalize p(i)=i+1 modulo11. First-pair frames
are t=(0 r), 1<=r<=5, and b=t p. The first-pair centralizer is trivial.
For every twist, t_k must join two DIFFERENT cycles of the third entry c.
Thus its three-cycle partition is a proper3-colouring of the chord graph
whose edges join x to x+r. The producer grows restricted-growth colours;
the verifier independently partitions subsets and orders each cycle.

In the double-transposition case, the first product is again full. A double
transposition geodesic below p is a noncrossing matching on four circularly
ordered support points. There are2 binom(11,4)=660 labelled choices, or60
rotation classes. The independent verifier instead enumerates ALL three
matchings and selects geodesics by literal index. The centralizer is trivial:
a nonidentity rotation has odd orbit length dividing11 and cannot preserve
a four-point support. Each translated matching must form a tree on the
three cycles of the third entry. The two edges must join different colours
and collectively use all three colours. Only five frames admit any such
partition; each has33 partitions and2640 cyclic orderings. The other55
frames are excluded by complete partition enumeration, not skipped.

The fourth entry is calculated. Only those of the required index are in
scope; all are separately witnessed nonsmooth. Missing rows, false words,
duplicate rows, changed counts and omitted zero-survivor frames are rejected
by corruption tests. A surviving bounded local test would have been saved
as PENDING, not called negative. There are none.

**Conclusion.** There is no identity-free source-genus-three all-smooth braid
orbit in degree11, for ANY transitive permutation group. Together with the
general bound d>=2g+3 and the completed degree9/10 censuses, a target has
degree at least12.
This does not exclude genus0/1/2 or other genera in degree11.

## 5. Degree15: remove mixed profiles before enumerating completions

Consider the saturated star budget(4,10,10,10). Every ordinary entry has
five cycles. The possible active types of the index4 entry are

    5; 4 2; 3^2; 3 2^2; 2^4.

Their active support sizes are5,6,6,7,8. For the first four types,
ceil(15/S)=3. Formula(3.1) forces all five ordinary cycles to have length
at least3; their lengths sum to15. Thus ALL three ordinary entries must
have type3^5. This is an elementary exclusion of every mixed ordinary profile
in these four special types, not just the18 earlier Hodge-selected cases.

The pure4 2 case was already exhausted,5600 necessary S15 states. New exact
checks complete the other three multi-cycle pure cases:

| Active special type | First-pair frames | Completions/frame | Necessary S15 states | Pending |
|---|---:|---:|---:|---:|
| 3^2 | 1 | 2800 | 2800 | 0 |
| 3 2^2 | 9 | 2800 | 25200 | 0 |
| 2^4 | 6 | 2800 | 16800 | 0 |

All44800 have literal negative words. Each frame set is checked by exact
Frobenius connection coefficients, literal passports and distinct rotation
classes with trivial pair centralizer. The independently count-validated2800
completion library is REUSED, not regenerated. A full product establishes
transitivity; full A15 generation is neither required nor asserted.
The special5 pure case has the previously written isotriviality argument;
it is not a new negative passport or a new positive discovery here.

For the2^4 special, (3.1) permits ordinary lengths2 through7. The seven
possible five-cycle types are

    7 2^4; 6 3 2^3; 5 4 2^3; 5 3^2 2^2;
    4^2 3 2^2; 4 3^3 2; 3^5.

Its pure3^5 triple is now excluded, but mixed triples in this scope remain
open. The bound is not a declaration that any of them is all-smooth.

## 6. A much narrower next triangle target in degree14

The unique degree14 triangle index budget is(6,6,7,13): from a>=6 and
b+c<=13, Riemann--Hurwitz leaves no other ordered possibility.

The index7 entry MUST be fixed-point-free. Otherwise order it first, put
an index6 entry second and the other index6 entry third. Smoothness forces
the first product to be a full14-cycle. Both the first and third entries
have fixed points, so the lemma in Section1 gives a contradiction. An index7
fixed-point-free permutation on14 points necessarily has type2^7.

Write the tuple(c,a,b,e), with c of type2^7, a,b of index6 and e a full
14-cycle. Let A,B be the fixed-point sets of a,b. Both have at least two
points. Twist c,a by p=ca. The extra inverse moves Q1^-1,Q2^-1 then put
(p^k a p^-k,b) first. Its required additive index is12. Two common fixed
points would give at least three product cycles and index at most11.
Consequently

    |p^k(A) intersect B| <= 1 for EVERY k,
    |A| |B| <= 14.                                      (6.1)

The second line sums the first over14 shifts. At least one of |A|,|B| is
therefore at most3. The corresponding entry has type2^6 1^2 or3 2^4 1^3.
This is a useful passport reduction, not a completed degree14 search.

The fully involutory subcase

    (2^7), (2^6 1^2), (2^6 1^2), (14)

is an especially small proposed exact target. The first two involutions form
a connected alternating14-point path, so there is one first-pair frame up
to S14 conjugacy. It has a two-element centralizer (path reversal), unlike
the odd-degree frames used above. There are

    14! / (2^6 6! 2!) = 945945

third involutions before filtering the fourth entry. This population has
NOT been enumerated here. Counts modulo its residual centralizer require
care; labelled occurrences must not be misreported as distinct S14 states.

## 7. What this does and does not say about global nonexistence

We now have more proved exclusions, not a proof in every degree. Important
gaps include the unsaturated degree12 budget(1,9,9,9), whose additive first
product has TWO cycles, and the degree13 star budgets. Neither is covered by
the full-cycle rotation argument. Higher-degree triangle budgets can also
escape the half-rank fixed-point corollary.

The productive next construction is to impose the cyclic forest constraints
on partitions BEFORE choosing cycle orientations and calculating the fourth
entry. This is what the new small exact enumerators do. It tests all outer
twists at once, unlike a handful of local random braid tests, and imposes no
block system. It is still followed by literal rejection or full closure.

For the continuous search, retire genus3 degree11 and replace the now-closed
degree13 triangle focus with degree12(1,9,9,9). Keep generic higher-degree
and different quotient mechanisms. The scalar and elementary reductions
do not classify isotriviality or change public data. No new GPU benefit is
claimed: these small irregular exact enumerations are already fast on one
CPU; RAM is used only for small frame/partition/word libraries.
