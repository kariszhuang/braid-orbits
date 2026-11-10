# Complete degree-thirteen source-genus-three exclusion

October23,2026. Identity-free transitive product-one quadruples, full six-move
Hurwitz orbit modulo simultaneousS13 conjugation. This is a finite-degree
result, NOT global nonexistence or a theorem about all genera in degree13.

## Scalar coverage

The sorted indices sum30. Each class pairing must have one pair sum<=12.
After the source-reviewed C225 transposition obstruction removes index1,
the complete necessary list is

    (2,8,10,10), (2,9,9,10), (3,9,9,9), (6,6,6,12).

C209 excludes the last triangle budget in every cycle type. C229 completes
the first and third via saturated normal forms and literal negatives; see
../overnight_genus3_20261023/actions/DEGREE_THIRTEEN_SATURATED.md and its
D13_SATURATED_CHECKS.json. The remaining budget is(2,9,9,10). Its special
entry has either one3-cycle or two disjoint transpositions. No fullA13/S13
generation assumption is needed anywhere.

Arrange the tuple(a,c,b,e) with indices(2,10,9,9). Smoothness forces p=ac
to be a full13-cycle, since its required additive index is12 and the
complementary pair sum18 cannot be additive. Relabel p=(0 1 ...12).
For Q1(a,c)=(aca^-1,a), Q1^(2k) replaces a by a_k=p^k a p^-k. Following
by Q2^-1 tests a_k*b. Bringing e to position2 with Q3^-1,Q2^-1 instead
tests a_k*e, retaining a complementary c' conjugate to c. These are actual
four-point braids, not independent rotations of product components.

## One3-cycle: congruence before enumeration

Every a_k/e pair must be additive: its index sum11 is the only eligible
side. Thus the three support points of a_k lie in distinct e-cycles.
If C is any e-cycle of length ell, summing over all13 shifts gives
3ell<=13, hence ell<=4. Since3ell<13, some shift misses C entirely.

At that frame a joins exactly three of the four e-cycles into ONE full
product cycle M of length13-ell, leaving C untouched. All three support
points of a lie in M. Twist the (a,e) pair around ae, then pair its rotated
a with the complementary c', which still has exactly three cycles. This
new first pair must be additive of index12 and therefore transitive. Its
three a support points hit each c'-cycle exactly once.

Averaging those hits over the full rotation of M gives

    |M| = 3*|M intersect C'_j|   for every c'-cycle C'_j.

Thus3 divides13-ell. Every e-cycle has ell<=4 and ell congruent1 modulo3,
so e has type4^3*1. At a frame missing its fixed singleton, M has length12;
each c'-cycle contains four points of M. Distributing the remaining point
forces c' and hence c to have type5,4,4.

Since c=a^-1p, a cuts p into cyclic arcs of lengths4,4,5. Up to a rotation
there is one marked frame a=(0 4 8), c=a^-1p. The orientation of a must be
the splitting orientation. Every e-cycle support must have no two points
whose difference equals a difference within{0,4,8}. In the step4 cyclic
ordering this is proper colouring of the square of C13. Choose e's fixed
point (13 choices). The remaining12 points have three colours forced
periodically, up to colour names; each colour class has size4. Orient each
4-cycle freely, giving13*(3!)^3=2808 necessary labelled e choices.

For each, b=p^-1e^-1 is forced.1183 have the wrong b-index. All1625 index9
choices have a literal nonsmooth word of length<=3. An independent Rust
audit enumerates ALL75075 unordered supports of type4,4,4,1 (full labelled
permutation class16216200), filters the13 compatible supports and every
orientation, and replays every negative. Its sorted2808-permutation byte
set exactly equals the Python producer's, SHA256
ac7938e1c7db978c28145dcf4afe2f84766c176e35a82ff3221563e69f34e21c.
There are no unresolved cases.

## Two transpositions: bounded supports and complete negative certificates

In each a_k/e forest, any e-cycle can meet at most one endpoint from each
a edge. Summing the maximum two hits over all13 shifts gives4ell<=26,
so every e-cycle has length<=6. There are four e-cycles, so at least one
has length ell<=3. Its average intersection with the four-point support
is4ell/13<1; some shift misses it.

The two forest edges then join the other THREE e-cycle vertices in a tree.
Thus ae has one full active cycle M of length13-ell>=10, and the untouched
cycle. All four support points of a are in M. Twisting ae and testing
against c' must join all three c'-cycles, so each is hit at least once at
every rotation. Averaging gives

    4*|C'_j intersect M| >= |M| >=10.

Consequently every c'-cycle and hence every c-cycle has length>=3. The
allowed c types are7,3,3;6,4,3;5,5,3;5,4,4. Enumerating all2145 double
transpositions a,1430 have the required full-product cut count,325 meet
these c-length bounds, and quotienting by rotations gives25 marked frames.

For each frame enumerate every partition into FOUR e-cycle supports of
length1..6. There are2341339 such unordered partitions. Retain those for
which both rotated a-edges form a forest on these support labels for ALL13
shifts. Across25 frames72579 support partitions pass. Every possible cycle
orientation is included, giving6435312 (a,e) records. Force b=p^-1e^-1.

-3343470 choices have ind(b)!=9.
-3090763 index-correct choices have an explicit negative Q1^(2k),Q2^-1.
-The remaining1079 ALL have independently replayed negative words<=4 moves.

No cap or unclosed orbit is counted negative. The standalone Python
full-coverage auditor independently enumerates edge pairs and restricted-
growth four-colour words rather than the Rust matchings/subset partitions,
reconstructs every record, proves exact equality of the6435312-record set,
checks no duplicates, and replays EVERY negative. It finishes with zero
pending cases; sorted-record SHA256
2b04487f0c4931d478b6ed0649079dc5bed10cb68758ae36cfe83e7a94b24ad4.

## The failed audit was investigated, not waived

The earlier exact-set assertion correctly failed because the two programs
chose different canonical LABELS for the same25 cyclic-rotation classes:
Rust minimized the permutation array; the auditor minimized the edge list.
An isolated0.01s regression reproduced the mismatch. Independent orbit
normalization proved the class sets agree. The audit alone was corrected to
compare the same permutation-array marking, retaining its different edge-
and colour-based construction. The full143-second-scale audit was then
rerun and passed; no producer candidates were deleted or results overwritten.
The frozen producer remained unchanged. Regression source and receipt hashes
are bound in FULL_D13_CHECKS.json.

## Conclusion

Every scalar budget and every special cycle type in degree13 is now covered.
There is no identity-free transitive source-genus3 all-smooth full braid
orbit in degree13. Combined with inherited9--12 results and the general
d>=2g+3 bound, a target must have d>=14. This does not exclude all degree14
or unbounded degrees, and supplies no non-isotrivial positive example.

The new congruence and double-transposition coverage reductions were
independently adversarially reviewed by the resumed geometry agent.
Human or formal proof verification is not recorded.
