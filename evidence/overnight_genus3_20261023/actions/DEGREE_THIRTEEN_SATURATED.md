# Degree thirteen: only the unsaturated index2 budget remains open

October23,2026. Identity-free transitive product-one quadruples of source
genus3, full six-move braid orbit modulo simultaneous S13 conjugation.
This is NOT a complete degree13 exclusion or a global theorem.

## All scalar budgets and exact remaining scope

After C225 removes any transposition entry, the complete scalar list is

    (2,8,10,10), (2,9,9,10), (3,9,9,9), (6,6,6,12).

C209 excludes the entire last, triangle budget. The results below exclude
the first and third. Therefore the ONLY open degree13 budget is
(2,9,9,10), whose index2 entry may still be a single3-cycle OR2^2.
The single-cycle theorem requires TWO partners with the saturated cycle
count. It does NOT remove the single3-cycle in this middle budget, where
only one ordinary entry has three cycles. Do not make that extrapolation.

## Saturated normal forms

Let the special index be s=2 or3. Arrange the tuple(a,c,b,e) so that c and
e have s+1 cycles. Their indices are13-s-1 and the remaining b-index is
6+s. Smoothness forces p=ac to be a full13-cycle. Relabel p=(0 ...12).

If a is a single(s+1)-cycle, C209 would require(s+1) to divide13. Neither
3 nor4 does. Thus the active special types left are2^2 for s2, and3·2 or
2^3 for s3. For2^2 the active support has size4, two active cycles; C209
bounds each c/e cycle length between4 and6. Three cycles summing13 must
therefore have type5,4,4. For3·2 and2^3 the active supports are5 and6;
the lower bound is3, and four cycles summing13 must have type4,3,3,3.
This also forces the third ordinary entry's profile when s3, though the
constructor does not need to impose that extra restriction on b.

Enumerate the entire conjugacy class of a, form c=a^-1p, keep these exact
c types, and quotient ONLY by rotations preserving p. No free division
by a centralizer size is used. This gives the complete marked frames below.

## Full cyclic forest constraints and forced complements

For each frame enumerate every unordered partition into e-cycle supports
of the forced sizes. For each shift k, the incidence graph of the active
a_k cycles and prospective e-cycles must be a tree. An equivalent test
expands each a_k cycle into star transpositions and verifies that its s
edges join distinct components of the e-cycle labels, with no cycle.
This depends only on supports, not on orientations of e's cycles.

Keep supports satisfying this condition for ALL13 shifts, then enumerate
ALL cyclic orientations. Set b=p^-1 e^-1. If ind(b)!=6+s, it is outside
the budget. Otherwise every shift must satisfy ind(a_k b)=6+2s. Any failure
has the literal negative word Q1^(2k),Q2^-1, length at most25. The opposing
pair c,e has sum greater than12, so no complementary smoothness equality
can rescue that endpoint. There are zero survivors, not capped survivors.

| Active type | All labelled a | Full-p frames | Profile frames | Rotation frames | e-supports passing all shifts | Oriented e | Wrong b-index | Literal negatives |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
|2^2|2145|1430|65|5|117|101088|61880|39208|
|3·2|25740|6435|117|9|26|1248|832|416|
|2^3|25740|8580|156|12|26|1248|832|416|

Support partitions before filtering number45,045 in the first row and
200,200 in each other row. Oriented e counts are necessary marked labelled
choices, not distinct quotient states or braid components.

## Independent complete audit

The Rust producer constructs a by edge matching or explicit3-cycle/edge
choices and support partitions by subset masks. The separate Python
auditor enumerates the WHOLE active class in S4/S5/S6 first, embeds it on
every support subset of13 letters, and uses restricted-growth colour words
for e supports. It checks forests by explicit component sets, not the
producer's DSU. Cycle orientations use a different permutation generator.

Every one of the103,584 (a,e) records is reconstructed, matched with no
duplicates to the producer's complete binary set, and each of40,040
index-correct cases is independently replayed to a nonsmooth endpoint.
No BFS, cap, geometric isotriviality classification or assumed fullA13/S13
generation is needed. All unresolved counts are zero.

Sources and raw immutable hashes are in D13_SATURATED_CHECKS.json.
The remaining unsaturated budget is NOT represented by these normal forms;
it must retain its own independent search or proof.
