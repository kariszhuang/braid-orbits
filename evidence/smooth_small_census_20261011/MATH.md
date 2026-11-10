# Exact scope and elementary implication

Product-one and transitivity are essential hypotheses. Human mathematical
review is pending. No sampling is used in this finite target.

## Smooth implies zero defect without topology

Let abcd=1 on n points. Put C=cycles(ab), L=orbits(<a,b>),
R=orbits(<c,d>). Joining each permutation cycle by a tree gives a graph
with ind(a)+ind(b) edges, n vertices and L components. Consequently
ind(a)+ind(b)>=n-L. Also L<=C because ab belongs to <a,b>.

If ind(ab)=ind(a)+ind(b), then n-C>=n-L, so L=C. The bipartite
incidence graph with the L left and R right pair-orbits as vertices has
one edge per ab-cycle. Each cycle lies in both corresponding pair-orbits
because ab=(cd)^-1. Whole-tuple transitivity makes this graph connected:
a union of graph components would be preserved by all four permutations.
Thus C>=L+R-1, forcing R=1. Defect C-L-R+1=0. Exchanging the two halves
proves the other smooth equality. This holds for every degree and every
transitive generated group, including identity entries.

## Why the exhaustive reduction covers every smooth tuple

For a smooth tuple at least one half is index-additive. The rotation
(a,b,c,d)->(c,d,a,b) preserves product one, transitivity, genus and defect.
It is achieved modulo simultaneous conjugation by applying the positive
word sigma1 sigma2 sigma3 twice (moves listed chronologically). Indeed
one such word gives the cyclically shifted tuple conjugated by a.
Therefore every everywhere-smooth braid component has a state with an
index-additive first pair. It is enough to enumerate all such first pairs.

Every product p in S_n is conjugate to the canonical permutation of exactly
one integer partition of n. For each canonical p, iterate every a in S_n;
b=a^-1 p is forced. Test ind(a)+ind(b)=ind(p), with no other exclusions.
This literally enumerates every geodesic first pair up to simultaneous
conjugation. It does not assume that entries are 3-cycles or involutions,
or that the group is primitive, alternating or symmetric.

Canonical pair keys handle disconnected actions. On each pair orbit, root
at every point and label by forward breadth-first traversal with generators
a,b in this fixed order. Minimize the labeled image tuples over roots.
Sorting these component keys is a complete isomorphism invariant for the
ordered pair action; it neither identifies distinct nonisomorphic pairs
nor omits relabeling-equivalent pairs. Forward traversal reaches the orbit
since the generators are finite permutations.

For every retained pair representative, independently enumerate every
c in S_n, and set d=(abc)^-1. Transitivity is the only admission test for
the implication experiment. For the everywhere-smooth orbit search, exclude
tuples with any identity entry, recording those separately. Every admitted
nonidentity occurrence must resolve to a nonsmooth word, exact smooth
six-move closure, or an explicit cap/timeout. No cap proves nonexistence.

## Labeled counting and independent checks

For a connected ordered pair action, an equivariant automorphism is
uniquely determined by the image of a root; its size is the multiplicity
of the minimal rooted key. Equal component types may be permuted. Thus the
whole-pair centralizer has order product(A_type^m_type * m_type!), and its
simultaneous S_n orbit has n!/centralizer representatives. Tests compare
this formula and every canonical pair against literal S4 conjugation.

Multiplying each frame's transitive third-entry count by its pair-orbit
size counts every labeled product-one transitive tuple with geodesic first
half exactly once. Let this total be N_L and the corresponding total with
BOTH halves geodesic be N_B. Half rotation gives N_R=N_L. Hence ALL smooth
labeled transitive tuples number 2*N_L-N_B. These are quadruples, not
distinct braid components. Independent Python repeats all third entries,
checks every implication occurrence and every literal negative word, and
validates any smooth closure. Full literal S4 triples independently confirm
the reduction, inclusion-exclusion and zero-defect implication.

An integer-partition product has product(Catalan(cycle_length)) geodesic
factorizations. This is a redundant mass regression only; the constructor
still tests every a and never uses this formula to discard one. Pair-orbit
masses also reproduce the literal conjugacy-class-weighted pair totals.

## Why the genus-zero full-cycle components are smooth everywhere

This is a sufficient condition in every degree, not a universal converse.
If the covering genus is zero and d is an n-cycle, then
ind(a)+ind(b)+ind(c)=n-1=ind(d^-1)=ind(abc). The transposition-length triangle
inequality gives

    ind(abc) <= ind(ab)+ind(c) <= ind(a)+ind(b)+ind(c) = n-1.

Both inequalities must be equalities; hence ind(ab)=ind(a)+ind(b).
If the full cycle is in another position, cyclically rotate the product-one
tuple and apply the same argument to the half not containing it. Braid moves
preserve genus and the multiset of entry cycle types, so every braid translate
is smooth. These are rational-curve cases (over C all smooth genus-zero curves
are isomorphic to P1), distinct from the positive-genus existence question.
The finite census finds that all its identity-free everywhere-smooth components
in degrees5--8 have genus zero and contain a full cycle. That finite converse
is not asserted in other degrees.
