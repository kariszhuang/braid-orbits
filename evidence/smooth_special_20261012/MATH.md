# Elementary smoothness reductions and a transversal-block theorem

Original derivations for this mission; human review pending, no novelty claim.
Permutation composition acts rightmost first. Fixed points count as cycles.

## 1. Cycle-incidence forest criterion (exact, all permutations)

For a,b in S_n form a bipartite multigraph: vertices are the cycles of a and
the cycles of b; each point is an edge between its two containing cycles.
Connected components are exactly the orbits of <a,b>. Write L for their
number and C for cycles(ab). Always L<=C. Its cycle rank is

    beta = n-cycles(a)-cycles(b)+L.

If indices add, cycles(a)+cycles(b)=n+C, hence beta=L-C<=0. Nonnegative
cycle rank forces beta=0 and L=C: the graph is a forest.
Conversely a tree has one boundary cycle under the cyclic orders induced by
a,b: successively remove leaf vertices (a fixed point of one permutation),
contract their edge and insert/delete that point in the other permutation's
cycle. This preserves the single face until the one-edge tree is reached.
Thus a forest has C=L, and beta=0 gives ind(a)+ind(b)=ind(ab).
This is an exact O(n) DSU test; no genus or transitivity hypothesis required.
For a quadruple, smoothness is equivalent to at least one complementary
pair-cycle incidence graph being a forest. Computing product cycles remains
a separate equally-linear control, not an asserted universal timing gain.

## 2. Exact slack and cheap passport obstruction

Index is transposition length, so ind(ab)<=ind(a)+ind(b). Product one gives
ind(ab)=ind(cd). Therefore smoothness is exactly

    ind(ab)=min(ind(a)+ind(b),ind(c)+ind(d)).

The nonnegative slack between these quantities is even by permutation parity.
For indices i1<=i2<=i3<=i4, an everywhere-smooth orbit necessarily satisfies

    min(i1+i4,i2+i3)<=n-1.

Proof: every ordering of the four conjugacy classes is braid-realizable.
For every complementary partition, some pair sum must be <=n-1. Of the three
minimum pair sums, the largest is min(i1+i4,i2+i3). This filter alone cannot
certify existence or smoothness. With a full n-cycle entry and genus g>0,
put it next to each other entry: its pair is never geodesic (identities
excluded), so the opposite pair sum n-1+2g-i must be <=n-1. Hence every other
entry has index at least2g. This excludes many positive-genus full-cycle
passports before constructing permutations, not all of them.

## 3. Commuting-pair and abelian exclusions

If a,b commute, their action on any joint orbit of size m is a regular
abelian quotient. Let their local orders be s,t and the product order u.
The local slack equals m*(1-1/s-1/t+1/u). If both act nontrivially,
s,t>=2 makes this strictly positive. If either acts trivially, it is zero.
Thus commuting permutations are index-additive exactly when their moved
supports are disjoint. A transitive abelian permutation group is regular,
so each nonidentity element is fixed-point-free. Consequently no identity-free
tuple generating a transitive abelian group is even smooth at one state.
Do not mistakenly exclude nonabelian imprimitive groups.

## 4. New braid-stable sufficient condition

Partition n=rk points into r blocks of size k. Assume all four entries
preserve these blocks and induce translations in the cyclic group C_r.
Exactly one entry a acts nontrivially inside ONE block and trivially outside.
The other three induce unit translations mod r and have cycle type r^k.
Assume product one and transitivity.

Then the entire braid orbit is smooth, even if a has several cycles or order
different from three. Reason: all entry conjugates preserve the blocks.
The induced quotient C_r is abelian, so braid moves merely swap quotient
translations. Conjugation keeps a supported inside one block. A unit
translation entry of type r^k has each cycle meet every block exactly once.
Thus the points in the support block lie in distinct cycles of that entry.
Each cycle of a joins those distinct cycles, and different cycles of a join
disjoint sets. Equivalently its cycle-incidence graph with that entry is a
forest. Their indices add. At every state, whichever half contains a is
geodesic, hence the tuple is smooth. Defect is zero by the earlier theorem.

If J=ind(a), the covering genus is

    g=1 + (k*(r-3)+J)/2.

Thus r=3 allows new small types2^2 (J2),5 (J4),3^2 (J4),7 (J6), etc;
r=5 or7 gives higher genera even with a single small3-cycle. This is not
the earlier genus-zero/full-n-cycle mechanism. It generalizes the known
three-block commutator mechanism, not proof of an unrelated construction.

Even r is impossible under these hypotheses: each of the three unit
translations is odd modulo2, but their sum must be0 modulo r. Likewise r
must divide the sum of the translations. These are quotient obstructions,
not impossibility results for other block actions.

## 5. Explicit, automatically transitive construction

Let r>=3 be odd, choose q,q+1 units mod r, and x=(1..k), y any permutation.
Give c block translation q and within-block voltages

    v0=x, v_q=y, v_(-q)=(yx)^-1, other v_j=1.

The return product along c from block0 is (yx)^-1*y*x=1; hence c has type r^k.
Let b be the untwisted translation1. Its product b*c has translation q+1.
Starting at block0, let V be its within-block return map after r steps.
Take a=V^-1 on block0 and identity elsewhere, and d=(abc)^-1. The return
map of abc is a*V=1; hence abc and d have type r^k. If a is nonidentity the
preceding theorem applies. In the q=1 case, V=y*(yx)^-1*x, so

    a=x^-1*y*x*y^-1.

There is no group-generation bottleneck for transitivity: b^-q*c preserves
every block and restricts to the full k-cycle x on block0. Together with b
it makes all rk points lie in one orbit, independently of y or a.
Exact group order and full ambient generation remain separately measured.
The preserved nontrivial block system rules out full A_n/S_n here.

Isotriviality does not follow merely from smoothness, fixed genus, or finite
quotient closure. The initial construction left geometry untested; section7
below supplies a separate post-run geometric argument using the block quotient.

## 6. Genus-zero smoothness is exactly complementary-pair transitivity

Additional post-run interpretation. Put C=cycles(ab), L=orbits(<a,b>),
R=orbits(<c,d>). For a transitive product-one genus-zero tuple, the genus
decomposition forces both pair genera and defect to vanish. Consequently
C=L+R-1. Index additivity on the left is equivalent to C=L, hence R=1;
on the right it is equivalent to C=R, hence L=1. Therefore

    genus0: smooth iff at least one complementary pair is transitive.

This does not imply genus zero is smooth everywhere. It explains why a full
cycle suffices, and why genus-zero weak-failure components are nonsmooth.

## 7. Post-run geometric distinction: quotient genus, not source genus

This section supersedes the initially untested geometric status above and
the frozen launch plan. It is a mathematical interpretation requiring human
review, not a computer computation of curve equations or moduli invariants.
Work over C. Normalize the three branch points with nontrivial block action
to0,1,infinity; let the fourth point t move. After a finite marking/base
change, the block quotient is a fixed cyclic degree-r cover

    E -> P^1,

branched fully at those three fixed points. Its genus is (r-1)/2 by
Riemann-Hurwitz. The original cover factors as

    C_t -> E -> P^1,

with degree k on the first map. The r^k cycle types have exactly the same
ramification index r as the quotient at0,1,infinity, so C_t->E is unramified
there. The local entry has identity block action and nontrivial action only
in one block, so C_t->E is ramified over exactly ONE point P_t above t,
with local ramification type equal to the small entry. P_t moves nonconstantly
on the fixed E. This also rederives g(C_t)=1+k*(r-3)/2+ind(a)/2.

For r=3, E is elliptic. Translate P_t to a fixed origin. Degree-k covers of
that fixed once-punctured elliptic curve with the specified monodromy form
a finite set (finite permutation representations of its finitely generated
fundamental group). On a connected family this discrete covering class is
constant after finite base change. The source curve is therefore isotrivial.
This argument applies even when the source has genus3 or4: the intermediate
curve's genus is the relevant distinction.

For r>=5, E has genus>=2. If the family were isotrivial, trivialize its
source after a finite base change, obtaining maps from a fixed smooth C to
the fixed E while the branch point P_t still moves. The de Franchis theorem
says there are only finitely many nonconstant maps C->E. Their branch loci
are a finite set, contradicting nonconstant P_t. Thus these families are
NON-ISOTRIVIAL. This conclusion uses the actual invariant block quotient
and localized ramification, NOT merely the genus/state count of the source.

Primary source for de Franchis finiteness: Alan Howard and Andrew J. Sommese,
On the theorem of de Franchis, Ann. Scuola Norm. Sup. Pisa (4),10(3),1983,
429-436; introduction on429 and Theorem1 on433. Source opened October12:
https://www.numdam.org/item/ASNSP_1983_4_10_3_429_0.pdf

Historical context, NOT a novelty claim: Jason Starr's August5,2011 15:49
comment already describes using a fixed higher-genus Belyi curve to produce
varying source curves in complete four-branch families:
https://mathoverflow.net/questions/71898/curves-which-are-not-covers-of-p1-with-four-branch-points
The exact permutation tuples and state counts here are computed separately;
do not call this geometric principle new to the literature.
