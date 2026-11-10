# A global dihedral obstruction and the degree-fourteen triangle

October 22, 2026. C212. Mathematical proof draft; independent human proof
review is unrecorded. No global genus-three nonexistence theorem is claimed.
The existing C209 six-worker search, sealed inputs and public labels are
unchanged. No completion enumeration or braid-orbit search was run here.

Permutations act rightmost first. Fixed points count as cycles. A state is
modulo simultaneous S_d conjugation. All-smooth means the FULL Hurwitz orbit.
Write i(x)=d-cycles(x). An additive pair has a cycle-incidence FOREST;
in particular two cycles, one from each entry, cannot share two letters.

## 1. A degree-independent exclusion

**Theorem 1.** A positive-genus all-smooth orbit cannot contain a state
whose first two entries are involutions with a full-cycle product.

Any selected pair of entries can be brought to the front by inverse
Hurwitz moves without changing those two entries, so this also excludes
such a pair in nonadjacent positions.
There is no alternating/symmetric generation, primitivity, block or
geometric hypothesis. The result concerns every positive source genus.

### The index argument makes both cross pairs additive

Suppose the state is (a,b,c,e), with a,b involutions and p=ab a d-cycle.
If necessary swap the last two branch classes using Q3^(-1), leaving the
first pair unchanged, so r=i(c)<=s=i(e). The reflection normal form below
gives i(a)+i(b)=d-1, and Riemann--Hurwitz consequently gives

    r+s=d+2g-1.                                             (1.1)

Normalize p(x)=x+1 on Z/dZ. The identity a p a=p^(-1) shows that

    a(x)=u-x,       b(x)=u-1-x.                              (1.2)

If d=2n+1, both involutions have index n. If r>n, both pair sums in the
matching (a,c)|(b,e) exceed d-1, contradicting smoothness. Thus r<=n.
Equation(1.1) then makes both complementary sums i(a)+s,i(b)+s strictly
greater than d-1. Every smooth cross pair involving c must therefore be
additive on the c side.

If d=2n, the indices of a,b are n-1,n, in some order. Suppose r>=n. Pair
the index-n involution with c and the index-(n-1) involution with e. The
first sum is at least2n. Since s>=r, equation(1.1) gives s>=n+g. The
second sum is at least(n-1)+(n+g)=2n+g-1>d-1. Hence neither side can be
additive. Thus
r<=n-1. Now s>=n+2g, and again both complements of the c cross pairs
have index sum>d-1.

The conclusion in either parity is that for EVERY integer k,

    (p^k a p^(-k),c) and (p^k b p^(-k),c) must be additive. (1.3)

Indeed Q1^(2k) conjugates a,b by p^k, leaving c,e fixed. Q2^(-1) brings
the first cross pair to the front; Q1^(-1) followed by Q2^(-1) brings
the second one. The complementary indices just calculated remain fixed.

### Reflections cover every pair of letters

The conjugates in(1.3) are reflections with constants

    u+2k and u-1+2k modulo d.

For odd d, the first set alone contains every constant. For even d, the
two sets cover the two parities. Consequently for ANY distinct letters
x,y, one of these reflections has the transposition(x y).

Choose x,y in the same nontrivial cycle of c. That cycle and the displayed
reflection transposition share two letters. Their incidence graph has
parallel edges, contradicting the forest required by(1.3).

Finally c cannot be identity: then e=p^(-1), the total index is2d-2 and
g=0. This contradicts the assumed positive genus and proves the theorem.

This is stronger than a fixed-point alignment argument: it aligns a
TRANSPOSITION with two letters of an external cycle, not merely one fixed
point with another. It excludes a whole mechanism in every degree.

### The full transitive-dihedral exclusion

**Corollary.** In a positive-genus all-smooth orbit, any two involution
entries generate an INTRANSITIVE permutation group.

If two involutions generate transitively, their two-colour edge graph is
a single path or an alternating cycle. The path has full-cycle product
and is excluded by Theorem1. In the cycle case both involutions are
fixed-point-free, their index sum is d, and their product has two cycles
and index d-2. They cannot be additive. Smoothness would force the
complementary pair to be additive with sum d-2, giving total index2d-2
and genus0. That is again impossible in positive genus.

Also, FOUR involution entries can never give a positive-genus all-smooth
orbit, in any degree. Their index sum is at most2d. Positive genus would
therefore require g=1 and equality, so all four are fixed-point-free.
Every pair sum is then d>d-1, and no tuple is even smooth initially.

## 2. A general difference-set formulation

Let p be a full d-cycle, labeled as x->x+1. Define

    D_p(a)={y-x !=0 : x,y occur in the SAME a-cycle}.

If (p^k a p^(-k),c) is additive for every k, then

    D_p(a) intersect D_p(c) is empty.                         (2.1)

For a common difference y-x=v-u, choose k=u-x. The shifted a-cycle and
the c-cycle both contain u,v, violating the forest condition.

Equivalently each c-cycle is an independent set in the Cayley graph on
Z/dZ whose nonzero steps are D_p(a). This packages ALL cyclic twists into
one invariant graph, rather than testing just their averaged densities.
It is necessary, not sufficient for additivity or all-smooth closure.

For example, if D_p(a) contains every nonzero element of a subgroup H,
each c-cycle meets each H-coset at most once and has length at most[d:H].
If D_p(a) contains all nonzero differences, c would have to be identity.
Theorem1 uses the union of the two reflection difference sets instead.

## 3. A seven-edge plane-tree lemma

**Lemma.** Let a,c be permutations on14 letters with

    c of type2^7,       i(a)=6,       ca a full14-cycle.

If a,c are additive, label p=ca by x->x+1. Then every nonzero EVEN
difference modulo14 occurs between two letters of some a-cycle.

### The tree and its contour

The additive incidence graph is a connected tree. Each c-vertex has degree
two. Contracting these seven vertices gives a plane tree T with seven
edges and eight vertices, whose vertices are the a-cycles. The cyclic
order at each vertex is its a-cycle; p=c a is the boundary contour.

Each step of this contour crosses a tree edge. The bipartition of T is
therefore precisely the even/odd parity of p-labels. In particular all
letters in an a-cycle have the same parity.

We prove that contour returns at distances2,4,and6 up to sign always
exist. These are all the nonzero even differences modulo14 up to sign.

*Distance2.* Follow a pendant edge out and back. Its neighbor has another
incident edge, so the two contour letters at that vertex are distinct.

*Distance4.* Take a longest path v0,v1,v2,..., with v0 a leaf. If v1 has
degree two, the branch from v2 through v1 to v0 has two edges; traversing
it gives the return. If v1 has degree at least three, all neighbors of
v1 other than v2 are leaves, or the path was not longest. Two of those
leaf branches are consecutive in the plane order, again giving a
two-edge excursion and return distance4.

*Distance6 up to sign.* Choose a centroid of the eight-vertex tree.
Each component after removing it has at most four vertices. Its branch
edge counts are positive integers summing to7, each at most4. A branch
of size3 or4 gives return distance6 or8, respectively, and8=-6 modulo14.
If there is no such branch, all sizes are1 or2. Their odd sum forces a1.
If a2 occurs, some1 and2 are consecutive in the cyclic order and give
an excursion of size3. Otherwise three consecutive1's do so. This proves
the lemma for every plane embedding, without enumerating trees.

The return-distance assertion is special to these small trees. For a
tripod with three arms of length3 (nine edges), the half-differences4,5
are absent. It would be false to extend this lemma to every odd size.

## 4. The ENTIRE genus-three degree-fourteen triangle is excluded

**Theorem 2.** No transitive product-one quadruple of source genus three
and degree14 with indices(6,6,7,13) has an all-smooth full braid orbit.

This is the only possible degree14 TRIANGLE index pattern by the earlier
scalar reduction. It is NOT the only degree14 pattern overall: star
patterns remain open.

### Proof

Call the index7 entry c, the index6 entries a,b, and the index13 entry e.
The last is a full14-cycle. First c has no fixed point. If it had one,
order the tuple(c,a,b,e). Its first product is forced to be full because
7+6=13 while6+13>13. Twisting that pair can align a c fixed point with
a b fixed point (b fixes at least two letters). The resulting(c_k,b)
cannot have full product, contradicting smoothness. This is C209's
fixed-point lemma. A fixed-point-free permutation of index7 on14 points
has seven cycles and must be2^7.

Now normalize the full product p=ca. For every k the pair(a_k,b), with
a_k=p^k a p^(-k), must be additive: its sum is12 and its complementary
sum is20>13. By Section2, b-cycles cannot contain any difference in
D_p(a). Section3 supplies EVERY nonzero even difference. A b-cycle with
at least three letters has two of the same parity, a contradiction.

Thus b has only1- and2-cycles. Its index6 makes its type2^6 1^2, an
involution. Smoothness also forces cb to be a full14-cycle, since7+6=13
and the complementary sum is19. Theorem1 applied to c,b is now a
contradiction. QED.

The proof allows arbitrary a,b cycle types at the start. No involution
restriction, primitivity or full S14 generation was imposed on them.
With C209's exclusion of the degree13 triangle, a genus3 TRIANGLE target
now needs degree at least15. Minimum-index star cases in degrees12--14
have NOT thereby been excluded.

## 5. Three involution entries: a uniform degree bound

**Theorem 3.** If a positive-genus all-smooth orbit has three involution
entries, then

    d>=4g+4;       if d is odd, d>=4g+5.                       (5.1)

Write their indices r1,r2,r3 and their fixed-point counts f1,f2,f3,
so ri=(d-fi)/2; let F=f1+f2+f3. For the fourth index s,

    s=(d+F+4g-4)/2.                                         (5.2)

For each pair of involutions, its complement has sum

    s+ri=2d+2g-2-rj-rk >=d+2g-2 >d-1.

Thus every involution pair must be additive. It cannot have a full
product by Theorem1. Therefore ri+rj<=d-2, or fi+fj>=4, for every pair.
Summing gives F>=6. If d is odd all fi are odd, so F>=7.
The bound s<=d-1 in(5.2) now gives(5.1).

At equality d=4g+4, F=6 and all fi=2: each fi<=F-4=2 and the sum is6.
The fourth entry is full. Hence the boundary passport is uniquely

    (2^(2g+1) 1^2)^3 (4g+4).                              (5.3)

This does not say that(5.3) actually exists as an all-smooth orbit.
For g=3 it is the degree16 passport(2^7 1^2)^3(16). In particular all
three-involution genus3 candidates of degree at most15 are excluded by
the proof, not just the candidates that happened to be sampled.

### A transversal consequence used at the boundary

For an additive involution pair a,b, its union graph consists of paths,
and each path is one cycle of p=ab. Suppose all twisted cross pairs
(a_k,c),(b_k,c) are additive. On EACH path, these rotations of a,b run
through reflections covering every pair of its letters. Thus every
c-cycle meets each path at most once. The twist k may depend on the
chosen path and pair; no simultaneous phase alignment is assumed.

In(5.3), put n=d/2=2g+2. Each pair has two paths. The third involution
has n-1 transpositions, all joining the two paths, and two fixed points.
Their sizes must therefore be

    (n,n) or(n-1,n+1).                                    (5.4)

## 6. The first remaining three-involution genus-three target

For the degree16 passport(2^7 1^2)^3(16), every pair product in a
hypothetical all-smooth component has cycle type8^2 or7 9 by(5.4).
Let theta be the fraction of full-orbit states with first product type7 9.
The three involutions have reciprocal-cycle sum11/2 each; the full entry
has1/16. The C202 normalized Hodge identity consequently gives

    12L=16-265/16+3*((1-theta)/4+theta*16/63),
    L=1/64+theta/1008,
    1/64<=L<=67/4032.                                     (6.1)

Thus ANY all-smooth component here would be non-isotrivial. This part
uses C202's geometric Hodge/compactification argument, with its original
human-review qualification. It is a conditional consequence, NOT a
construction, classifier update or discovery.

### If primitive, the group must be S16

A full16-cycle makes the generated group transitive. Jones's full-cycle
classification (Theorem1.2, case k=0) says that a primitive group not
containing A_d has prime degree, projective degree(q^h-1)/(q-1), or
degree11 or23. Degree16 is none of these: for h=2 it needs q=15, not a
prime power; for odd h>=3 the sum is odd, and even h>=4 with odd q is
at least1+3+9+27=40 (even q gives an odd sum). Hence a primitive group
here contains A16. The full16-cycle is odd, so the group is S16.

This last group-theoretic consequence uses the classification of finite
simple groups through Jones. It is separate from the elementary exclusions.

### What an imprimitive target would have to look like

Proper block sizes divide16 and are2,4,8. Two blocks of size8 are
impossible: every involution has fixed points and hence fixes both blocks,
while the full16-cycle swaps them, violating quotient product one.
Block size2 gives a persistent degree-two deck involution on the source.
The previous C204 genus3 extra-automorphism obstruction rules that out
for a non-isotrivial complete family; this is inherited geometric input.

Thus an imprimitive target must admit four blocks of size4. Its three
induced involutions each have a fixed block, so their quotient indices
are at most1. The fourth entry induces a full4-cycle of index3. If h is
the quotient genus, Riemann--Hurwitz gives

    2h = sum(three quotient indices)-3.

Nonnegativity forces all three indices to be1 and h=0. The quotient
passport is(2 1^2)^3(4). Each source involution has four transpositions
over the swapped block pair; on its two fixed blocks it has respectively
one and two internal transpositions. This is a rational quartic quotient,
not the localized elliptic mechanism of the existing successful examples.

The remaining degree16 target is therefore either S16, or imprimitive
with this rational four-block description (subject to the stated C204
geometric input). Neither case is proved to exist or be impossible here.

## Scope and provenance

The new elementary arguments use only the forest criterion, Hurwitz
twists, Riemann--Hurwitz and a short plane-tree contour proof. Prior
C209 supplies the triangle scalar reduction and fixed-point lemma. C202
supplies the separately qualified Hodge formula; C204 supplies the
generic-automorphism obstruction. Their sealed files were not changed.

Primary group-theory source: Gareth A. Jones, *Primitive permutation groups
containing a cycle*, Bull. Aust. Math. Soc.89(2014),159--165,
https://doi.org/10.1017/S000497271300049X,
https://arxiv.org/pdf/1209.5169, Theorem1.2(k=0). Metadata and the primary
theorem were inspected October22,2026; independent human review is pending.

The scientific-writing and paper-lookup workflows helped separate proved
arguments, inherited geometry and unsuccessful global routes. Workflow
reference, not mathematical evidence: Timothy Kassis, Vinayak Agarwal,
Yuhuan He, Darshil Patel, Aubrey M. Brueckner, *Scientific Agent Skills:
A Library of Procedural Knowledge for Research Agents*(2026),
https://arxiv.org/abs/2609.00065 (API versionv2 inspected).
