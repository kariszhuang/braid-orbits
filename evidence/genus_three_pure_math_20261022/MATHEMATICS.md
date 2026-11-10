# Four first-pair forms, cyclic packing, and the remaining geometric obstruction

October 22, 2026. C210. Mathematical proof draft; human proof review is not
recorded. This is pure mathematical work, not a new enumeration, a production
filter, a deployment, or a change to the running six-worker search. The earlier
C209 finite certificates and their bound files are left untouched.

Throughout, permutations act rightmost first, fixed points count as cycles,
and states are modulo simultaneous S_d conjugation. We consider identity-free,
transitive product-one quadruples and the FULL braid orbit. Write

    ind(u) = d-cycles(u),
    g = 1-d + (sum_i ind(g_i))/2.

Smooth means ind(g1 g2) equals either ind(g1)+ind(g2) or
ind(g3)+ind(g4). All-smooth means smooth at every braid translate.

The strongest new result here is NOT global nonexistence. It reduces the last
mixed part of one degree-fifteen index budget to four explicit first-pair forms
and twenty possible passport triples, with generated group necessarily A15.
No assertion is made that any of those twenty triples has an all-smooth orbit.

## 1. Twisting when the product is not a full cycle

Let (a,b,c,e) be an all-smooth tuple, p=ab, and a_k=p^k a p^(-k).
The braid Q1^(2k) conjugates a,b by p^k and fixes c,e; Q2^(-1) then
places (a_k,c) first. Suppose

    r=ind(a)+ind(c)>0,       ind(b)+ind(e)>d-1.

Every product a_k c must therefore have index r. Set A=Fix(a), C=Fix(c).
If a_k,c shared at least d-r fixed points, their product would act on at
most r other points, hence have index at most r-1 (or index zero if none
remain). Consequently

    |p^k A intersect C| <= d-r-1                 for every k.       (1.1)

Let P range over the cycles of p, of lengths ell_P. Averaging over the
order of p, each point of P is visited equally often. Thus

    sum_P |A intersect P| |C intersect P| / ell_P <= d-r-1.          (1.2)

This is an exact orbitwise density bound, not an independence approximation.
For a full p it gives |A||C| <= d(d-r-1), including C209's fixed-point
obstruction at r=d-1 and its degree-fourteen bound at r=d-2.
For r=d-1, even without a full p, no p-cycle can contain both an a-fixed
point and a c-fixed point.

There is also a cycle-incidence version. Let T be the support of any one
nontrivial cycle of a, and D the support of any cycle of c. Additivity
requires the incidence graph to be a forest, so

    |p^k T intersect D| <= 1,
    sum_P |T intersect P| |D intersect P| / ell_P <= 1.             (1.3)

Two hits would be parallel edges between the same two cycle vertices.
In a full-product additive pair every partner cycle must be hit by the
total active support S, giving in addition average(S,D)>=1. When BOTH
the twisting product p and the tested additive products a_k c are full
d-cycles, (1.3) and that lower bound give, for the longest active a-cycle
length Lmax,

    ceil(d/|S|) <= |D| <= floor(d/Lmax).                            (1.4)

This upper bound can sharpen C209's floor(q*d/|S|); both remain necessary,
not sufficient. For non-full p one must use the separate densities in
(1.2)--(1.3), not replace them by the full-cycle averages.

## 2. A two-product-cycle obstruction in the degree-fourteen triangle

Consider the involutory subcase

    c: 2^7,       a,b: 2^6 1^2,       e: 14.

Smoothness forces ca and cb to be full cycles and ab to be additive.
The union of the a- and b-matching edges is therefore a forest with two
path components; p=ab has exactly two cycles, one on each path. Let their
lengths be u,v, with u+v=14. A path of odd size has one a-fixed endpoint
and one b-fixed endpoint; a path of even size has both endpoints fixed
by the same one of a,b. An isolated vertex counts as an odd path.

Suppose u,v are odd and coprime. Let x,y be the a-fixed points on the two
p-cycles. Since c must connect the two path components for transitivity,
it has a transposition (z,w) with one endpoint in each component.
The Chinese remainder theorem supplies k with

    p^k x=z,       p^k y=w.

Then a_k fixes both z,w. The product a_k c restricts to (z w), so it cannot
be a full fourteen-cycle. But the braid from Section 1 puts (a_k,c) first;
its pair sum is 13 and its complementary sum is 19. This is nonsmooth.

**Conclusion.** In this subcase an odd-odd first-pair partition must have
gcd(u,v)>1. Since u+v=14, its ONLY possibility is 7+7. The odd partitions
1+13, 3+11, and 5+9 are excluded by proof. Even-even path partitions are
not excluded by this argument. Neither is 7+7. This does not complete the
degree-fourteen triangle census.

More generally, the argument works for two odd path components of any
even degree when the outer involution must have a full product with each
twisted low-index involution: coprime path lengths are impossible.

## 3. The remaining degree-fifteen special type

Now take source genus three and indices (4,10,10,10), with special entry

    a: 2^4 1^7.

The other entries b,c,e have five cycles each. C209 proves that each of
their cycle lengths is between 2 and 7. There are exactly seven such types:

    7 2^4; 6 3 2^3; 5 4 2^3; 5 3^2 2^2;
    4^2 3 2^2; 4 3^3 2; 3^5.                         (3.1)

Two ordinary indices sum to 20>14. Thus, at every state, the special entry
paired with any ordinary entry has additive index 14 and a full product.
The three ordinary entries are fixed-point-free. Alternatively, the last
fact follows at once from the full-cycle fixed-point lemma: a fixes seven
points, so any ordinary fixed point could be aligned with one of them.

The pure ordinary triple (3^5)^3 has already been finitely excluded in C209.
Every other ordinary triple in (3.1) has at least one ordinary 2-cycle.
Choose such an entry as c, and a different ordinary entry as b. Normalize

    p=ab: x -> x+1 on Z/15Z,       b=a p,
    F=Fix(a), |F|=7,       S=supp(a), |S|=8.

All arithmetic labels below are 0 through 14 for this cyclic normal form;
add one to obtain the project's usual 1-through-15 labels.

### 3.1 A 2-cycle forces a maximum independent set

Take a 2-cycle {x,x+delta} of c. For every k the pair (p^k a p^(-k),c)
must have a full product. Its incidence tree must hit this c-cycle, so
both endpoints cannot be fixed by that shifted a. Equivalently,

    delta is not in F-F.                                           (3.2)

The step-delta graph on Z/15Z has h=gcd(delta,15) odd cycles of length
15/h. An independent set has at most (15-h)/2 vertices. Since F has seven,
h must be 1. In a fifteen-cycle a seven-point independent set has six
gaps of length two and one gap of length three. After translating labels,

    F = {0,2delta,4delta,...,12delta},
    F-F = Z/15Z minus {delta,-delta}.                              (3.3)

Consequently EVERY 2-cycle of c, and also of e, has cyclic step +/-delta.
The assertion for e uses twists of the same first pair, leaving e fixed;
one can first swap c,e and then move e next to a. It is not an assertion
that b stays fixed under those twists.

Up to replacing delta by -delta and translating F, the only choices are
delta=1,2,4,7. The sign replacement is legitimate because the seven-term
progressions for delta and -delta differ by a translation.

### 3.2 The step-one case is impossible

If delta=1, c has a transposition (x,x+1). Product one gives

    e=c^(-1) p^(-1),       e(x+1)=c^(-1)(x)=x+1.

This contradicts fixed-point-freeness of e. Thus delta!=+/-1.

### 3.3 The step-seven case is impossible

For delta=7, F={0,9,10,11,12,13,14}, seven consecutive points on the
p-circle. Since b=a p and a fixes this interval, b follows seven successive
p-edges. It therefore has a cycle of length at least eight. That contradicts
the upper bound seven in (3.1).

One can avoid that upper bound: all eight support points are consecutive,
and any noncrossing perfect matching of them has an adjacent p-edge, giving
a fixed point of b. The noncrossing fact is proved in the next subsection.

### 3.4 Noncrossing chords leave three plus one possibilities

The equality ind(a)+ind(b)=ind(p) means that applying the four disjoint
transpositions of a to p splits a cycle at every step. Their chords on
the p-circle cannot interlace: interlacing chords would cause a join after
the preceding split. Conversely noninterlacing disjoint chords perform
four splits. Hence a is a noncrossing matching of S. An adjacent p-edge
(x,x+1) in a would make b(x)=a(x+1)=x, and is forbidden.

For delta=2,

    F={0,1,4,5,8,9,12},       S={2,3,6,7,10,11,13,14}.

The support consists of four adjacent pairs. On eight cyclically ordered
support vertices v0,...,v7, exclude the four edges v0v1,v2v3,v4v5,v6v7.
Recursing on the partner of v0 leaves exactly these three matchings:

    (v0 v3)(v1 v2)(v4 v7)(v5 v6),
    (v0 v7)(v1 v2)(v3 v4)(v5 v6),
    (v0 v7)(v1 v6)(v2 v5)(v3 v4).

Indeed v0 may meet v1,v3,v5,v7. The first is forbidden; v3 gives the
first row; v5 leaves the forbidden edge v6v7; v7 gives the last two rows.
This is a four-chord proof by cases, not a scan of permutations.

For delta=4,

    F={0,1,2,3,8,9,10},       S={4,5,6,7,11,12,13,14}.

A chord internal to either four-point interval forces an adjacent p-edge
inside that interval, by noncrossing and recursion. Thus every chord goes
between the two intervals. Noncrossing then forces the unique nested
matching (4 14)(5 13)(6 12)(7 11).

We obtain exactly FOUR NECESSARY FIRST-PAIR FORMS, with p fixed as +1:

| delta | a | b=a p cycle type |
|---|---|---|
| 2 | (2 7)(3 6)(10 14)(11 13) | 6 3 2^3 |
| 2 | (2 14)(3 6)(7 10)(11 13) | 4 3^3 2 |
| 2 | (2 14)(3 13)(6 11)(7 10) | 5 3^2 2^2 |
| 4 | (4 14)(5 13)(6 12)(7 11) | 5 4 2^3 |

These are not four braid states or four completed quadruples. Completing
c,e and full braid closure are still necessary. Their pair centralizers
are trivial: a rotation commuting with a would preserve its eight-point
support, impossible for a nonidentity rotation of a fifteen-cycle, whose
point-orbit lengths are 3,5,or15. The four forms are distinct because their
b cycle types are distinct.

### 3.5 Every remaining ordinary type is one of four

The preceding proof allowed arbitrary b in (3.1), not just b of type 3^5.
In particular if b=3^5 and c has a 2-cycle there is NO possible first pair.
So a mixed triple cannot contain 3^5. Once that type is removed, all three
ordinary entries have 2-cycles. Choose any of them as b and a different
one as c, and repeat the argument. Every ordinary entry must therefore have
one of the FOUR types in the table.

There were binomial(9,3)=84 unordered ordinary triples with repetition
among the seven types. One pure triple was already excluded in C209.
The new proof excludes 83-binomial(6,3)=63 MIXED triples. Exactly twenty
unordered triples of the four remaining types are not excluded by this
argument. These are counts of passport triples, not tuples or components.
The conclusion is invariant under reordering the branch classes by braids.

### 3.6 The first pair already generates A15

This deduction needs no classification of primitive groups. In (3.3), F
and F+delta are disjoint, so S and S+delta meet in exactly one point.
The union of the four a-edges and four p^delta a p^(-delta)-edges has one
three-vertex path and six isolated edges. Their product therefore has type

    3 2^6.

Squaring gives a SINGLE 3-cycle in H=<a,p>=<a,b>.

H is primitive. Any block system invariant under its full fifteen-cycle
p is the congruence partition modulo 3 or modulo 5. The seven-term
progression F meets all of these residue classes. Hence the induced
involution a fixes every block. But its fixed-point counts in the blocks
are, respectively, a permutation of (3,2,2) or of (2,2,1,1,1).
An involution on an odd-sized block has an ODD number of fixed points,
contradicting either distribution. Both possible proper block systems fail.

Jordan's three-cycle theorem now gives A15<=H. Both a and p are even, so

    <a,b>=A15,       and hence <a,b,c,e>=A15.                       (3.4)

For a primary implementation of the exact standard theorem used, see
Mathlib's Jordan module, theorem
Equiv.Perm.alternatingGroup_le_of_isPreprimitive_of_isThreeCycle_mem.
The proof above has NOT itself been formalized in Lean.

Thus imprimitive/isotrivial constructions cannot supply completions in
this remaining mixed scope. This does not assert that primitivity alone
implies non-isotriviality in other scopes.

### 3.7 A quantitative non-isotriviality consequence

This subsection additionally uses the C202 arbitrary-profile Hodge formula,
with its geometric compactification hypotheses and pending human review;
the elementary reductions above do not depend on it.
Let L=deg(lambda)/D for an actual finite family base of degree D over
Mbar_0,4. Put rho(u)=sum over cycles of 1/length. All pair products in this
index budget have the type of a full fifteen-cycle, so Rbar=1/15.
Here rho(a)=9 and

    rho(6 3 2^3)=2,       rho(5 4 2^3)=39/20,
    rho(5 3^2 2^2)=28/15, rho(4 3^3 2)=7/4.

The formula becomes

    12 L = 31/5 - rho(b)-rho(c)-rho(e),
    1/60 <= L <= 19/240.                                          (3.5)

Therefore ANY all-smooth component among the twenty remaining triples
would be non-isotrivial. No such component is constructed here. Its full
all-smoothness cannot be replaced by a short walk, and the rational number
L is not itself the integer Hodge degree on an arbitrary coarse base.

## 4. A necessary geometric event: hyperelliptic fibres

This is a standard geometric obstruction, not a new literature theorem.
All-smoothness gives, after the usual finite base change and stabilization,
a smooth complete family f:X->B of source curves. For a non-isotrivial
genus-three family, its Hodge degree Lambda=deg(det f_*omega)>0.
Such a family cannot be entirely hyperelliptic: the hyperelliptic moduli
space is affine (finite branch-configuration quotient), so its complete
moduli image would be a point. Thus its general fibre is nonhyperelliptic.

Let E=f_*omega (rank3), V=f_*omega^2 (rank6). Multiplication of canonical
forms gives nu:Sym^2 E->V. It is an isomorphism away from hyperelliptic
fibres. Its determinant vanishes precisely on the hyperelliptic divisor.
For a smooth family, Grothendieck--Riemann--Roch gives

    c1(V)=13 lambda,       c1(Sym^2 E)=4 lambda.

Therefore

    sum_{b in B} ord_b det(nu) = 9 Lambda > 0.                     (4.1)

The equality counts intersection multiplicities, not distinct fibres.
Every hypothetical target MUST acquire at least one SMOOTH hyperelliptic
fibre somewhere on its compactified family base. A proof that a candidate
family cannot meet that divisor would prove its nonexistence.

The hyperelliptic fibre may occur over a Hurwitz boundary point, after
rational tails contract. A nonhyperelliptic generic representative, or
checking covers at a few cross-ratios, does not rule it out. Nor does its
hyperelliptic involution have to commute with the covering map: it need
not be visible in the permutation centralizer. This is different from
the earlier obstruction to an automorphism persisting on EVERY fibre.

For the degree-fifteen remainder, (3.5) would give the normalized count

    3/20 <= (hyperelliptic intersection degree)/D <= 57/80.

This is a bridge between the permutation data and a concrete missing
geometric proof, not an automatic tuple classifier. The standard facts
and the hyperelliptic divisor formula are in de Jong, Section9,
Proposition9.1 and Remark9.7 [E2].

## 5. Why some attractive global arguments do not finish the proof

1. Complete smooth genus-three families can vary; a genus-two argument
   cannot simply be reused in M3. Flapan's Theorem9.1, case(2a), even
   constructs complete families with a fixed two-dimensional rational
   Hodge summand and varying four-dimensional complement [E1]. Thus a
   fixed elliptic factor alone is not a nonexistence proof.
2. A simple moving branch point gives a Schiffer-type infinitesimal
   deformation, but its rank-two Higgs kernel in genus three need not
   be a fixed rank-two local subsystem. Gonzalez-Alonso--Torelli explicitly
   distinguish the Higgs kernel from the unitary flat part [E3]. Their
   semistable examples are not asserted to be everywhere smooth. A new
   argument using completeness AND no singular source fibres is needed.
3. Nonexistence of complete curves in the NONPROJECTIVIZED stratum of
   holomorphic differentials cannot be applied to a family supplying only
   a projective differential. One must also control scaling, zeros of
   that scaling, and changes of zero type at the Hurwitz cusps [E4].
4. A finite upper degree cutoff has not been proved. Degree-twelve's
   two-product-cycle star budget (1,9,9,9), the other star budgets, and
   higher degrees remain open. No failure rate from a biased search is
   converted here into a probability of universal nonexistence.

## 6. The precise next mathematical questions

**Finite obstruction, first priority.** Can product one and the four frames
of Section3 force a later nonforest incidence graph in all twenty remaining
ordinary passport triples? In the first three frames all external 2-cycles
have p-step +/-2; in the last they have step +/-4. Use this simultaneous
constraint for c and c^(-1)p^(-1), not another random forest start. A proof
of escape would complete this mixed special-type slice, not the entire degree.

**Geometric obstruction.** Can one show that a genus-three all-smooth
Hurwitz family cannot meet the smooth hyperelliptic divisor? Section4 shows
exactly why that would settle GLOBAL nonexistence. At present this is an
unproved question, and possible hyperelliptic stabilized cusp fibres must
be included. It is stronger than saying the generic curve is a plane quartic.

**Unsaturated twists.** Use (1.2)--(1.3) on the separate cycles of the pair
product. The full-cycle proof cannot decide the degree-twelve collision
budget merely by pretending its two product cycles form one rotation.

**No symmetry assumption.** The four-frame argument forced A15 without
assuming blocks, a special quotient, or a primitive monodromy target. Its
combination of repeated twists, extremal independent sets, noncrossing
chords and a three-cycle witness is the useful general pattern to extend.

## References and provenance

Literature was opened at the primary full-text sources below on October22,
2026; metadata was resolved with the arXiv Atom API and its strict parser.
These are source-inspected draft dependencies; no accountable human source
verification or proof approval is recorded. No literature-priority claim is
made for the new elementary derivations. See SOURCES.json for exact queries.

- E1: Laure Flapan, *Monodromy of Kodaira Fibrations of Genus3*,
  arXiv:1709.03164v2, Theorem9.1, especially realized case(2a), and Section8.
  https://arxiv.org/pdf/1709.03164v2
- E2: Robin de Jong, *On the height of Gross--Schoen cycles in genus three*,
  Research in Number Theory4,38(2018), DOI10.1007/s40993-018-0131-0,
  Section9, Proposition9.1 and Remark9.7. https://arxiv.org/pdf/1711.02990v2
- E3: Victor Gonzalez-Alonso and Sara Torelli, *Families of curves with
  Higgs field of arbitrarily large kernel*, DOI10.1112/blms.12437;
  arXiv:1812.05891v2, Introduction and Lemmas2.2--2.3.
  https://arxiv.org/pdf/1812.05891v2
- E4: Dawei Chen, *Complete curves in the strata of differentials*,
  Comptes Rendus Mathematique361(2023),371--373,
  DOI10.5802/crmath.416, theorem and the following definition of P_g(mu).
  https://arxiv.org/pdf/2208.10717
- E5: Mathlib contributors, *Mathlib.GroupTheory.GroupAction.Jordan*,
  the primitive-subgroup/three-cycle theorem (Wielandt13.3),
  https://leanprover-community.github.io/mathlib4_docs/Mathlib/GroupTheory/GroupAction/Jordan.html

Research-workflow reference (not mathematical evidence): Timothy Kassis,
Vinayak Agarwal, Yuhuan He, Darshil Patel, Aubrey M. Brueckner,
*Scientific Agent Skills: A Library of Procedural Knowledge for Research
Agents* (2026), arXiv:2609.00065. https://arxiv.org/abs/2609.00065

The scientific-writing procedure influenced the separation between proved
reductions, conditional geometric consequences, and open proof steps; the
paper-lookup procedure supplied metadata/provenance, not mathematical proof.
