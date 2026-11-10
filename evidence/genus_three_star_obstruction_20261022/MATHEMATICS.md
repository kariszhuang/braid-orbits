# Closing the degree-fifteen mixed star case

October 22, 2026. C211. Proof draft; accountable human proof review is not
recorded. This is pure mathematical work. The six-worker C209 search, its
configuration, scientific binaries, sealed evidence and public data are
unchanged. No completion enumeration or new braid-orbit search was run.

Permutations act rightmost first. Fixed points count as cycles. All-smooth
means smooth at EVERY state of the FULL braid orbit, modulo simultaneous
S_d conjugation, not just at a chosen representative.

## 1. The new theorem

**Theorem.** A transitive product-one degree-fifteen quadruple of source
genus three, with indices (4,10,10,10), special entry of type 2^4 1^7,
and at least one ordinary entry not of type 3^5, cannot be all-smooth.

This closes ALL of the twenty mixed passport triples left open in C210.
It is an elementary proof, independent of Hodge degrees, isotriviality
criteria, a group-generation assumption, and the earlier finite negative
certificates. The previously completed pure triple (3^5)^3 remains a
separate finite result; it is not counted as part of the new proof.

### Dependencies from C210, with the key argument recalled

Write the tuple (a,b,c,e), with a the special involution and F=Fix(a).
Every special/ordinary pair has additive index 14 and a full product;
ordinary/ordinary index sums are 20>14. Normalize p=ab to x->x+1 on
Z/15Z. Put S=supp(a), |S|=8, |F|=7. All ordinary cycles have lengths
between 2 and 7, so the seven possible ordinary types are

    7 2^4; 6 3 2^3; 5 4 2^3; 5 3^2 2^2;
    4^2 3 2^2; 4 3^3 2; 3^5.

If an external ordinary entry c has a 2-cycle with cyclic step delta,
outer twists of a,b force delta not in F-F. The step-delta graph is a
union of gcd(delta,15) odd cycles. Its independent-set bound forces delta
to be a unit and F to be a maximum seven-point independent set. Hence,
after a cyclic translation,

    F={0,2delta,...,12delta},       F-F=Z/15Z minus {+/-delta}.

Every external 2-cycle of c or e consequently has step +/-delta. The
possibilities delta=1,2,4,7 up to sign reduce to 2 or4: step1 gives a
fixed point in e=c^(-1)p^(-1), and step7 gives an ordinary cycle longer
than seven. Noncrossing matching chords for a then give exactly the four
first-pair forms in C210 Section3.4. We use that proved finite proof by
cases, not an unrecorded search. Labels here are 0,...,14.

A mixed all-smooth tuple cannot contain an ordinary 3^5 entry: choose
that entry as b and a different ordinary 2-cycle entry as c; none of the
four first-pair forms has b of type3^5. Thus ALL three ordinary entries
have 2-cycles, and branch reordering puts every special/ordinary pair
among the four forms. This explicitly supplies the hypothesis needed
when comparing an external c with its own first-pair form below.

## 2. An incidence-degree invariant removes nineteen of twenty triples

For an ordinary cycle D define its incidence degree to be |D intersect S|.
This is a simultaneous-conjugacy invariant of the PAIR (a,ordinary entry).
In the cycle-incidence tree for an additive pair, it is precisely the
degree of the ordinary-cycle vertex, with the four active a-cycles
contracted into edges. The four C210 first-pair forms have these data:

| Ordinary type | Number of 2-cycles with incidence degree TWO | Incidence tree |
|---|---:|---|
| 6 3 2^3 | 2 | path |
| 5 4 2^3 | 3 | path |
| 5 3^2 2^2 | 2 | path |
| 4 3^3 2 | 0 | star |

For an EXTERNAL ordinary entry c in the p-frame, all its 2-cycles have
step +/-delta. Because F is a maximum independent set, S contains exactly
ONE edge of the step-delta fifteen-cycle. Thus c has AT MOST ONE 2-cycle
with both endpoints in S, i.e. of incidence degree two.

But (a,c) must itself be one of the same four first-pair forms, since the
other ordinary entries still contain 2-cycles in a mixed tuple. Its
intrinsic incidence degrees must match the table. The first three rows
are impossible: they require two or three internal 2-cycles, not at most
one. Hence c has type

    U = 4 3^3 2.

The same reasoning applies to e, and then to b by choosing a different
ordinary entry as the first partner. All ordinary entries must have type
U. This leaves just U^3 of the twenty necessary triples from C210.

No orbit counts enter this argument. In particular, it compares an
intrinsic invariant of (a,c) with constraints imposed from a DIFFERENT
first pair; it never incorrectly holds b fixed during a twist of a,b.

### The same bivalent-vertex bound in other genera

More generally take q>=2, degree d=4q-1, special type2^q 1^(2q-1), and
three ordinary entries with q+1 cycles each. Their indices are3q-2,
and Riemann--Hurwitz gives source genus q-1. Suppose at least two of
the ordinary entries have 2-cycles. For ANY chosen first partner j,
another ordinary entry has a 2-cycle. The same twist argument makes
Fix(a), of size(d-1)/2, a maximum independent set of a unit-step odd
d-cycle. Its complement contains exactly one edge of that step cycle.

For an ordinary entry with a 2-cycle, choose a DIFFERENT first partner so
that this entry is external; the bound just proved applies to its own
2-cycles. Entries with no 2-cycles satisfy the conclusion vacuously.
Therefore every special/ordinary incidence tree has AT MOST ONE ordinary
2-cycle vertex of degree two. This is an all-degree, all-genus necessary
condition in the stated extremal-involution family. The degree-fifteen
table above exploits its especially restrictive q=4 case.

## 3. Star centres and a general intersection lemma

Normalize the remaining first pair to

    a=(2 14)(3 6)(7 10)(11 13),
    b=(0 1 14)(2 6 10 13)(3 4 5)(7 8 9)(11 12),
    p=ab: x->x+1,
    F={0,1,4,5,8,9,12},
    S={2,3,6,7,10,11,13,14}.

Every pair (a,j), j=b,c,e, has the star incidence tree in the last table
row. Let T_j be the support of its 4-cycle, called its CENTRE. It contains
one endpoint of each of the four a-transpositions. The other four cycles
of j each contain exactly ONE point of S: three cycles contain two fixed
points and one contains one fixed point. Thus both j and j^(-1) map

    T_j to T_j,       S minus T_j into F.

Put B=T_b={2,6,10,13}, R=S minus B, C=T_c, E=T_e. Since p=a b,

    p(B)=R,       p(R) is contained in F.

Using p=e^(-1)c^(-1), count points x in C for which p(x) remains in S:

    |B intersect C| = |C intersect E|.

Using p^(-1)=c e, count points x in E for which p^(-1)(x) remains in S:

    |B intersect E| + |C intersect E| = 4.

Write k=|B intersect C|=|C intersect E|. All three centres are transversals
of the four a-pairs. If D,H are the sets of pairs where C,E choose the
opposite endpoint to B, then

    |D|=4-k,       |H|=k,       |D symmetric_difference H|=4-k.

It follows that |D intersect H|=k/2. Thus k is even. The choice k=4 is
impossible because then D is empty, so

    k=0 or2.                                                     (3.1)

**General centre lemma.** The same proof works with q a-transpositions
and star centres of size q, provided all other ordinary cycles each
have length at least two, contain one active point, and product a b c e=1.
It gives

    |B intersect C|=|C intersect E|=k,
    |B intersect E|=q-k,
    k even,       0<=k<=2q/3.                                    (3.2)

This lemma does not depend on degree fifteen or on a genus assumption.
It is a restriction on this stated star-incidence structure, not on
arbitrary all-smooth quadruples.

## 4. Disjoint centres are impossible

If k=0, then C=R and E=B. Product one gives

    c(B)=p^(-1)(B)={1,5,9,12},
    e^(-1)(R)=p(R)={0,4,8,12}.

The unique 2-cycle of c consists of one point of B and one of F with
step +/-2. Its fixed-point endpoint must lie in c(B), so the ONLY
possibility is (10 12). Similarly, the unique 2-cycle of e has one point
of R and one of F; its fixed endpoint lies in e^(-1)(R), leaving only
(14 12). But then

    c e(14)=c(12)=10 !=13=p^(-1)(14).

Contradiction. Consequently k=2.

## 5. Two-point intersections give four, not twenty-four, configurations

The membership triples of points of S in B,C,E now identify S with the
eight vertices of a cube {0,1}^3. Take the centres to be the coordinate
faces x=0,y=0,z=0, respectively. Each a-pair flips all three coordinates.
The transversal intersection counts imply that every membership triple
occurs exactly once, so these are actual distinct labels.

For x in C, p(x) stays in S exactly when x lies in B. Since
p=e^(-1)c^(-1), this forces

    c^(-1)(C intersect B)=C intersect E.

For x in E, p^(-1)(x) stays in S exactly when x lies in R. Thus

    e(E intersect R)=E intersect C.

Each of these equations maps a two-element subset of a four-element
centre onto another two-element subset meeting it once. There is exactly
one FOUR-CYCLE doing this: if the subsets are {u,v},{u,w}, the cycle must
be (v u w z). A fixed point at u would not be a four-cycle.

In cube labels the two central cycles therefore are

    c_4=(000 001 101 100),
    e_4=(110 100 000 010).

For x=000,001, use p=e^(-1)c^(-1) and b=a p. This gives

    b(000)=001,       b(001)=011.

Since b has a four-cycle on B, its order there must be

    b_4=(000 001 011 010).

Hence the entire cube labeling is determined by t=000 in the four-point
set B. In the original labels, if b_4=(2 6 10 13), then

    C={t,b t,a b^2 t,a b^3 t},
    E={t,b^3 t,a b^2 t,a b t},
    c_4=(t,b t,a b^3 t,a b^2 t),
    e_4=(a b t,a b^2 t,t,b^3 t).                                 (5.1)

There are only FOUR cases, t=2,6,10,13. No choice of a completion's
long-cycle orientations remains in the centres.

## 6. Four explicit contradictions

For y in E minus C, product one forces

    c(y)=p^(-1)e_4^(-1)(y).                                     (6.1)

If (s f) is the unique 2-cycle of e, then s lies in S minus E, f lies
in F, f-s=+/-2 modulo15, and

    c(s)=f-1,       c(f)=s-1.                                    (6.2)

To list possibilities in the following table, use the known c_4 when
s lies in C, and c(s) in F otherwise. These exhaust the candidates
because a noncentral cycle has exactly one active point. Entries c(y)
in the fourth column come directly from (6.1).

| t | c_4 | e_4 | Forced external c-values | Possible e 2-cycle after these necessary constraints |
|---:|---|---|---|---|
| 2 | (2 6 11 7) | (3 7 2 13) | c(3)=12, c(13)=1 | (14 1) |
| 6 | (6 10 14 11) | (7 11 6 2) | c(2)=5, c(7)=1 | (3 1), (3 5), (14 12) |
| 10 | (10 13 3 14) | (11 14 10 6) | c(6)=9, c(11)=5 | (7 5), (7 9) |
| 13 | (13 2 7 3) | (14 3 13 10) | c(10)=12, c(14)=9 | (11 9) |

For completeness the fixed-point neighbours at cyclic distance two are

    2:{0,4}, 3:{1,5}, 6:{4,8}, 7:{5,9},
    10:{8,12}, 11:{9}, 13:{0}, 14:{1,12}.

This small neighbour list and (6.2) let one check the last column by hand;
it is not an unexplained computation or an enumeration of completions.

### t=2

The only e 2-cycle is (14 1). Equation(6.2) gives c(1)=13, but (6.1)
already gives c(13)=1. Thus c has the 2-cycle (1 13), whose cyclic step
is three, not two. Contradiction.

### t=6

First (3 5) as the e 2-cycle would give c(5)=2 while c(2)=5, producing
another forbidden step-three 2-cycle.

The c 2-cycle is forced to be (13 0). Indeed its active endpoint is in
S minus C={2,3,7,13}. The forced values c(2)=5,c(7)=1 exclude endpoints
2 and7. Endpoint3 would need fixed neighbour1 or5; but a c 2-cycle (3 f)
would require e(f+1)=3. The known e_4 has e(2)=7 and e(6)=2, ruling out
both choices. Only endpoint13 with fixed neighbour0 remains.

From c(13)=0,c(0)=13 and e=c^(-1)p^(-1) follow

    e(1)=13,       e(14)=0.

These respectively rule out the remaining possible e 2-cycles (3 1)
and (14 12). Contradiction.

### t=10

If the e 2-cycle is (7 9), (6.2) gives c(9)=6 while c(6)=9, again a
forbidden step-three 2-cycle. If it is (7 5), then c(5)=6 while c(11)=5.
The c-cycle containing 11->5->6 contains TWO points of S minus C,
namely11 and6. Neither its central four-cycle nor any of its leaf cycles
can do this. Contradiction.

### t=13

The only e 2-cycle (11 9) forces c(9)=10. Together with c(14)=9 this
places the two distinct external active points14,10 in the same c-cycle.
The star leaf structure forbids that. Contradiction.

All four cases fail. This proves the theorem.

## 7. Consequences and exact scope

1. The twenty mixed ordinary passport triples left in C210 are now all
   excluded. With C210's preceding63, all83 mixed triples in this special
   degree-fifteen index budget are excluded by elementary arguments.
2. Together with the PREVIOUS finite pure-case certificates in C209,
   every all-smooth genus-three degree-fifteen tuple of indices
   (4,10,10,10) must have special type5 1^10 and ordinary triple(3^5)^3.
   That remaining passport is isotrivial by the earlier uniform cubic
   Hodge identity. This corollary combines the new proof with established
   finite results and geometric input; do not call the entire corollary
   a computation-free classification.
3. A non-isotrivial genus-three all-smooth orbit with MINIMUM INDEX FOUR
   must have degree at least16. In the small-index star regime, minimum
   index s gives d>=2s+2g+1; equality forces all three other indices
   s+2g. For s=4,g=3 equality is exactly the budget just closed. The
   triangle regime has minimum index at least2g=6, so contributes no
   minimum-four exception. This is NOT a universal lower degree sixteen:
   minimum indices1,2,3 and other degree-fourteen triangle cases remain.
4. Why genus3 was unusually rigid here: at the minimum star degree
   d=6g-3 with special involution2^(2g-2), its fixed-point set has2g+1
   elements. It equals the largest possible independent set on the odd
   d-cycle, (d-1)/2=3g-2, EXACTLY when g=3. In larger genus this specific
   packing argument has slack g-3 and needs additional ideas.

No global nonexistence theorem, new all-smooth component, or exclusion of
the complete degree fifteen has been obtained. The search is unchanged.
The global curve-lifting/homology formulation is developed separately in
[GLOBAL_FRAMEWORK.md](GLOBAL_FRAMEWORK.md).

## Proof provenance

The four-form lemma, its exact populations and hypotheses are in
../genus_three_pure_math_20261022/MATHEMATICS.md, Section3 (C210).
The old pure finite negatives are bound in C209; their files were not
modified or rerun. The present centre and leaf arguments are written
derivations, not claims of literature priority or human-reviewed theorems.

The scientific-writing workflow helped separate the new elementary proof,
its inherited dependencies, the conditional corollaries and the remaining
global question. Workflow reference, not mathematical evidence: Timothy
Kassis, Vinayak Agarwal, Yuhuan He, Darshil Patel, Aubrey M. Brueckner,
*Scientific Agent Skills: A Library of Procedural Knowledge for Research
Agents* (2026), arXiv:2609.00065, https://arxiv.org/abs/2609.00065.
