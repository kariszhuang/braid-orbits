# Fixed cyclic quotient and one moving branch value

Work over C. A tuple describes a branched cover, not by itself a varying
curve. Non-isotriviality below concerns its associated connected Hurwitz
family with four distinct branch values, after finite marking/base changes.
States use simultaneous conjugation by S_d.

## The criterion

Let n=rk with r odd. A transitive product-one tuple has an invariant partition
into r blocks of size k. Its block monodromy is cyclic: one entry a induces
identity and acts nontrivially inside just one block; the other three induce
unit translations on Z/r and each has cycle type r^k. Assume a is nonidentity.

The Hurwitz family factors through a fixed intermediate curve E of genus
h=(r-1)/2, and the degree-k map C_t -> E has exactly one branch value P_t.
This point moves. Consequently:

- For r>=5 the source family C_t is **non-isotrivial**.
- For r=3 the source family is **isotrivial**.

These are sufficient conditions, not an equivalence for arbitrary smooth
components, and not a claim that every source curve from different components
or parameter values is distinct.

## 1. Why an invariant block system gives an actual curve factorization

Normalize the three large branch values to 0,1,infinity and write t for the
fourth. This may require ordering the branch values, a finite cover of the
original parameter space. On a small parameter disk keep their monodromy
markings and the block system fixed.

For U=P^1 minus {0,1,infinity,t}, let Gamma=pi_1(U). Let H be the stabilizer
of a point in the transitive monodromy action, and K the stabilizer of its
block. Then H <= K <= Gamma. The corresponding covering spaces have degrees
k and r, giving C_t^o -> E_t^o -> U. Extend them across punctures with the
local power maps determined by permutation cycles. The resulting compact
smooth Riemann surfaces give an actual holomorphic factorization

    C_t --degree k--> E_t --degree r--> P^1.

This is the subgroup description of covering spaces together with Riemann
existence, not a supposition inferred from the sizes of cycles. See Debray's
notes, Lecture11, Theorem11.2 and the immediately following construction.

## 2. Why E is fixed, including an explicit equation

The small entry acts identically on blocks, so E_t -> P^1 is unramified at t.
Its remaining monodromies are translations by u,v,w, each a unit modulo r,
with u+v+w=0 modulo r. They generate the regular cyclic group C_r. Removing
the extra unramified marked fiber over t leaves the same three-point cover
for every t. Thus identify E_t with a fixed E.

More explicitly E is the smooth projective normalization of

    y^r = x^u (x-1)^v,

where 1<=u,v<=r-1 and gcd(r,u)=gcd(r,v)=gcd(r,u+v)=1. A positive loop around
0 or1 multiplies y by the corresponding rth root of unity. This gives exactly
the required cyclic quotient monodromy. Its genus satisfies

    2h-2 = -2r + 3(r-1) = r-3,

so h=(r-1)/2. The equation contains no t. It is an equation for the quotient
E, **not** an equation for the source C_t. Ordering the large branch values
in each state can change the presentation, without altering fixedness.

## 3. Why C_t -> E has exactly one, genuinely moving, branch value

At each large branch value, E -> P^1 has one ramification point of index r.
Every point of C_t above it has index r as well, because the entry has type
r^k. Ramification indices multiply in a composition, so every relative index
of C_t -> E there equals r/r=1: no relative branching is hidden there.

Over t the quotient has r unramified points, one for each block. The local
entry is identity on all but one block. Consequently C_t -> E is unramified
at all but one of those r points, and ramified at the remaining point P_t.
Its ramification indices are exactly the cycle lengths of a inside that
block (including any fixed points). This includes a double transposition:
two ramification points upstairs, but only **one branch value** downstairs.

The fixed map x:E -> P^1 satisfies x(P_t)=t. Since the four-branch Hurwitz
parameter t varies on an open disk, P_t cannot be constant. A marking or a
finite dominant base change can choose its branch, but cannot make t constant.
It is precisely this moving point, not the source's genus or state count,
that drives the non-isotriviality argument.

If J=ind(a), the relative Riemann-Hurwitz formula independently gives

    2g(C_t)-2 = k(2h-2) + J,
    g(C_t) = 1 + (k(r-3)+J)/2.

## 4. Non-isotriviality when h>=2: finiteness proof

Suppose all C_t are isomorphic to a fixed curve C. Finite source
automorphisms allow local trivialization (or finite marking/base change),
but the following contradiction does not even require choosing those
isomorphisms continuously. For every t choose an isomorphism C -> C_t.
Compose it with C_t -> E to obtain a nonconstant degree-k map f_t:C -> E.
Its unique branch value is P_t, because precomposing by an isomorphism does
not change the branch locus in E.

For fixed compact Riemann surfaces C and E with genus(E)>=2, de Franchis
finiteness says there are only finitely many nonconstant holomorphic maps
C -> E. They provide only finitely many branch loci, whereas x(P_t)=t
gives infinitely many distinct branch values. Contradiction. Target
automorphisms cannot remove this variation: E is fixed and its automorphism
group is finite. Source genus is also >=2 by the preceding formula.

Primary reference: Howard and Sommese, *On the theorem of de Franchis*
(1983), introduction p429 and Theorem1 p433. Finiteness for a fixed target
is all that is needed. The paper's PDF text extraction misreads some >= signs;
the criterion uses genus **at least two**, not strictly greater than two.

## 5. Local rigidity cross-check (independent proof of the obstruction)

If an isotrivial family is locally trivialized, the degree-k maps form a
holomorphic family f_t:C -> E with fixed source and target. Differentiating
in a parameter direction gives a holomorphic section of f_t^* T_E.
The degree of this line bundle is k(2-2h), which is negative when h>=2.
A negative-degree line bundle on a compact curve has no nonzero holomorphic
section: a nonzero section would have an effective zero divisor of that
negative degree. Hence all parameter derivatives vanish and f_t is locally
constant. Its branch value cannot move. This yields the same contradiction
and explicitly includes h=2. The audit records the negative degree but does
not pretend to machine-check complex geometry.

## 6. Elliptic control: why the same argument must not be used for r=3

For h=1 the tangent degree is zero and E has translations. Translate P_t
to a fixed origin O on E. The composed covers have the same fixed target
E and branch value O. Their restrictions to E minus {O} are degree-k
unramified covering spaces. This punctured elliptic curve has a finitely
generated fundamental group, with only finitely many homomorphisms to S_k.
Thus only finitely many covering types occur. Each has a unique holomorphic
completion for this fixed E. A connected holomorphic family cannot vary in
this finite set of source isomorphism types; its source is isotrivial.
Translation changed the map, not the abstract source. This explains why
higher genus of C alone does not prove non-isotriviality.

## Evidence boundary

The independent audit checks every saved state in all52 closed components,
including all21,648 states of the ten higher-genus cases. Their group orders
and literal six-move closure come from the prior independently verified C138
certificates, bound here by exact table/catalog hashes. The new checks verify
the block factorization hypotheses, quotient equations, relative branching,
two genus formulas, smoothness/defect, graph connectivity and inverse edges.

The classification is a theorem-based deduction over C, not merely absence
of a numerical counterexample. It is not a formal proof-assistant certificate,
an explicit source-curve/period-matrix computation, recorded human approval,
or a literature novelty claim. Ten distinct braid components need not give
ten disjoint loci after forgetting the covering map. Dated C138 evidence is
preserved; this mission adds the explicit bridge and independently audits it.
