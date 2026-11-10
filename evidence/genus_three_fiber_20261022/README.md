# C196: solve elliptic return equations instead of guessing

Goal remains **one independently closed all-smooth source-genus-three component
with a separate proof of non-isotriviality**. No such example is claimed here.
The public dataset, public classifier, nine weak counterexamples, open TeX and
backup state are unchanged. This is an optimization and a new search population.

## Why this is not merely a faster repetition

The prior elliptic lane allowed moving relative index only 1, 2 or 3 and
required positive ordinary excess. In genus three the total relative excess is
four. It therefore omitted the case where **all four units occur above the
moving fibre, spread over all its points**, with zero ordinary excess.
The familiar isotrivial construction instead has just one active block.

For three blocks distribute moving indices as (2,1,1), and for four blocks use
(1,1,1,1). Use unramified lifts of the elliptic three-branch quotient. Prescribe
identity returns on all ordinary quotient cycles. This does not assume that
the resulting braid orbit is smooth: every surviving candidate must still pass
the complete six-move closure and independent audit.

Templates are marked quotient actions, not an assertion that all elliptic
quotients or all tuples in a degree have been enumerated:

| Name | beta | gamma | block count |
|---|---|---|---:|
| C3 | (012) | (012) |3|
| C4 | (0123) | (0123) |4|
| S4 | (0123) | (012)(3) |4|

Beta is the second entry on blocks, gamma the third, and delta=beta*gamma.
The quotient indices sum to twice the block count, so its genus is one.
Writing internal permutations A_j and the two remaining voltage variables u,v,
the final return equation is **v X v^-1 = Y**, where composition is rightmost
first and:

| Template | X | Y |
|---|---|---|
| C3 | A2 u^-1 | (A0 u A1)^-1 |
| C4 | A3 u^-1 | (A0 u A2 A1)^-1 |
| S4 | A2 u^-1 | (A0 A3 u A1)^-1 |

The other return equations are eliminated along the quotient incidence tree.
For fixed u there is a solution precisely when X,Y have the same cycle type.
Construct every possible conjugator by matching equal-length cycles and choosing
cyclic offsets. The sampler chooses such matchings/offsets; it does not enumerate
all conjugators or sample uniformly from all source tuples. Compared with
guessing independent u,v and rejecting the root return, this removes a severe
construction rejection. It is an intentional change of proposal law, not a
same-stream micro-optimization.

The conjugator implementation uses bounded stack arrays and stable cycle-root
ordering, with no per-cycle heap allocation. Against the literal cycle-list
version it preserves the resulting permutation and final RNG state on every
S5 pair and 30,000 larger fixtures. Cache use, outer-block closure and the other
sampling lanes remain as in C195. No speculative GPU filter is introduced;
the previous checked Metal path was slower end to end.

## Independent finite coverage and exclusions

With internal degree five, normalize A0 up to conjugacy, enumerate the other
local transpositions, require geodesy, and quotient local frames by simultaneous
S5 conjugacy. The C3 index-two block can be placed at block zero by cyclic block
rotation; its two possible types are a 3-cycle and two transpositions. A pure
full-cycle beta has only a common internal conjugation left after gauge fixing.
This explains the normalization without deleting any parameter in these marked
families. It is not normalization of every genus-three quadruple.

The producer **literally enumerated both remaining S5 voltage permutations**.
A separate Python program independently builds its own S5 multiplication table,
normalizes local frames, solves the independently written return equations,
reconstructs every qualifying source tuple and compares its exact parameter
rank and images against the journal. The native auditor independently checks
every recorded negative braid word.

| Marked family | Source degree | Normalized local frames | Voltage parameters | Qualifying tuples | Non-smooth | Unresolved |
|---|---:|---:|---:|---:|---:|---:|
| C3 |15|10|144,000|2,000|2,000|0|
| C4 |20|25|360,000|5,000|5,000|0|
| S4 |20|25|360,000|5,000|5,000|0|

These are 864,000 parameters and 12,000 qualifying occurrences, **not distinct
components**. They retire just these finite families from the new full-fibre
lane. They do not retire degrees 15/20 from other lanes. Independent coverage
receipts are in FINITE_COVERAGE.json; raw journals remain under
`/Volumes/T7/Research/genus3_full_fiber_20261021/exact_{C3,C4,S4}`.

A separate 50-million-draw degree-24 solved pilot produced 10,513,281
genus-qualified occurrences, 137 final closure calls, no positive component and
no unresolved cap. All 137 persisted late negative words pass the independent
auditor. Earlier arithmetic/short-word rejections in this random pilot are
not claimed independently enumerated. This is negative evidence, not a census.

## A geometric reason to target full fibres

This is a separate sufficient non-isotriviality argument, not a claim that a
smooth orbit has been found. It applies over C, and its analytic reasoning
is written out here for mathematical review; finite tests do not verify it.

**Proposed criterion.** Suppose a product-one transitive quadruple has an
invariant partition into r blocks of size m. For a special entry s, require:

1. Its quotient permutation is the identity, and the quotient curve E has genus
   one. Hence the quotient cover pi:E -> P1 has its three possible branch values
   fixed at 0,1,infinity as the special branch value t moves.
2. The special entry acts nontrivially inside **every** block.
3. Every other entry has relative ramification excess zero.

Then the associated family of source curves is non-isotrivial. The arithmetic
checker additionally restricts to the present source-genus-three target.

**Proof.** After a finite base change the three-branch quotient cover pi is
fixed: its monodromy is finite data with fixed branch values. For t away from
0,1,infinity the relative cover f_t:C_t -> E has no branch value outside
pi^-1(t), by condition 3. Condition 2 makes every point of this fibre a branch
value. Thus its *reduced branch-support divisor* is exactly

    D_t = pi^*[t].

In particular all D_t are linearly equivalent degree-r divisors on fixed E.

For fixed compact C,E and degree m, there are only finitely many maps C -> E
up to translation. To see this, choose a nonzero invariant differential omega
on E. The periods of f^*omega lie in the period lattice of omega. The norm
integral of f^*omega is m times that of omega, so these differentials lie on a
fixed sphere in the finite-dimensional space of holomorphic differentials on C.
The period map is injective and continuous, so its intersection with that
lattice is finite. Maps with the same pullback differential differ by a
translation, since their difference has zero differential. This proves the
fixed-degree assertion without assuming finiteness of all maps to E.

If the source family were isotrivial, locally trivialize its sources as one
fixed C. The normalized maps into E then lie in a finite set, so locally

    f_t = T_{v(t)} composed with f_0,   D_t = T_{v(t)} D_0.

On an elliptic curve, translating a degree-r divisor by v changes its class in
Pic^r(E) by r*v. Because [D_t] is constant, r*v(t) is constant. Multiplication
by r has finite kernel over C, so v(t) is locally constant. This would make
D_t constant, contrary to pi(D_t)={t} with t nonconstant. Contradiction.

The Picard-group identity used in the last step follows from Proposition 4.10
and Remark 4.11(c) in J. S. Milne, *Elliptic Curves*, 2nd edition, printed
pages 32--33: https://www.jmilne.org/math/Books/EC2.pdf. The finite-map argument
also appears in our existing `docs/professor-report/quotient-isotriviality.tex`.
Neither citation asserts that this new combinatorial population contains an
all-smooth orbit. No such existence conclusion is being drawn.

This explains both the targeting and the distinction from the known isotrivial
one-active-block case: that case has one moving branch point, whose divisor
class is not fixed. Here the entire fibre is active, fixing the divisor class
and obstructing the translations that made that earlier family isotrivial.

## Validation, measurement and lifecycle

VALIDATION.json binds the final clean 45-file Rust snapshot, native binary and
independent auditor. The unrelated dirty main.rs is not part of that build.
There are 27 passing relevant Rust tests, 25 passing Python checks and 66
independently audited controls containing 25,192 states; ordinary and macro BFS
agree on their complete state sets. The corruption regression intentionally
prints an audit error while passing. An initial test command named a nonexistent
test module; the corrected command passes and is retained in the test receipt.

BENCHMARK.json records seven randomized-order paired blocks in each of degrees
18 and 24, five million proposals per mode, with the same final binary. The raw
diagnostic flag guesses v; the admitted path solves its conjugacy equation.
The proposal laws differ, so the measured ratio is **genus-qualified candidates
per wall-second**, not same-candidate CPU speed or a discovery multiplier.
Bootstrap intervals describe these timing blocks, not uncertainty about the
existence of the desired example. Exact final measurements are in that receipt.

| Degree | Raw median qualified/sec | Solved median qualified/sec | Mean paired ratio | Bootstrap 95% interval |
|---|---:|---:|---:|---:|
|18|6,386|905,409|141.856|140.047--143.456|
|24|2,259|600,098|267.536|264.331--270.639|

The earlier heap-allocating prototype measured92.132x/171.292x and remains in
the raw benchmark folder. Those are historical intermediate numbers, not the
final deployment measurement. Separate final matched literal/metadata runs have
identical candidate/resolution journals, counters and RNG state.

The previous optimized root finalized stopped after 102 audited jobs,
13,359,279,529 proposal draws and zero new closed/unresolved occurrences.
Its six epoch-25 cursors and learned negative words are preserved. The successor
retains ancestral reporting identities even when the immediate predecessor
found zero new cases; these identities deduplicate reporting, never certify a
positive tuple or justify excluding a new parameter.

The new schedule spends three quarters of slot zero on solved full-fibre
degrees 18/21/24; the remaining quarter retains the prior anchored elliptic
population. The other five slots retain independent block-free, mutation,
rational and alternative cycle-type populations. There is one shared six-slot
CPU pool, no competing old campaign and no automatic publisher. A new closed
case outside the established isotrivial criterion stops the pool for separate
geometry review. A cap/timeout stays unresolved. See LAUNCH.json, PLAN.json and
the internal `.run-control/genus-three-current.json` pointer for the successor;
live STATE is authoritative, not a historical PID or the timing fixtures.

Reproduce finite checks with `python3 scripts/verify_genus_three_fiber.py FOLDER`.
Benchmark with `python3 scripts/benchmark_genus_three_fiber.py --help` only after
stopping the main pool. Source implementation is `src/smooth_genus_three.rs`.

## Software-workflow reference

Experimental-design/performance procedures informed strata, blocked timing
comparisons and explicit evidence limits. Procedural reference only, not a
source for the mathematics: Kassis, Agarwal, He, Patel and Brueckner (2026),
*Scientific Agent Skills*, arXiv:2609.00065, revision 2 (September 2), metadata
checked October 21: https://arxiv.org/abs/2609.00065.
