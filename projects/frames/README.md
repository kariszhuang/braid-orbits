# Collision frames of a single moving branch point (October 30, 2026)

Lemma D and the reduction below are new. The experiment is exact for the orbits it reports, but it uses random
samples, not a census.

Code: `src/main.rs` (crate `frames`, `cargo build --release`). Run:
`./target/release/frames --k K --r R --trials T --seed S [--cap N] [--verbose]`. Raw output is in `out/`.

## Question

In the genus-1 quotient lane, consider a SINGLE moving branch point a of C_b -> E (block size k), with fixed
branch points P_1..P_r, r >= 1. The family is all-smooth iff every collision of a with every P_i, along every
arc reached in the braid orbit, is additive. A collision is the "frame" (sigma', tau), and it is additive iff
ind(sigma' tau) = ind sigma' + ind tau, equivalently the cycle hypergraph of the pair is a forest.

* Strong frames conjecture (old README): the frames at P over the orbit are ALL pairs (x sigma x^-1, tau) with
  x in G_E, the monodromy group of C -> E.
* Weak frames conjecture (WF): the frames at P over the orbit include a non-additive pair whenever some
  G_E-conjugate pair is non-additive.

## Lemma D (arc transitivity; proof)

Let S = E minus {P_1..P_r} with r >= 1, and let Push be the image of pi_1(S, a) in Mod(E; a, P_1..P_r)
(Birman exact sequence, valid because chi(S) < 0). Then Push acts transitively on isotopy classes of simple arcs
from a to P_i.

Proof. Let gamma and gamma' be two such arcs. Cutting along either arc leaves homeomorphic surfaces, so some phi
in Mod(E; a, P's) maps gamma to gamma'. The stabiliser of gamma surjects onto Mod(E; P's): collapse gamma into
P_i, and isotope any class to be the identity near P_i. Choose psi in Stab(gamma) with image equal to the inverse
of phi's image. Then phi psi lies in the kernel of the forgetful map, which is Push, and phi psi(gamma) = gamma'.

**Corollary D1.** Fix a state rho of the orbit.
* The frames at P_i over the whole braid orbit are exactly the frames of rho along all simple arcs a -> P_i.
* They are invariants of the full Mod(E; a, P's)-orbit of rho, the "big" Hurwitz component in which E varies too.
  The B_4-orbit used by the search lies between the Push-orbit and the Mod-orbit, so all three have the same frames.

**Corollary D2 (reduction).**
* For r >= 1 and a single moving point, the family is non-isotrivial by de Franchis.
* So "no single-moving all-smooth orbit with r >= 1" is EQUIVALENT to WF for every orbit with frames at
  some P_i.
* If G_E is 2-transitive and sigma, tau are nontrivial, a non-additive G_E-conjugate pair exists: move two points
  of a sigma-cycle onto two points of a tau-cycle.
* By Theorem B, minimal non-isotrivial genus-1 quotients have G_E = A_k or S_k for k >= 9, and for k <= 8 G_E lies
  in a short list of 2-transitive groups plus D_5.
* Hence WF for 2-transitive G_E would close the single-moving-block case in every degree and for every k.
  Exhaustive search (`../quotient_exhaust`) already closes it for k <= 8.

## Experiment

Model: degree 3k tuples on the m = 3 genus-1 dessin. The single-moving-point covers are generated with random
sigma, a random fixed branching tau_2 / tau_3 of total index <= 4, and random handles. The tuple is accepted when
it is transitive, has source genus 3 and has exactly r nontrivial fixed holonomies.

For each accepted tuple:
* Enumerate the full six-move orbit (smooth or not) modulo S_d, keeping the original block labels.
* At every state, for both neighbours of the block-trivial entry s, record the frame (sigma' on block j,
  tau = t^L on block j), up to simultaneous S_k-conjugacy (canonical form valid for intransitive pairs).
* Compare with F_full = {(x sigma' x^-1, tau): x in G_E}, where G_E is the stabiliser of block j restricted to
  block j (Schreier generators).

Results for k = 5, 6 (seeds 51..63, 25 trials per (k, r); for k=6, r=3 only one tuple was found, since that case
needs four transpositions, which Gabai-Kazez already excludes):

| | orbits | weak WF holds | strong conjecture holds | primitive G_E |
|---|---|---|---|---|
| k=5, r=1 | 25 | 25 | 24 | 25 |
| k=5, r=2 | 25 | 25 | 25 | 25 |
| k=5, r=3 | 25 | 25 | 25 | 25 |
| k=6, r=1 | 25 | 25 | 17 | 25 |
| k=6, r=2 | 25 | 25 | 25 | 24 |
| k=6, r=3 | 1 | 1 | 1 | 1 |
| k=7, r=1 (seed 701) | 10 | 10 | 8 | 10 |

The two k = 7 strong failures again have r = 1 and sigma, tau both 2^2 (G_E = A_7, N = 38592 and 11136).

All 126 orbits contain non-additive frames (none is all-smooth). The strong conjecture FAILS in 9 orbits, all
with r = 1:
* k = 5, G_E = A_5, sigma and tau both 3-cycles: the orbit (N = 1620) contains the cancelling frame (tau^-1, tau)
  but misses (tau, tau) and the frame with product of type 2^2.
* k = 6, sigma and tau both 2^2 (G_E = S_6, N = 8640; and G_E = A_6, N = 3168): the orbit misses the frames whose
  product has type (4,2).
* k = 6, sigma and tau both 3-cycles (G_E = S_6, N = 12960, 5 times): 4 of 6 frames observed.

All 9 (k=5,6) failures, and the 2 at k=7, have r = 1 AND sigma, tau of the same cycle type. No failure with r >= 2 or with different types.
The one imprimitive G_E in the sample (k = 6, r = 2, |G_E| = 24) satisfies both forms.

Interpretation. With two branch points the big Hurwitz space can be disconnected, with components separated by
lifting-invariant type data, and a component then sees only part of the double cosets. So the right statement
is the weak one. In every case it holds, and the cancelling frame (tau^-1, tau) appears whenever sigma is
conjugate to tau^-1.

Larger k: see `out/k7_*`, `out/k8_*`, `out/k9_*` (sweep2.sh; caps reported as CAPPED, never as evidence).

## Open

Prove WF for 2-transitive G_E. Two possible routes:
1. Degenerate E inside the big Hurwitz component (pinch a handle curve). The a-P frames then contain those of a
   genus-0 four-point problem (sigma, tau, A, A^-1-conjugate) for every nonseparating curve alpha.
2. Show that every component of the two-point Hurwitz space of genus-1 covers contains a cancelling or
   overlapping degeneration.

## E-level decider for a single moving point (October 30, evening)

The criterion is new. The computations are exact and exhaustive in the scopes stated.

**Criterion E.** Write the cover C_a -> E by its monodromy on the n-punctured torus: handles A, B and branch
elements c_1..c_n with [A,B] = A B A^-1 B^-1 = c_1 c_2 ... c_n in S_k. One c is the moving point a, the others are
the fixed points P_i (r = n - 1 <= 3, since the total index is 4).
* By Lemma D, the collisions of a over the point-pushing orbit are exactly those over the whole mapping-class-group
  orbit.
* Hence the family is all-smooth iff in every state of that orbit, every adjacent pair (c_i, c_{i+1}) containing a
  is additive: ind(c_i c_{i+1}) = ind c_i + ind c_{i+1}. That is, the preimage of the two-point disk is a union of
  disks.
* For r >= 1 any all-smooth orbit is automatically NON-isotrivial (de Franchis). So a single surviving orbit
  would be a target.
* For r = 1 the condition reads: ind [A,B] = 4 in every state. Equivalently, for EVERY embedded one-holed torus
  T in E minus {a, P}, the cover restricted to T is a connected genus-3 origami (the one-point cover obtained by
  sliding a into P never loses genus).

Moves used (all genuine mapping classes; automorphisms of pi_1 preserving the peripheral classes):
* HA: B -> BA, and HB: A -> AB (handle twists).
* PC: (A, B, c_1, c_j) -> (A, c_1 B, A c_1 A^-1, c_1 c_j c_1^-1), which pushes puncture 1 through a handle.
* BR_i: braids of the punctures.
* All of these and their inverses. Composing PC n times gives the frame ((s^-1 A)^n s (s^-1 A)^-n, t), where s and t
  are the branch elements c_1, c_2 (sigma and tau earlier), so "twist the handles, then push" yields conjugators
  that are powers of s^-1 A' for primitive handles A'.

Rigour:
* A non-additive state reached by literal moves is a rigorous negative.
* An orbit with no non-additive state is only a CANDIDATE. The generators are not proved to contain the whole
  push group, so a candidate must be confirmed by the degree-3k engine. No candidate ever occurred in genus 3.
* Starts: sigma is fixed to a class representative, tau_1 runs over C(sigma)-orbit representatives with
  (sigma, tau_1) additive, the other taus over their whole classes, A over all of S_k, and B over all solutions of
  B A^-1 B^-1 = A^-1 g (one conjugator times the centraliser). Every state of an all-smooth orbit is conjugate to
  one of these starts, so the enumeration is complete. ALL transitive G_E are included, imprimitive ones too.

Code: `src/bin/elevel.rs` (modes `--exhaust-r1`, `--exhaust-r2`, `--exhaust-gen`, `--control`, `--depth`,
`--push-lemma`); scripts `r1_sweep.sh`, `r2_sweep.sh`, `queue_*.sh`; raw output `elevel_out/`.
* Binaries were built from successive versions of the same file. The move, canonical-form and BFS code is
  identical across versions; later versions only add modes and raise the array width K from 16 to 32.
* SHA-256: target/ e9d938f2..., target2/ 3cb42341..., target3/ d4c05548..., current source 197b53d7....

Validation:
* Positive control (genus 4, the C -> D -> E mechanism with D of genus 2 constant): in `--control`, 100/100 tuples
  give all-additive orbits. Exhaustive `--exhaust-r2 --k 6 --sig 3,3 --tau1 2 --tau2 2` finds 432 all-additive
  starts, all imprimitive, in orbits of 1296 states, and refutes all 9504 primitive starts. `--exhaust-gen`
  reproduces these numbers exactly.
* Consistency: for k <= 8 the decider finds no genus-3 survivor, in agreement with the independent degree-3k engine
  (single moving block closed for k <= 8).

**Results (genus 3, exhaustive, all transitive groups):**

| r (fixed points) | class patterns | k | starts | negative | all-additive |
|---|---|---|---|---|---|
| 1 | all 10 pairs (sigma, tau) of total index 4 | 4..8 | 373,391 | all | 0 |
| 1 | all 10 pairs | 9 | 867,616 | all | 0 |
| 1 | all 10 pairs | 10 | 2,302,102 | all | 0 |
| 2 | all 6 triples (index 2,1,1 in every order) | 4..7 | 459,254 | all | 0 |
| 2 | all 6 triples | 8 | 1,187,178 | all | 0 |
| 2 | all 6 triples | 9 | 3,350,973 | all | 0 |
| 3 | four transpositions | 4..9 | 4,210,224 (k=9: 2,768,256) | all | 0 |
| 1 | all 10 pairs | 11 | 5,455,290 | all | 0 |

So, independently of the engine, NO all-smooth genus-3 family with one moving branch point over a fixed elliptic
curve exists for block size k <= 11 with one fixed branch point, k <= 9 with two, and k <= 9 with three. Refutation is always fast: at most 171 BFS states.

**Large k (random, not exhaustive).** `--depth` samples random additive starts and reports the exact number of
moves to the first non-additive state:
* k = 8, 12, 14: 150 samples per pattern, every one refuted within 6 moves.
* k = 20, 30: 80 samples for each of the 11 patterns (r = 1, 2, 3), 1,760 samples, every one refuted within 9 moves.
* The depth grows slowly with k (typical 3-5 at k = 30).

**Correction (October 30, later the same day).** Conjecture SM is not new. Earlier proof notes (October 25,
elliptic survivor classification) already exclude every one-moving-value fixed-E
configuration in ALL degrees: a moving value with a single simple point by the ONE_MOVING_SIMPLE_TARGET
(Gabai--Kazez slide) theorem; zero profiles (4), (3,1), (2,2) directly; (3)/(2^2) by the revised theta proof;
a moving order-two zero with two fixed simple values by an earlier collision-graph theorem; a moving (2^2) with fixed (2),(2) by the
four-simple separation theorem; and (2^2)/(2^2) by the synchronized-pair theorem. The exhaustive runs above are therefore an INDEPENDENT
COMPUTATIONAL CHECK of those notes for small block size, not new coverage of the global question. They agree
with every case, including the synchronized (2,2)/(2,2) case, which is negative for every k <= 11.
The genuinely open fixed-E cases have at least two moving values (multi-moving blocks).

Conjecture SM (single moving point): for every k, no all-smooth genus-3 family has exactly one moving branch point
over a fixed elliptic curve. By Corollary D2 this is equivalent to WF. The evidence is the exhaustive k <= 10
(r = 1), the k <= 8 engine (all r) and the random samples to k = 30.

A proof could go through the one-holed-torus form: show that some embedded one-holed torus in E minus {a, P}
carries a disconnected, or lower-genus, restricted cover. In the equal-cycle-type case, an arc giving the
cancelling frame (tau^-1, tau) is the same thing as an embedded one-holed torus whose handle images COMMUTE.

## Genus-0 quotient lane, single moving point (P^1-level decider, `--p1`)

Setting: C -> D = P^1 of degree k over a FIXED rational quotient (a genus-0 dessin of degree m), with r fixed
branch points (dessin vertices, so r <= m + 2) and ONE moving point a. Lemma D holds for P^1 minus r >= 3 points,
so the collisions are invariants of the full braid orbit of (c_a, c_1, ..., c_r), with c_a c_1 ... c_r = 1. A
collision of a with P_i is smooth iff one side of the pinched circle is all disks: either the pair is additive,
or the complementary disk is.

**Index lemma (proof).** For source genus 3 the total index is 2k + 4. If the pair side were not additive, the
complementary side would have to be, so its total index would be at most the index of the boundary element
(at most k - 1). The complement has index 2k + 4 - (ind c_a + ind c_i), which is at least k + 5 because
ind c_a + ind c_i <= 2(k - 1). That is impossible.

Hence every collision is smooth only on the PAIR side, and ind c_a + ind c_i <= k - 1 for every fixed point P_i
reached. With s = ind c_a this forces r >= (2k + 4 - s)/(k - 1 - s): at least 5 fixed points for k = 5, at least 4
for k = 6..8, and r = 3 only for k >= 9. That last case is the original four-point problem, where the known
bound d >= 2g + 3 = 9 applies.

The genus-0 quotient lane therefore needs either large block size k or many fixed branch points (r >= 4,
i.e. dessins of degree m >= 2). The degree-3k engine never covered it beyond k = 4, and skipped the m = 2 quotient
(z = id). Unlike the genus-1 lane, survivors would NOT be automatically non-isotrivial: isotriviality
must be tested on the lifted four-point orbit.

Validation of `--p1` (`--p1-check`): the engine's three all-smooth controls are recognised exactly. They are
degree 9 genus 2 (N = 12), degree 20 genus 6 (N = 240) and degree 24 genus 4 (N = 192). The counts are for moving
entry 1 and equal the engine's N.

Random index-constrained samples (`p1_sweep.sh`, `elevel_out/p1_random.txt`; not a census): k = 6, 7, 8 with r = 4, 5
and k = 9 with r = 3, 4 gave 1,600 smooth starts, ALL refuted, 0 survivors. Further (k, r) pairs are appended by
the running sweep.
