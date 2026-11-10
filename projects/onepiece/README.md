# One-piece orbits and bad clusters (November 2026)

Three questions: (1) search for all-one-piece orbits
instead of all-smooth ones; (2) why most catalogue examples have a genus-2
intermediate cover; (3) in genus 2, how few non-one-piece clusters a non-isotrivial
orbit can have. Code: `src/main.rs` (crate `onepiece`); independent check: `check.py`.
Raw runs: `data/one_piece_20261105` (binaries and sources hashed in
HASHES.txt). Small-degree outputs are in `out/`.

## Conventions

Composition rightmost-first, `(x*y)[i] = x[y[i]]`, `g1 g2 g3 g4 = 1` (`--orbit`
also accepts the other order and reverses the tuple). States are
taken modulo simultaneous conjugation in S_d. Moves: sigma_j (a,b) -> (a b a^-1, a).

## 1. Testing one-piece at a state (O(d))

At the cusp {12|34} of (g1,g2,g3,g4), with p = g1 g2:

1. E = number of cycles of p (the nodes).
2. V12 = number of orbits of <g1,g2>, V34 = number of orbits of <g3,g4>.
3. Defect = E - V12 - V34 + 1. The dual graph is connected (transitivity), so it is a
   tree iff defect = 0. Most states fail here, and this needs no genus computation.
4. For each orbit O, 2 g_O = |O| + 2 - c_O(x) - c_O(y) - c_O(p), where c_O counts cycles
   inside O and (x,y) is the pair on that side.
5. One-piece iff defect = 0 and at most one g_O > 0.

Equivalent to Lemma G1 of `../geosmooth` (stable limit smooth). Being
combinatorially smooth (the index equality) implies one-piece, but not
conversely.

**Clusters.** sigma1 and sigma3 do not change <g1,g2>, <g3,g4>, p, or the cycle types
of g1, g2 on each orbit. So every invariant of the {12|34} cusp, including being
one-piece, is constant on the orbits of <sigma1, sigma3>, called clusters. The code
asserts this on every orbit it reports. Geometrically, a cluster is a point of the
compactified Hurwitz curve over the boundary of M_{0,4}. A bad cluster is a point
where the curve meets the boundary of M_g-bar, i.e. a singular fibre.

## 2. The census algorithm (`onepiece --census --d D --genus G`)

Goal: every transitive, identity-free orbit of degree D and genus G that is one-piece at
every state.

* **Normal form.** Any entry can be moved to position 1 by braid moves without changing
  its conjugacy class. So every orbit has a state where g1 has the smallest index of the
  four entries; conjugating, g1 is a fixed class representative. Conjugating by the
  centralizer C(g1), g2 runs over C(g1)-orbit representatives only. g3 runs over all of
  S_D, and g4 = (g1 g2 g3)^-1. Every orbit is met.
* **Filters, cheapest first** (each is a necessary condition for all-one-piece):
  - index budget: ind g1 <= the other three indices, and the indices sum to 2D - 2 + 2G;
  - <g1,g2> has at most one positive-genus orbit (checked once per g2);
  - tree test without computing g4: <g3,g4> = <g3,p>, so its orbits are the
    components of the cycles of p joined by g3. A union-find over the E cycles
    gives V34, which must equal E - V12 + 1;
  - then g4, transitivity, and the full one-piece test at four cusps: {12|34}, {23|41}
    (the rotation), and {12|34} of sigma2^(+-1) of the state.
* **Closure.** Survivors are looked up in a shared hash set of canonical states;
  unseen ones get a breadth-first six-move closure. The closure stops at the first
  state that is not one-piece. All visited states are stored, so no orbit is explored
  twice. Canonical form: relabel by BFS from base points, trying only base points with
  the least (g1,g2,g3)-cycle-length signature.
* **Orbit report.** Clusters, bad clusters (split into defect > 0 and two
  positive-genus pieces), bad-cluster sizes, smooth-state count, primitivity and all
  block sizes, and the boundary counts delta_0, delta_1 of the stable model.
* **Other modes.** `--all-orbits` closes every orbit of the given degree and genus,
  with no smoothness filter. `--all-orbits --need-good` closes every orbit with at least
  one one-piece cluster. This is exhaustive for that class, because the smallest-index
  class can be brought to position 1 inside any cluster using sigma1 and rotation by
  two. `--orbit "a;b;c;d" --clusters` prints one orbit cluster by cluster.

**Stable model and delta.** Each node over a cycle of length e of p has thickness 1/e.
Stabilize by deleting genus-0 leaves and merging genus-0 bivalent pieces (thicknesses
add). Remaining edges are non-separating (delta_0) or separating (delta_i). Summing
over states and the three cusp types gives delta on the Hurwitz curve.

**Validation (tool-checked).**
* O9, O10, A, B reproduce Nitiwit's statements. O9 and O10 are one-piece everywhere
  (O10 is smooth at 8 of 12 states); A and B are not.
* Degrees 4-7: the all-one-piece census reproduces `geosmooth` (the two degree-6
  genus-1 orbits).
* `check.py` (brute force over S_d^3, union-find, no canonical forms, separate code)
  gives the same orbits and the same (clusters, bad clusters) histogram at d = 4, 5,
  genus 2.
* `--need-good` reproduces exactly the subset of the full runs at d = 5, 6, 7.

## 3. Results

### 3a. All-one-piece orbits (exhaustive)

| degree | genus 1 | genus 2 | genus >= 3 |
|---|---|---|---|
| 4, 5, 7 | 0 | 0 | 0 |
| 6 | 2 | 0 | 0 |
| 8 | 9 | 0 | 0 (genus 3..7 all run) |
| 9 | 18 (all imprimitive, blocks of size 3) | **1 (O9)** | 0 (genus 3..8 all run) |
| 10 | not run | **6** (all imprimitive, blocks of size 5, isotrivial; includes O10) | 0 (genus 3..9 all run) |

All 11 orbits are imprimitive (blocks of size 2, 3 or 4) and of genus 1, hence
isotrivial. One degree-8 orbit (4 states, types 3^2, 3^2, 2^4, 3^2) is one-piece
everywhere but has NO combinatorially smooth state. So "one-piece everywhere" does not
imply "smooth somewhere", and a census cannot be seeded only from smooth states.
Degree 8 genus 3 completes the run that was stopped on Oct 30.

### 3b. Genus 2: how few bad clusters (exhaustive, every orbit)

| degree | genus-2 orbits | one-piece everywhere | fewest bad clusters | example |
|---|---|---|---|---|
| 4 | 6 | 0 | 3 | (12)(34), (1324), (12)(34), (1423): 3 of 3 bad |
| 5 | 38 | 0 | **1** | (23)(45), (12)(34), (13542), (14325): group D5, 1 of 3 bad |
| 6 | 209 | 0 | 2 | (23)(56), (12465), (12534), (12534): 2 of 4 bad |
| 7 | 452 | 0 | 6 | e.g. (2346)(57), (12)(3576), (24)(3567), (1453)(26) |
| 8 | 2,376 | 0 | 3 | (3 4)(7 8), (1 2)(5 6), (1 3 5 7 2 4 6 8), (1 7 6 3 2 8 5 4): 3 of 5 bad |
| 9 | 4,513 with a good cluster (+1 capped) | 1 (O9) | 3 | 18 states, 3 of 5 bad (two such orbits); v5 need-good census, see section 7 |
| 10 | (sampled only) | | **1** | (5 8)(6 9), (2 9)(5 6), (1 6 3 9 10 8 4 2 7 5), (1 5 10 9 7 2 3 6 4 8): 1 of 5 bad, blocks of size 5 |

Degree 9 genus 2: the all-one-piece census finds exactly ONE orbit, O9. Degree 9 genus 3: none.

Random sampling of the same shape (two double transpositions, two full cycles; script
`sample_family.py`, orbits capped at 20,000 states) found at best 8 bad clusters at
d = 12 (1,500 tries) and 22 at d = 15. This is biased sampling, not a census: the
one-bad-cluster orbits are special, and random tuples are generic.

So the count does not grow with the degree in this range: degrees 5 and 10 both have
an orbit with a single bad cluster (one boundary point of M_2-bar). In that orbit the bad cluster has 8
of the 12 states and defect 2, i.e. two non-separating nodes. Geometrically, every
non-isotrivial genus-2 orbit has at least one bad cluster, since M_2 contains no complete
curves. Whether the minimum grows with d is open; the data so far do not suggest it.

## 4. Why a genus-2 intermediate cover (catalogue v17)

Among the 16,449 catalogue components that are non-isotrivial or unknown, 15,485
(94%) have a recorded genus-2 block quotient. Among the 240,111 isotrivial ones, only
86 do; their smallest quotient is almost always genus 1.

Mechanism: in every case the genus-2 quotient curve D is FIXED along the family, for
one of two reasons.
* 6,362 of the 15,485 have a THREE-point genus-2 quotient: the special entry acts
  inside the blocks, so it is the identity on the quotient (index patterns (0,5,5,6),
  (0,5,6,7), (0,4,5,5), ...). D is then a Belyi curve, which is rigid.
* For the rest, smoothness descends to block quotients, so D is an all-smooth genus-2
  family. All-smooth genus 2 is isotrivial (M_2 contains no complete curves), so D is
  again fixed.

The source is then a cover of a fixed curve D, branched over points that move.
* If D has genus 1, translations of D move the branch point back, so the family is
  isotrivial (Proposition 1). These are the 240k isotrivial entries.
* If D has genus >= 2, Aut(D) is finite, so the moving branch point gives a
  non-isotrivial family (the fixed-quotient theorem).
* Genus 2 is the smallest genus that does this, and genus-2 Belyi curves already occur
  in degrees 5-10. Higher-genus bases need larger degree, which the degree-bounded
  searches reach less often. Some lanes also lift the degree-9 genus-2 orbit O9 by
  design.

## 5. Speed-ups (v3/v4, Nov 5 late; `census_fast` in src/main.rs)

Profile of v1: 99% of the time was the loop over all d! choices of g3 (about 40 ns
each); only ~6% of g3 passed the cheap tests.

* **Half-cusp tables (the main gain).** Whether the {12|34} cusp is one-piece depends on
  the (g3,g4) half only through p = g1 g2. For each cycle type mu of p, once:
  every g3 with g4 = (p_mu g3)^-1 whose side <g3,g4> has at most one positive-genus
  orbit and total genus 0 or G, bucketed by (mu, number of orbits, genus > 0) and
  sorted by min(ind g3, ind g4). Each item (g1,g2) is conjugated so that p = p_mu;
  the tree condition E = V12 + V34 - 1 picks the bucket, and the min-index condition
  is a prefix of it. The index budget then holds automatically.
* **New prune.** An item whose <g1,g2> side has genus strictly between 0 and G is
  impossible (two positive pieces) and is skipped.
* **Orbit ownership.** Concurrent closures of one orbit are resolved by closure id:
  the larger id aborts on meeting a state owned by a smaller id, so each orbit is closed
  once (d7 need-good: 17.6M closed states fell to 5.2M).
* **Canonical forms** try only base points in the rarest (g1,g2,g3,g4)-cycle-length
  class. **Reports** compute the {12|34} cusp once per cluster and the {23|41} cusp once
  per <sigma2>-orbit. The one-piece test is allocation-free and runs before canonicalizing.

Validation: v3/v4 output equals v1/v2 output (orbits, states, clusters, bad clusters,
cluster sizes, delta, smooth states, defect states, primitivity, block sizes) on
all 22 one-piece runs d4..8, the need-good runs d5..7 and the every-orbit runs d5..7;
d9 genus 2 (O9) and genus 3 (none) agree with v1.

| run | v1 (8 threads) | v4 (10 threads) |
|---|---|---|
| d9 genus 3 | 1,201 s | 32 s |
| d9 genus 2 | 932 s | 115 s |
| d9 genus 5 | about 25 min (estimate, stopped) | 2 s |
| d7 genus 2 need-good | 29 s | 8.4 s |
| d8 genus 1 | 9.6 s | 11.5 s (no gain: seed canonical forms dominate) |

Binaries: data/one_piece_20261105/onepiece_v4 (HASHES.txt; its summary
lines say "v3"). Chain run_v4.sh: d9 need-good genus 2, then d10 genus 3, 2, 4..9, then
d10 need-good genus 2.

## 6. Nov 6 status

* Degree 10 all-one-piece census (v4) complete for genus 2..9: genus 2 has 6 orbits
  (12, 12, 24, 24, 24, 48 states), all imprimitive with blocks of size 5 and isotrivial,
  including O10; genus 3..9 have none. So through degree 10, every orbit of genus >= 2
  that is one-piece everywhere has genus 2, is isotrivial and is imprimitive.
* The v4 degree-9 need-good run (the genus-2 bad-cluster question) ended WITHOUT a summary
  line after about 14 h, most likely from memory (v4 stored every state of every
  closed orbit). Its 774 orbits (fewest bad clusters 4) are partial and O9 is missing
  from them, so they are not a census. Kept as good_g2_d9_v4_INCOMPLETE.jsonl. The
  chain's "exit=0" was a script bug ($? read after $(date)).
* The v4 chain was killed at about 22:22 Nov 6 when its parent application restarted (it
  was a child process of the app). good_g2_d10 was partial; kept as good_g2_d10_v4_KILLED.jsonl.
* v5: the shared table keeps only states that could be seeds (g1 of minimal index and
  {12|34} one-piece): 160k instead of 7.0M stored states at d7, 371 MB peak RSS, d7
  need-good 4.9 s. Output identical on d6/d7 need-good and d8 genus 1..3.
* run_v5.sh runs in its own session with launchd as parent (survives app restarts):
  good_g2_d9, then degree 11 genus 3, 2, 4..10.

## 7. Nov 7: degree-9 genus-2 cluster census (v5, complete)

* `good_g2_d9.jsonl` (T7), summary complete:true, 5,273 s on 10 threads, 106,299
  closures, 137M stored seed-like states. Scope: every genus-2 orbit with AT LEAST ONE
  one-piece cluster. 4,513 such orbits; ONE orbit hit the 20M-state cap and is
  unresolved (tuple recorded in the file).
* Fewest bad clusters among non-isotrivial orbits: 3 (two orbits, 18 states, 5
  clusters each; one has 1 defect + 2 separating bad clusters, the other 3 defect).
  Next: 4 (23 orbits). Only O9 has 0. Sequence d4..9: 3, 1, 2, 6, 3, 3.
* Orbits with NO good cluster are outside this mode. In d4..8 every such orbit has at
  least 3 clusters, so none beats the minimum there; at d9 this is not checked.
* 4,216 of the 4,513 are primitive (all non-isotrivial).

## 8. Nov 7: one-piece descent and the genus-1 lift search

### 8.1 Lemma (one-piece descends to block quotients)

SUPERSEDED by DESCENT.md (Nov 7): full proof by cycle counting (positive pieces) and
surface/graph homology (defect), no admissible covers; stronger monotone form; exhaustive
checks over S2 wr S3 and S3 wr S2 plus 30,000 random tuples. The sketch below is kept.

Let a transitive tuple preserve a block system, with quotient tuple on the blocks (the
intermediate cover X -> Z -> P1). If a state is one-piece at {12|34}, so is its quotient.
Hence the quotient of an orbit that is one-piece everywhere is one-piece everywhere.

Sketch. At the cusp the admissible cover gives a finite surjective map X0 -> Z0 of nodal
curves whose components are the pieces (orbits of <g1,g2>, <g3,g4> on letters, resp.
blocks) and whose nodes are the cycles of g1g2. (i) Tree: the defect is the first Betti
number of the dual graph = toric rank of Pic^0. Pullback Pic^0(Z0) -> Pic^0(X0) has
finite kernel (norm o pullback = multiplication by the degree), and a torus has no
nonzero map to the compact Pic^0(X0) of a tree-like curve; so toric rank(Z0) = 0.
(ii) A genus-0 piece only maps onto genus-0 pieces, so every positive-genus piece of Z0
is the image of the unique positive-genus piece of X0: at most one.
This is different from the older all-smooth "quotient descent" (pair additivity). Random check: 12,385 random block-preserving
tuples (S_2 or S_3 wr S_3..5) that are one-piece at {12|34}; every quotient is one-piece.

### 8.2 Consequence for genus 3 (block size m, quotient genus h)

h = 3 forces m = 1. h = 2: the quotient orbit is genus 2 and one-piece everywhere, hence
isotrivial (no complete curves in M_2), and Riemann-Hurwitz 4 = 2m - 4 + R gives m = 2,
R = 0: an unramified double cover of a fixed curve, isotrivial. h = 1: the quotient is a
genus-1 orbit that is one-piece everywhere (the census lists all of them up to degree 9:
2 at d6, 9 at d8, 18 at d9; NONE is all-smooth, so all lie outside the older elliptic
lanes, which used quotients with an identity entry), the elliptic curve is fixed, and
the returns have total index 4. Theorem A (the paramodular obstruction) predicts m >= 5 for a
non-isotrivial example. h = 0: no constraint from the lemma.

### 8.3 Lift search (src/bin/lift.rs)

Over a fixed base state, a lift is given by S_m voltages on the 3b edges x -> a_i x;
spanning-tree voltages are gauged to 1. The return of each base cycle has index e_c and
sum e_c = 2g - 2 - m(2h - 2). A pattern prescribes the cycles with nontrivial return and
their cycle types. Tietze elimination solves every trivial-return relation; designated
ones get a class-valued variable. For every genus-1 base and pattern this leaves exactly
two free permutations (A, B) and one class condition on the last designated return. The
first free one runs over class representatives (global conjugation); the second is
solved as a conjugator (a coset of a centralizer) for each element t of the prescribed
class, instead of being enumerated over S_m. Each lift is then checked directly on the
permutations (genus, transitivity, return types by assertion), one-piece at the three
cusps of the state, and closed with early exit; closed orbits are de-duplicated.

Checks:
* O9 (degree 9, genus 2), m = 2, genus 4: exactly the four known degree-18 genus-4
  all-smooth components (1,296 states each, non-isotrivial) found earlier by the
  all-smooth search. Tietze and brute force (all 2^19 voltage
  assignments) give the same 432 seeds.
* All 29 genus-1 bases, m = 2: Tietze = brute force (identical seed sets). m = 3 on
  four bases: plain enumeration = conjugator solver (identical seed sets).
* m = 2 gives 26 genus-3 orbits that are one-piece everywhere (degree 16, over the
  degree-8 bases), all isotrivial, as required for degree-2 maps to a fixed elliptic
  curve. m = 3 gives 30/682/30/30 such orbits over B00/B04/B11/B20, all isotrivial.
* m = 4 on B00: conjugator solver = plain enumeration (identical seed sets, 621 orbits).
* lift_v2 adds a cache of states visited by failed closures (128-bit hashes); B00 m = 4
  376 s -> 49 s, same 621 orbits and seeds.

### 8.4 Runs (T7 data/one_piece_20261105/lift, RUN.log, HASHES.txt)

* m = 3, all 29 bases (lift_v1): 3,366 genus-3 orbits one-piece everywhere (degrees 18,
  24, 27), every one isotrivial. No all-smooth condition is imposed.
* m = 4 then m = 5, all bases (lift_v2, 4 threads, alongside the d11 census): running.
* Orbit-level de-duplication (`lift --dedupe FILE...`): every recorded orbit is re-closed
  from its tuple (asserting one-piece everywhere, no cap) and identified by its least
  canonical state modulo simultaneous S_d relabeling. m = 2: 26 distinct of 26. m = 3:
  3,366 of 3,366. m = 4: 67,985 distinct of 67,991; the 6 repeats are degree-32 orbits
  found over two different degree-8 bases (two genus-1 block quotients). m = 5 (B00..B02):
  26,491 of 26,491. Degree-24 cross-check (m = 4 over d6 vs m = 3 over d8): no overlap.
  No duplicates inside any file.
* Genus-2 bases (O9 and the six degree-10 orbits), m = 2, genus 3 (unramified double
  lifts, excess 0; out/lift/G*_m2_g3.jsonl): O9 gives 15 orbits at degree 18 (12 states,
  ALL states smooth: all-smooth), the d10 orbits give 30 at degree 20 (none all-smooth).
  Note: section 8.2 excludes quotient genus 2 only for NON-isotrivial targets; as
  combinatorial objects, one-piece genus-3 orbits with genus-2 quotients exist.

## 9. Nov 8: speedups (both searches switched; outputs unchanged)

* onepiece_v6 (census): NEON one_piece for M = 16 (composition = one `tbl`, orbit/cycle
  labels by min-and-jump). `--selftest`: 9,000,000 states agree with the scalar test.
  Census d9 genus 3: 151 s (v5) -> 57 s (v6), identical scanned/seeds/passed counts. Exact
  output matches (cmp.py) on d6 g1, d8 g1, d9 g1, d9 g2 (O9) and d7 genus-2 need-good (312
  orbits). New `--start-item S` resumes a run (item list is deterministic). d11 genus 3 was
  switched: v5 reached printed item 11,200,000 of 23,466,359 (found 0, capped 0; kept as
  op_d11_g3_v5_part.*), v6 resumes at 11,000,000 (op_d11_g3_resume.*), then d11 genus 2.
* lift_v3: allocation-free one_piece / canon / move (`--selftest`: 250,935 tuples, 1.5M moves
  identical to the reference code), a 64-step random-walk probe before each closure (only
  rejects on a state that is not one-piece), sampled assertions. B03 m = 4: 403 s -> 193 s,
  identical 7,165 orbits (--dedupe). On the live m = 5 run about 4-5x faster per pattern.
