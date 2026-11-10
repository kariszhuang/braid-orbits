# Geometric versus combinatorial all-smoothness, and the pillowcase picture

October 30, 2026. Mathematics plus an exact small-degree census.

## 1. Two notions of "smooth at a cusp"

Fix a state (g1,g2,g3,g4), product one, transitive, source genus g >= 1, and the cusp {12|34}. The
admissible cover has one source component for each orbit of <g1,g2> (over the bubble) and of <g3,g4>
(over the main component), and one node for each cycle of p = g1 g2. Call this bipartite graph the
**cusp graph**. A vertex v has genus g_v given by Riemann--Hurwitz on its orbit:
2 g_v = |v| + 2 - c_v(x) - c_v(y) - c_v(p).

* **Combinatorial smoothness** (the index equality, used by every engine in this project):
  ind(g1 g2) = ind g1 + ind g2 or ind(g3 g4) = ind g3 + ind g4. Equivalently, one side of the cusp graph
  consists only of genus-0 leaves (one node each). See `../../evidence/smooth_mission_20261011/MATH.md`.
* **Geometric smoothness**: the stable limit of the source curve at this cusp is smooth.

**Lemma G1 (criterion).** Geometric smoothness holds iff the cusp graph is a tree with exactly one
positive-genus vertex.

Proof. Stabilizing the admissible source (forgetting the map) repeatedly contracts genus-0 components with
at most two special points. In a tree with one positive-genus vertex v0, every other vertex lies in a
genus-0 subtree hanging off v0. Such a subtree contracts leaf by leaf, so only v0 survives, smooth. A cycle
in the graph leaves a non-separating node. Two positive-genus vertices leave a separating node. Either
way the stable limit is singular.

Combinatorial smoothness implies geometric smoothness. The other side is then a single vertex (the
transitive-pair lemma), which carries all the genus. **The converse fails.** A genus-0 vertex with
at least two nodes (a "bridge") is allowed geometrically. The core may also lie over the BUBBLE. The
limit curve is then a three-point cover of degree |v0| < d, and the degree-d maps degenerate.

**Index remark.** On the non-core side, Riemann--Hurwitz gives ind x + ind y = d + E - 2V, where E =
cycles(p) and V is the number of vertices on that side. Hence ind x + ind y >= ind p, with equality iff
the side is all leaves. So geometric smoothness does NOT force a "low pair" with index sum <= d - 1 at
each cusp. The STAR/TRIANGLE index-budget dichotomy, and everything built on it, is specific to the
combinatorial notion.

## 2. The gap is real for whole orbits (exact, independently checked)

`geosmooth` (Rust, this folder) enumerates every transitive identity-free four-tuple seed up to
conjugacy. It closes the full six-move orbit and tests Lemma G1 at {12|34} for EVERY state (every cusp of
every state is the {12|34} cusp of some state). Seeds use a class representative for g1, so every orbit
is met; the second version also restricts to a minimal-index first entry and a two-cusp prefilter, both
valid for any all-smooth orbit. `geo_check.py` is an independent Python check (orbit closure from
`../common/braid.py`, and a separate dual-graph test).

| degree d | genus | orbits fully closed | geometrically all-smooth | of which combinatorially all-smooth |
|---|---|---|---|---|
| 4 | all | 7 | 0 | 0 |
| 5 | all | 37 | 0 | 0 |
| 6 | all | 259 | **2 (genus 1, isotrivial)** | **0** |
| 7 | all | 732 | 0 | 0 |
| 8 | 3 | stopped after 18 of 20 g1 classes (all but [2,2,1^4] and [2,1^6]); 0 orbits met (`out/exhaust_d8_g3.log`). Partial, not a census. | | |
| 8 | 2 | not run: genus-2 sources are isotrivial in every degree (M_2 has no complete curves), which covers the geometric notion too | | |

The two degree-6 orbits are geometrically smooth at every cusp, but each contains non-additive states:

    [(1 2 3 6 5 4); (1 3 5); (1 3 5); (1 4)(2 3)(5 6)]   N = 24, 4 combinatorially non-smooth states
    [(1 2 5)(3 4 6); (1 2 5); (1 3)(2 4)(5 6); (1 4)(2 6)(3 5)]   N = 12

Both are confirmed by `geo_check.py`. In the first one, at a bad state, the genus-1 core is the (3,3,3)
cover on three sheets over the bubble. The other side is ONE genus-0 component of degree 6 with four
nodes; it contracts, together with three one-sheet leaves. A genus-1 family that is one-piece everywhere is
isotrivial, so these show the gap; they are not targets. There is no geometrically all-smooth orbit of genus
2 or 3 through degree 7.

## 3. When the two notions agree

**Theorem G2.** Suppose the tuple has a block system whose quotient cover D -> P1 has positive genus and is
itself combinatorially all-smooth. Then the tuple is geometrically all-smooth iff it is combinatorially
all-smooth. This holds in particular when some entry acts trivially on the blocks: such a quotient has at
most three branch values, and its side containing that entry is always additive.

Proof. At a cusp, map each source vertex of C to the vertex of D's cusp graph below it. D is
combinatorially smooth, so one side of D's graph is all genus-0 leaves. The other side is a single vertex
v_D, which carries D's positive genus. Every C-vertex on v_D's side lies over v_D and covers a
positive-genus curve, so it has positive genus. Geometric smoothness allows only one positive-genus vertex.
So C has exactly one vertex on that side, the core. Every C-vertex on the other side meets only the core.
Being a tree, each has exactly one node, and each has genus 0. That side of C is therefore additive. QED

Check (`test_g2.py`, randomized, not a proof). 20,070 cusps were tested along random walks of random
tuples, with m = 3..5 blocks of size 2..4, a positive-genus quotient and a block-trivial fourth entry.
1,162 were geometrically smooth; every one was combinatorially smooth. Control (`control_g2.py`): random
transitive positive-genus tuples WITHOUT the block hypothesis give 166 geometric-only cusps out of
29,295. So the test can detect the gap.

Consequences.
* Every fixed-E (elliptic block quotient with a block-trivial entry) lane, including `elevel`, the October 25
  fixed-E theorems and the paramodular reductions, is unaffected. So are all 1793 catalogued genus-3
  components, which all have such a quotient.
* The geometric relaxation can only matter for primitive groups, or for orbits whose positive-genus
  quotients are not combinatorially smooth (for example rational-quotient towers). These are exactly the
  regions where the combinatorial searches found nothing.
* Purely geometric obstructions (genus 2 / hyperelliptic sources, Theorem A on elliptic factors of exponent
  <= 4, size-two blocks) apply to both notions.

**Open geometric question.** Is there a non-isotrivial genus-3 orbit that is geometrically, but not
combinatorially, all-smooth? It would give a complete curve in M_3 from a four-point Hurwitz family. It
would NOT answer the project's question, which uses the index equality. None exists for d <= 7.

## 4. Pillowcase / Teichmueller dictionary (framework, standard facts)

Let q0 be the pillowcase quadratic differential on P1, with simple poles at the four branch values, and
q = f^* q0. Then (C,q) is a pillowcase cover: a square-tiled half-translation surface made of 2d unit squares.
A pole sits at each fixed point of a g_i, and a zero of order e - 2 at each cycle of length e >= 3. The
four-point Hurwitz curve is the Teichmueller curve of (C,q), with the poles forgotten. The braid action is
the affine (Veech group) action, and the cusps are the rational directions. In direction theta the
pillowcase is one cylinder, whose two ends are saddle connections alpha_theta, beta_theta pairing the four
corners. The cylinders of C are the components of the preimage, one for each cycle of the matching
product.

* Combinatorially all-smooth <=> in every rational direction, the preimage of one end is a forest, with the
  side fixed by the index budget. Equivalently, every cylinder is capped on that common side by a flat
  disk around a tree of saddle connections. Geometric smoothness only asks that every cylinder core bound
  a disk (possibly nested, possibly on either side).
* The capping disks always fill exactly half the area (d of 2d squares). Gauss--Bonnet on each capping disk
  says (#poles) - sum(zero orders) = 2 inside it. Hence each tree contains at least two fixed points,
  which is the inherited fixed-point capacity bound.
* TRIANGLE (high corner infinity) <=> every saddle connection AVOIDING infinity lifts to a forest. Equivalently,
  for every simple arc beta from infinity to a low corner j, (g_inf, g_j along beta) generates a transitive
  genus-3 dessin, namely the cusp core. STAR (centre infinity) <=> every saddle connection THROUGH infinity
  lifts to a forest. Equivalently, every third entry c' in the B3 orbit of (g1,g2,g3) has moved space
  transverse to the FIXED space U_inf = im(1 - g_inf).
* The canonical double cover of (C,q) is an origami with an honest (abelian) Teichmueller curve. In every
  direction its cylinders come in homologous pairs over each capped cylinder, so the + part H^1(C) has
  trivial cusp monodromy.
* Genus-4 combinatorially all-smooth non-isotrivial orbits exist. So any global genus-3 obstruction must
  use something special to genus 3. The candidate is that every complete curve in M_3 meets the
  hyperelliptic divisor. Chen--Moeller-type arguments would need a reason why
  pillowcase fibres of an all-smooth family cannot be hyperelliptic; no such reason is known.

## 5. Relation to earlier work

* The single-moving-point fixed-E conjecture ("SM") of `../frames/README.md` is already covered by
  earlier notes (October 25, elliptic survivor classification), in all degrees.
  The exhaustive `elevel` runs (now k <= 11 with one fixed point, k <= 9 with two and with three, all
  negative) are an independent computational check of those notes, not new coverage.

## Files

* `src/main.rs` - census tool; `out/BINARY_HASHES.txt` binds the source and binary used for degree 8.
* `geo_check.py`, `inspect.py` - independent orbit check and cusp-graph printout.
* `out/exhaust_d7.txt` (first version, all genera), `out/exhaust_d8_g3.*`, `out/exhaust_d8_g2.*`.
* The d <= 6 rows were produced by the first version of `geosmooth` (all genera at once) and reproduced by
  the second version (per genus, minimal-index seeds).
