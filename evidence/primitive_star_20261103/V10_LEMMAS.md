# Twist census v10: in-search lemmas and centralizer symmetry (Nov 4, Claude; parent check only)

Setting as in PLAN.md: STAR normal form (a, b, c, e), abce = 1, a the unique smallest
index, ind b largest, ind e >= ind c. The search builds e one value at a time; c = p^-1 e^-1
with p = ab. v9 already checked T1, T1c, T2b (every j) and T3 during the search.

## 1. Every word gives a hyperforest condition

Track each entry symbolically as g X g^-1 (X in {a,b,c,e}, g a free-group word). A move
sigma_i sends (g_i X_i g_i^-1, g_{i+1} X_{i+1} g_{i+1}^-1) to
(conjugate of X_{i+1} by g_i X_i g_i^-1 g_{i+1}, the old i-th entry); sigma_i^-1 is similar.
After a word w, the cusp (12)|(34) pairs the a-conjugate g a g^-1 with some h Y h^-1, and
pairs the remaining two entries with each other.

- If the other half has index sum >= n, it is never additive (additive needs index of the
  product = sum <= n-1). This is the GATE; each lemma is used only when it holds for every
  allowed ind c.
- Then w fails iff (g a g^-1, h Y h^-1) is not additive, i.e. iff W = h^-1 g sends the cycles
  of a to a non-hyperforest on the Y-cycles. Additivity depends only on cycle partitions,
  and fixed points are leaves of the incidence graph.
- W may be multiplied on the right by powers of a and on the left by powers of Y.

Partial e is sound: (i) undefined points of W(A_i) are dropped, which leaves a
sub-hyperedge; (ii) for Y = e or c the current path components are used, and later arcs
only merge them. Both changes can only remove cycles of the incidence graph, never add
them. So a violation seen on a partial e holds for every completion.

Every derived W below was checked numerically against literal moves (random tuples in S_9,
50 each) in scratch scripts derive.py, fam.py and spec.py.

| lemma | words | W | Y | gate (other half) |
|---|---|---|---|---|
| T2c | 1^(2j),4 (j>0); 4^(2j+1) (j<0) | (ea)^j | c | b,e |
| R3 | 0,1^(2j) ; 0,4^(2j) | (ca)^j (j=1 is T3) | b | c,e |
| R4 | 3,1^(2j) ; 3,4^(2j) | (e^-1 b^-1)^j | b | c,e |
| R5 | 2,1^(2j) ; 2,4^(2j) | (c e^-1 c^-1 b^-1)^j | b | c,e (never pruned; off by default) |
| S131 | 1,3,1 | b c^-1 b^-1 | e | b,c |
| S404 | 4,0,4 | c^-1 b^-1 c | e | b,c |
| S0044 | 0,0,4,4 | a^-1 c a b | b | c,e |
| S0404 | 0,4,0,4 | a b^-1 a^-1 c^-1 | c | b,e |
| S0444 | 0,4,4,4 | c^-1 a^-1 c^-1 | e | b,c |
| S1131 | 1,1,3,1 | b c b^-1 c^-1 b^-1 | e | b,c |
| S1311 | 1,3,1,1 | b^-1 a b c b c^-1 b^-1 | c | b,e |
| S1331 | 1,3,3,1 | c^-1 b^-1 a^-1 b c^-1 b^-1 | b | c,e (never pruned; off) |
| S3111 | 3,1,1,1 | b^-1 a b c b^-1 | e | b,c |

The families were derived for j = +-1, +-2, +-3. For T2c and R4 the symbolic W agrees with
X^j after using abce = 1 (e.g. e(A) = c^-1 b^-1 (A), and c^-1 preserves c-cycles). In
general sigma_2^(2j) conjugates (x_2, x_3) by (x_2 x_3)^j, which gives the j-th power
directly. Families are checked for j up to --jmax (default 64) in both directions; a full
forward period makes the backward pass redundant. The singles are the most frequent short
failing words in the degree-16 census (WORDS files): after T2c (82.3% of all leaf
rejections) come 0,4,4 (8.0%), 1,3,1 (5.3%), 3,1,1 (1.8%) and 4,0,4 (1.3%).

## 2. Centralizer symmetry (exact)

Let H be a subgroup of the centralizer of <a,b> in S_n. For h in H, (a,b,c,e) and
(a, b, h c h^-1, h e h^-1) are simultaneously conjugate, so they lie in the same orbit up to
relabelling. Every test is conjugation-invariant: h commutes with a, b and p, so it permutes
the p-rotations of a.

Order independence: every test is a necessary condition on the completed e, so the set of
leaves below a partial assignment P is {valid e containing P}, whatever order points are
chosen in.

Rule: carry S(P) = {h in H fixing every point used by an arc of P}, a subgroup. At a branch
on the position x (resp. the preimage of a value y), keep only the values y (resp. x) that are
smallest in their orbit under Stab_{S(P)}(x) (resp. (y)). Children carry the pointwise
stabilizer.

Claim: for every valid e containing P there is g in S(P) with g e g^-1 explored below P.
Proof by induction. Let y = e(x) and choose g in Stab(x) with g(y) the orbit minimum. Then
g e g^-1 contains P together with x -> g(y). Induction at that child gives g' in S(child),
a subset of S(P), and g'g lies in S(P). Unexplored values cannot be used points: such a
point is fixed by the stabilizer, so its orbit is itself.

H is built from orbit automorphisms of <a,b> (base point -> t, extended along the
generators and checked) and swaps of isomorphic orbits. A generator whose closure would
exceed --sym-cap (5040) is dropped, so H is always a group. On the 1000 heaviest degree-16
frames the full centralizer would cut leaves by 27.75x (orders up to 256).

## 3. Lemma-driven branching

Each lemma records the unknowns at which its rotations stopped: a position e(x), or a
preimage e^-1(y) (letters E and c). After the T2b walk, the search branches next on the
first such unknown, either over values y or over preimages x. Order does not affect
correctness (section 2). A lemma is re-evaluated only when the new arc fills one of its
recorded unknowns or merges the e- or c-components it reads; otherwise its result is
unchanged.

## 4. Validation (tool-checked)

--brute recomputes, for every frame and every e in S_n, all v9 predicates plus every gated
lemma word replayed literally (tw_lemma_literal, independent of the hyperforest code). It
requires (i) enumerated subset of eligible, and (ii) equal sets of H-canonical forms.
- v10a/b/c/d degree 9 normal form: 190 frames ok, the known genus-2 component found (all 9
  surviving leaves lie in it; v9 had 482 leaves).
- degree 9 all partners: 360 frames ok, same component.
- degree 10 all partners: 1013 frames ok (v10a 39 leaves, v10c 26; v8 870,612).

## 5. Speed (one worker; frame-building setup subtracted)

| frame | v9 | v10 |
|---|---|---|
| d16 616872 (2.0M v9 leaves) | 122 s | 3.1 s (v10d) |
| d16 750258 (4.0M v9 leaves) | 173 s | 5.5 s (v10d) |
| d16 297310 (20M v9 leaves) | not timed | 1.7 s (v10b) |
| d17 1466654 | > 4 h, unfinished when stopped | 60 s (v10b) |

Ablation on 616872 (v10a): no symmetry 21.7M nodes against 8.0M; no lemmas 74.6M nodes and
367,270 leaves against 8.0M and 0; lemma-driven order (v10b) 1.45M nodes.

## 6. Shifted lemmas (v11, opt-in `--shifts`)

sigma1^(2m) conjugates (a,b) by p^m. Relabelling by p^-m turns the state into
(a, b, p^-m c p^m, p^-m e p^m). So every lemma also holds for the shifted pair
(c_m, e_m), with literal word sigma1^(2m) followed by the lemma's word. Example: the T2b walk
at shift m is the rotation family by X_m = (b p^-m c p^m)^-1 on the b-cycles. This was
checked on 3,000 random tuples in S_10, m in [-3,3], j in {+-1,+-2}, with no mismatch
(scratch fam3.py); --brute replays the prefixed words literally.

Effect on two d16 frames with m = +-1, +-2 (tested above depth 0.6 n, and on every leaf):
nodes fall 2.4x and 5.2x, but per-node cost rises about 3.5x. Wall time goes from 8.0 to
10.9 s on 750258 and from 4.2 to 2.8 s on 616872: not a clear gain. So shifts are OFF by
default, and v11d without --shifts reproduces v10d's counts exactly (frame 750258, all
columns).

## 7. Catalan check overflow (fixed in v11d)

The frame builder checks that each p-cycle type yields the product of Catalan numbers.
That expected value was computed through (2l)!, which overflows u128 for l >= 18: the d18
probe stopped at the 18-cycle with expected 83,074 against the correct 477,638,700. v11d
uses the recurrence C(k) = C(k-1)*2(2k-1)/(k+1). Degrees <= 17 were never affected
(34! < 2^128).
v11d validation: --brute d9 normal form, d9 all partners, d10 all partners all pass
(d10 counts equal v10c's); frame 750258 identical to v10d.
