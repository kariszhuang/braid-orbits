# Descent of one-piece cusps to block quotients

Nov 7, 2026. Proof written out in full; exhaustive and random checks in
`descent_check.py` (results in `out/descent_check.txt`).
This replaces the sketch in README section 8.1 (which used admissible covers and
generalized Jacobians); the proof below uses only permutation counting and the
homology of surfaces and graphs.

## Setting

Let `t = (g1, g2, g3, g4)` be permutations of a finite set `Omega`, with
`g1 g2 g3 g4 = 1` (composition rightmost-first) and `<g1,...,g4>` transitive. Suppose
`G = <g1,...,g4>` preserves a system of blocks of size `m`, and write `abar` for the
induced permutation of the set of blocks `Omegabar`. Then `tbar = (g1bar, ..., g4bar)` is
again a transitive product-one tuple. Entries of `tbar` may be the identity.

Fix the cusp `{12|34}`. Write `x = g1`, `y = g2`, `p = xy`, `z = g3`, `w = g4`
(so `zw = p^-1`). The **cusp graph** `Gamma(t)` has one vertex for each orbit of `<x,y>`
(side A) and of `<z,w>` (side B), and one edge for each cycle of `p`. That edge joins the
A-orbit and the B-orbit containing the cycle. Since `G` is transitive, `Gamma` is connected.
Its first Betti number is the **defect**:

    defect(t) = b1(Gamma) = cycles(p) - orbits<x,y> - orbits<z,w> + 1.

The **genus of a piece** `O` (an orbit on side A, say) is given by

    2 g_O - 2 = |O| - c_O(x) - c_O(y) - c_O(p),

where `c_O(a)` is the number of cycles of `a` inside `O`. The state is **one-piece** at
`{12|34}` if `defect = 0` and at most one piece has `g_O > 0`.

## Theorem

For every state `t` preserving a block system, with quotient `tbar`, at the cusp `{12|34}`:

* **(A)** If a quotient piece `Obar` has genus `gbar >= 1`, then every piece `O` over it has
  `g_O >= gbar`. Consequently `#{positive-genus pieces of tbar} <= #{positive-genus pieces of t}`.
* **(B)** The natural graph map `phi: Gamma(t) -> Gamma(tbar)` is onto on `H1( ; Q)`.
  In particular `defect(tbar) <= defect(t)`.

**Corollary.** If `t` is one-piece at `{12|34}`, so is `tbar`. The same holds at the
other two cusps of the state (apply the theorem to the rotated tuple and to `sigma2(t)`).

Braid moves preserve the group `G`, hence the block system, and they commute with taking
quotients. So the quotient of a braid orbit that is one-piece at every state is again one-piece
at every state.

## Proof of (A): counting cycles

Let `O` be an `<x,y>`-orbit and let `Obar` be the set of blocks meeting `O`. Then `Obar`
is an `<xbar,ybar>`-orbit. The map `O -> Obar` (point to its block) is equivariant and
onto, and `<x,y>` permutes its fibres transitively. So all fibres have the same size
`k = |O| / |Obar|`.

Let `a` be one of `x, y, p` and take a cycle of `abar` of length `l` in `Obar`. Its preimage in `O`
has `kl` points. Every `a`-cycle there maps onto this `abar`-cycle, so its length is a
multiple of `l`. Hence there are at most `k` such `a`-cycles, and summing over cycles,

    c_O(a) <= k c_Obar(abar).                                         (1)

Since `|O| = k |Obar|`, the genus formula and (1) give

    2 g_O - 2 = |O| - c_O(x) - c_O(y) - c_O(p) >= k (|Obar| - c(xbar) - c(ybar) - c(pbar))
              = k (2 gbar - 2).                                       (2)

If `gbar >= 1`, then `2 g_O - 2 >= k(2 gbar - 2) >= 2 gbar - 2`, so `g_O >= gbar >= 1`.

Distinct quotient pieces have disjoint sets of blocks, so the pieces of `t` above them are
disjoint, and every quotient piece has at least one piece above it. This injects the
positive-genus pieces of `tbar` into those of `t`. Side B is the same, with `(z, w, p^-1)`.
(Inequality (2) is Riemann-Hurwitz for the piece curves `X_O -> Z_Obar`, but the proof
above does not need that.) QED (A).

## Proof of (B): homology

Realize the tuples by branched covers `f: X -> P1` (degree `|Omega|`) and
`fbar: Z -> P1` (degree `|Omegabar|`), with branch values `b1..b4`. The block system
gives a factorization `f = fbar o pi`, where `pi: X -> Z` is a branched cover of degree `m`
between closed oriented surfaces. Let `N` be a closed annulus in `P1` that separates
`{b1, b2}` from `{b3, b4}`, contains no branch value, and has radial coordinate
`r: N -> [0,1]`.

1. **Pieces and annuli.**
   * The components of `f^-1(D12)` correspond to the `<x,y>`-orbits, where `D12` is the
     closed disc of `P1 - int N` containing `b1, b2`. Likewise the components of `f^-1(D34)`
     correspond to the `<z,w>`-orbits.
   * The cover is unramified over `N`, with monodromy `p` around it. So the components of
     `f^-1(N)` are annuli, one for each cycle of `p`. Each annulus has one boundary circle
     on a piece over `D12` and one on a piece over `D34`.
   * So `Gamma(t)` is exactly the dual graph of this decomposition of `X`. The same holds
     for `Z` and `Gamma(tbar)`.
   * The map `pi` sends each piece of `X` into a piece of `Z` and each annulus onto an
     annulus. In the permutation description, `O` goes to `Obar` and a `p`-cycle goes to
     the `pbar`-cycle of its blocks. These are exactly the vertex and edge maps of `phi`.
2. **Collapse maps.** Make `Gamma` a metric graph with unit edges oriented from side A to
   side B. Define `c: X -> Gamma` as follows:
   * a point of a piece goes to that piece's vertex;
   * a point `q` of the annulus of edge `e` goes to the point at distance `r(f(q))` along `e`.

   This is continuous. Define `cbar: Z -> Gamma(tbar)` in the same way. Because
   `r(fbar(pi q)) = r(f(q))`, the square commutes on the nose:

       cbar o pi = phi o c.
3. **`c` and `cbar` are onto on `H1`.** The pieces are connected. So any closed edge path in
   `Gamma` lifts to a closed loop in `X`: cross each annulus along an arc, and join
   consecutive arcs inside the piece between them. The image of that loop under `c` is
   homotopic to the edge path. So `c` is onto on `pi_1`, hence on `H1`. The same argument
   applies to `cbar`.
4. **`pi` is onto on `H1( ; Q)`.** `pi` is a map of degree `m != 0` between closed oriented
   surfaces. Every class in `H1(Z; Q)` has the form `[Z] cap alpha` by Poincaré duality.
   By the projection formula,
   `pi_*([X] cap pi^* alpha) = pi_*[X] cap alpha = m ([Z] cap alpha)`.
5. **Conclusion.** `phi_* c_* = cbar_* pi_*`, and the right side is onto `H1(Gamma(tbar); Q)`.
   So `phi_*` is onto, and `b1(Gamma(tbar)) <= b1(Gamma(t))`. QED (B).

## Remarks

* Neither part assumes one-piece. Both are monotonicity statements: going down to a block
  quotient never increases the defect or the number of positive-genus pieces. Every
  positive-genus quotient piece lies under positive-genus pieces only.
* The converse fails. Lifts of one-piece orbits are often not one-piece, as in the
  lift-search data, e.g. the degree-10 genus-2 bases at m = 2.
* **Equality data for the lifting problem.** Suppose a lift of a one-piece cusp is itself
  one-piece, and the quotient has a positive-genus piece `Obar`. By (A), the unique
  positive-genus piece of the lift lies over `Obar`, and it is the ONLY piece over `Obar`.
  Every piece over a genus-0 quotient piece then has genus 0. In (2), equality
  `2g_O - 2 = k(2gbar - 2)` holds exactly when (1) is an equality for `x`, `y` and `p`,
  i.e. the returns (relative ramification) inside `O` are trivial. These are the
  constraints a lifting rule must describe.
* The earlier "quotient descent" (an earlier note) is a different
  statement: additivity of a pair, i.e. combinatorial smoothness.

## Checks (`descent_check.py`)

Each tuple is checked for:

| check | statement |
|---|---|
| T0 | equal fibre sizes |
| T1 | inequality (1) for `x`, `y`, `p` |
| T2 | inequality (2) for every piece |
| T3 | part (A) count |
| T4 | defect inequality |
| T5 | exact rank of `H1(Gamma; Q) -> H1(Gammabar; Q)` equals `b1(Gammabar)`, over `Q` with a spanning-tree cycle basis |
| T6 | one-piece descends |

These run exhaustively over all transitive product-one tuples in `S2 wr S3` and
`S3 wr S2`, and randomly over `S_m wr S_n` for `(m,n)` in
`(2,4),(2,5),(3,3),(4,2),(3,4),(2,6),(5,2),(4,3)`. Results are in `out/descent_check.txt`.

Results (Nov 7): ALL PASSED.

| population | transitive tuples | one-piece | quotient defect > 0 (T5 nontrivial) | quotient has a positive-genus piece (A nontrivial) |
|---|---|---|---|---|
| S2 wr S3, exhaustive | 98,552 | 28,224 | 32,512 | 11,176 |
| S3 wr S2, exhaustive | 322,644 | 77,616 | 46,092 | 0 (degree-2 quotient pieces have genus 0) |
| random, 8 shapes | 30,000 | 1,079 | 12,625 | 8,021 |
