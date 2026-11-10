# Multi-active-point lifts (Nov 5, Claude; proofs parent-checked only, runs tool-checked)

Code: src/smooth_multi_active.rs (bin smooth_multi_active). Runs: /Volumes/T7/Research/
multi_active_20261105/{ctrl1,run1,run2,run3,run4}; binaries and source in HASHES.txt.

## Construction

Base: a transitive triple (x1, x2, x3) on m points, x1 x2 x3 = 1, group H = <x1, x2>; as
a 4-tuple (1, x1, x2, x3) it is trivially all-smooth. Geometrically this is a Belyi map
E -> P^1 of genus h. The source has n = m k points (beta, s):
- b, c, e lift x1, x2, x3 unramified: every return permutation is trivial;
- a fixes every fibre and acts on fibre beta by sigma_beta, nontrivial exactly on the
  active set S.

Geometrically, C -> E is branched over the points of S, all lying in the moving fibre over
t. Source genus: g = 1 + k(h - 1) + sum over S of ind(sigma_beta)/2. |S| = 1 is the old
localized mechanism (Proposition 1 and its generalization).

**Solving.** Gauge b's sheet labels to be trivial. Then c's return condition (on its
cycles) and e's (on the cycles of y = x1 x2) give one product equation per vertex of the
bipartite graph of (x2, y) cycles, with one edge per base point. The labels on a spanning
tree are solved leaf-first. That leaves exactly 2h - 1 + cyc(x1) free sheet labels, plus
one consistency equation at the root. There is no rejection sampling.

## Smoothness theorem

Let Gamma be the graph on base points joining beta and gamma when they lie in one cycle of
some H-conjugate of some x_i.

**Theorem 1.** If S is independent in Gamma, the whole braid orbit is all-smooth.
Proof: conjugation preserves "unramified", and keeps a in the kernel with active set
moved by H. Let y' be any entry and a' the conjugate of a in the same state. Every cycle of
y' lies over one base cycle C and meets each fibre over C exactly once, and C contains at
most one active point. So the a'-cycles at that point join distinct y'-cycles, and
different base cycles give disjoint pieces: the incidence graph is a forest. Hence a' is
additive with every entry, in every state.
Verification: across every run, every closure with an independent S came back all-smooth
(0 theorem failures).

Consequence: if H is 2-homogeneous, Gamma is complete and |S| <= 1. New examples need bases
with incomplete Gamma.

## Reduction to the old mechanism

Let N be the relation "beta = gamma, or beta and gamma are not adjacent in Gamma".

**Theorem 2.** If N is an equivalence relation, every lift of Theorem 1 is a single-block
localized lift for the coarser block system given by the N-classes.
Proof: N is H-invariant, so its classes D form a block system of the base. A cycle of an
x_i that met a class twice would put two N-related points in one cycle, i.e. make them
adjacent in Gamma. So each x_i is unramified over the quotient, and x_i^l fixes every point
of a class on whose class-cycle (of length l) it lies. Then b, c, e are unramified
over the coarse blocks (preimages of the classes), using their trivial returns on base
cycles. S is pairwise N-related, so it lies in one class, and a is active in exactly one
coarse block. That is the generalized localized/unramified criterion.

**Conjecture M.** For every transitive product-one triple, N is an equivalence relation
(equivalently, Gamma is complete multipartite).
- True for EVERY transitive product-one identity-free triple of degree m <= 10, by
  exhaustive enumeration over all of S_m. Gamma is incomplete for 0, 0, 0, 1, 0, 15, 2 and
  123 conjugacy classes of triples at m = 3..10, and N is an equivalence for all 141.
- None of those 141 bases is primitive. For a primitive base, Conjecture M says exactly
  that Gamma is complete.
- Proof for S_r acting on 2-subsets: to avoid the "disjoint" orbital, every x_i must be a
  transposition or a 3-cycle, which is too little index for a transitive triple. To avoid
  the "meeting" orbital, all three x_i must be fixed-point-free involutions, so x1 and x2
  commute and H is abelian. So Gamma is complete there.

If Conjecture M holds, Theorem 2 shows that multi-active unramified lifts never give new
structure.

## Runs (all closures by literal six-move closure, sm_scan)

| run | bases | k | active sets | lifts checked | closed all-smooth | class |
|---|---|---|---|---|---|---|
| ctrl1 | all, m = 3..6 | 2, 3 | single point | 1,861,114 | 330 | all localized_old; m=3,k=3 gives the known degree-9 component (12 states) |
| run1 | Gamma incomplete, m <= 8 | 2 | independent, >= 2 | 248 | 11 | all localized_old |
| run2 | Gamma incomplete, m = 9, 10 | 2 | independent, >= 2 | 6,160 | 98 | all localized_old |
| run4 | Gamma incomplete, m = 7, 8 | 3 | independent, >= 2 | 1,026,000 | 77 | all localized_old |
| run3 | all bases, m <= 7 | 2, 3 | ANY set of size >= 2 | about 1.87e10 | 19 (all at m = 6) | all localized_old; every other closure non-smooth (m = 7: 628,761 closures, none all-smooth) |

Run3 matters most: with arbitrary (dependent) active sets, every lift that stayed
all-smooth reduces to the old mechanism. So the independence condition of Theorem 1 is not
just sufficient: nothing outside it survived in this range. The active sets skipped by
the assignment cap (146 at m = 6, 5,471 at m = 7, both k = 3) are unresolved, not negative.
In run3 the v1 binary labels failed closures NEW_outside_localized; v2 labels them
not_all_smooth. Only closed_smooth rows are components.

## Conclusion

Idea 1 (several active points, unramified elsewhere) is fully implemented and validated.
Every example it produces up to these degrees reduces to the old localized mechanism, and
Theorem 2 plus Conjecture M explain why. New structure therefore needs ramification in
the other three entries as well, or a base on which Conjecture M fails (necessarily of
degree >= 11).
