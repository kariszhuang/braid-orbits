# Exhaustive fixed-quotient enumeration for source genus 3 (October 30, 2026)

Engine: `src/main.rs` (Rust, crate `quotient_exhaust`). Campaign: `campaign.sh`, status: `status.sh`.
Live output: `data/quotient_exhaust_20261030/<stage>/` (done.jsonl per job, closed.jsonl,
capped.jsonl, POSITIVE_CANDIDATES.jsonl). Frozen binaries and SHA256 in `.../bin`, `.../bin/v1`, `.../bin/v2`.
Stop: `touch data/quotient_exhaust_20261030/STOP` (current chunk finishes; resumable).

## Population (exact)

Four-point tuples (s, t2, t3, t4) of degree d = m k, source genus 3, transitive, with a block system of m blocks
of size k on which s acts trivially and whose quotient (x, y, z) is a dessin of genus h (h = 1: fixed elliptic
quotient E; h = 0: fixed rational quotient). Up to simultaneous S_d conjugacy and up to moving s to position 1.
Every such tuple is (s, t2, t3, t4) = ((id;sigma), (x;f), (y;g), inverse) with (sigma, f, g) a flat
S_k-connection on the dessin's 2-complex. Total relative index R = 4 - k(2h - 2).

Dessins: one per B3-orbit (moves not involving s act on (x, y, z) alone and reach the whole quotient orbit),
x, y, z all nontrivial. (For h = 0 this omits quotients branched at only two points; see limitations.)

## Enumeration (exhaustive, gauge-fixed)

* Primal spanning tree edges are set to the identity. Dual-tree edges are SOLVED from prescribed non-root face
  holonomies. The 2h handle edges are free.
* Every class placement (sigma's on blocks, holonomy classes on non-root faces) with total index <= R, root
  remainder checked.
* v2: stabilizer-chain DFS. Each free element is taken only up to conjugation by the stabilizer of the earlier
  choices, so tuples differing by a global conjugation are not repeated.
* v2: solved second handle. The root holonomy word has the form U H2^{+-1} W H2^{-+1} Y (checked symbolically on
  162 genus-1 dessins; the code falls back to brute force if not). So H2 ranges over one centralizer coset per admissible root element.
* v3: S_k multiplication table built by Lehmer rank in parallel (k = 8: 3.2 GB, 14 s).

Per tuple: random-walk smoothness filter (300 literal moves; a non-smooth state is a rigorous negative), then
exact six-move closure modulo S_d and an exact isotriviality test on every closed orbit. Caps are unresolved,
never negative.

## Validation

* `--selftest`: the isotriviality test reproduces catalogue controls (a degree-9 genus-2 isotrivial orbit with
  12 states, and non-isotrivial genus-6 and genus-4 orbits with 240 and 192 states). Face words versus actual
  lifted cycle structure:
  270 random connections, 0 mismatches (after fixing the bug below).
* `--crosscheck`: full brute force over all (sigma, f, g) in S_k^{3m}: identical canonical tuple sets in 11
  cases (h = 0, 1; m = 3, 4; k = 2, 3).
* `--crosscheck2`: solved-H2 versus brute-H2: identical sets, zero root mismatches (k = 3, 4, 5; m = 3, 4, 5).
* Positive controls: the isotrivial translation class (one moving point, no fixed branching) closes as
  isotrivial (k = 5, 6, 7, 8); the 9 catalogued isotrivial genus-3 components with k = 3, m = 6 match
  the 9 closures found.

## Bug found in an earlier sampler (not included here, never launched)

`eval_word(..., usize::MAX)` returns (left, right) with the full holonomy in `right`, but the code used `left`,
which is always the identity. That sampler therefore never produced tuples with fixed branching on its root
face. Its tuples were valid (the genus check uses the real tuple), but its class labels and coverage were
biased. The frozen crate is left unchanged. Do not use it for coverage claims.

## Mathematics used to choose the stages

1. **Single moving block is m-independent.** If exactly one block j has sigma_j != id, the family is a family of
   covers C -> E with ONE moving branch point a and r <= 3 fixed branch points. The braid orbit is the
   point-pushing orbit of a in E minus {P_i}; collisions with unbranched points are harmless. All-smoothness
   depends only on the topology of (E, P_1..P_r), so the m = 3 dessin (j = 0 curve; three face points)
   realizes every single-moving-block configuration for every k. Non-isotriviality with r >= 1 holds by
   de Franchis (finitely many maps C0 -> E up to translation).
2. **Paramodular Theorem A** (a separate note): minimal genus-1 block size
   k <= 4 forces isotriviality in every degree, and k = 5 does too unless the family has RM by Q(sqrt 5). So the
   k = 3, 4 lanes are consistency checks. Search priority is k >= 6.

## Results (exhaustive within stated scope; every closed orbit is isotrivial)

| stage | h | k | m | built tuples | walk-pass | closed (all isotrivial) | status |
|---|---|---|---|---|---|---|---|
| h1k3a | 1 | 3 | 3-7 | 162,641,671 | 119,075 | 9 | complete |
| h1k4a | 1 | 4 | 3-6 | 232,254,491 | 1,330,585 | 540 | complete |
| h1k5a | 1 | 5 | 3-4 | 5,742,392 | 13,417 | 244 | complete |
| h1k6a | 1 | 6 | 3 | 4,173,938 | 8,880 | 351 | complete |
| h1k7a | 1 | 7 | 3 | 18,133,789 | 137,280 | 1,737 | complete |
| h1k8a | 1 | 8 | 3 | 72,916,670 | 1,172,281 | 6,984 | complete |
| h1k6b | 1 | 6 | 4 | 27,667,652 | 318,558 | 2,084 | complete |
| h1k5b | 1 | 5 | 5 | 50,941,960 | 580,601 | 1,450 | complete |
| h1k7b | 1 | 7 | 4 | 118,428,568 | 3,351,165 | 10,400 | complete |
| h1k6c | 1 | 6 | 5 | 278,197,720 | 8,807,565 | 15,810 | complete |
| h0k2 | 0 | 2 | 3-8 | 56,066,717 | 0 | 0 | complete |
| h0k3a | 0 | 3 | 3-5 | 151,490,936 | 0 | 0 | complete |

(v1 counts for h1k3a, h0k2, h0k3a, part of h1k4a; v2/v3 counts are smaller per class because global-conjugation
duplicates are no longer generated. Counts are tuples, not distinct braid orbits. Use `status.sh` for live totals;
job lines are deduplicated by job key.)

Consequence: in the genus-1 lane, the single-moving-block case is closed for k <= 8 in ALL degrees (by
reduction 1; k = 8 completed October 30 afternoon), and with Theorem A for every k <= 4. Multi-moving-block
configurations are closed for the (k, m) pairs in the table. Stages h1k5c (k=5, m=6) and h1k3b are in progress;
see status.sh. Theorem C (the same note) shows that k >= 5 cannot be excluded by
Jacobian geometry alone, so these exhaustive lanes, not modular forms, are the tool for k = 5..8.

## Limitations

* h = 0 omits quotient dessins branched at two points (z = id); catalogued genus-0 quotients with k = 10, 12,
  m = 2 are therefore not reproduced. Genus-0 lanes explode for k >= 5 (10^11 to 10^18 candidates) and need
  another method.
* k <= 8 (multiplication table). k >= 9 needs a table-free implementation.
* Multi-moving-block configurations depend on m. Only the m values in the table are covered.
