# All-zero conjecture: counterexamples found (November 8, 2026)

Conjecture under test: a positive-genus braid orbit whose group is primitive and
whose defect is 0 at every state must be a C_p orbit. (Consequences: every positive-genus
one-piece-everywhere orbit is imprimitive, and isotrivial => imprimitive.)

## Result

FALSE in genus 1. Seven distinct braid orbits (COUNTEREXAMPLES.jsonl), all genus 1, identity-free,
primitive, non-abelian, defect 0 AND one-piece at every state, isotrivial (as genus 1 +
one-piece everywhere forces):

| degree | action | states | entry cycle types | generated group |
|---|---|---|---|---|
| 21 | PSL(3,4) on the 21 points of PG(2,4) | 28 (two orbits) | 2^8 1^5, 2^8 1^5, 2^8 1^5, 7^3 | PSL(3,4), order 20160 |
| 25 | product action on 5 x 5 | 12 (three orbits) | e.g. 2^10, 2^10, 2^10, 5^5 | order 7200, primitive |
| 25 | product action on 5 x 5 | 60 (two orbits) | e.g. 2^8, 10^2 5, 2^10, 2^10 | order 14400, primitive |

Example (PSL(3,4), 1-based, rightmost-first, g1 g2 g3 g4 = 1):
(1 11)(2 20)(3 9)(5 14)(7 21)(8 16)(10 13)(15 18); (1 20)(3 8)(4 12)(5 16)(6 17)(9 10)(13 14)(18 21);
(2 21)(3 18)(4 20)(5 19)(10 16)(11 14)(12 15)(13 17); (1 14 10 11 4 15 2)(3 16 18 12 20 21 7)(5 17 6 13 9 8 19)

No primitive non-abelian all-zero orbit of genus >= 2 was found in any group searched.

## How

* make_groups.py builds 127 group entries (PSL/PGL/PGammaL(2,q) q<=32, AGL1 subgroups and
  AGammaL(1,q), AGL/ASL(n,p), PSL(3,q) q<=4, S_k/A_k on 2- and 3-subsets, S5 wr S2 and
  (A5xA5).2 product action, M11, M12, controls C7, D5). Searching G covers every subgroup.
* src/bin/allzero.rs: every identity-free product-one positive-genus 4-tuple in G up to
  conjugacy (g1 a class rep of minimal index, g2 modulo C(g1)), defect-0 filter at the three
  cusps of the state, six-move closure with early exit. 111 groups complete; skipped:
  imprimitive entries and orders > 120000 (AGL(4,2), AGL(3,3), ASL(3,3), PGammaL/PSigmaL(2,32)).
  Run log ALLZERO_RUN.jsonl. Controls: C7 gives 2 all-zero orbits (cyclic), D5 none beyond C5.
* verify_cex.py (independent Python + SymPy): product 1, genus, primitivity, group order,
  independent closure with its own canonical form, defect and one-piece at every state: all
  10 reported records CONFIRMED, 7 distinct orbits (least canonical state).
* Many IMPRIMITIVE all-zero orbits also occur (cyclic C_{q+1} in PGL(2,q); order 500/1029/2058
  subgroups of AGL(2,5), AGL(2,7)); by the descent theorem their quotients are all-zero.

## Consequences

* The all-zero conjecture fails in genus 1; a corrected statement must assume genus >= 2.
* "Positive-genus one-piece everywhere => imprimitive" fails in genus 1 (these are one-piece).
* "Isotrivial => imprimitive" fails in genus 1 (genus 1 + one-piece => isotrivial).
  Genus >= 2 versions remain open.

## Nov 8 (later): speedups and the genus >= 2 attack

* allzero_v3: probe (64 random moves, reject on positive defect; sound) BEFORE the canonical
  form, NEON kernels for permutations of up to 64 points (composition = four `tbl4`), and
  orbits<g3,g4> computed as orbits<g3, g1 g2> (no inverse). On 10 validation groups the tuple
  counts, defect-0 counts and primitive all-zero counts are identical to the first run;
  AGL(2,7) 300 s -> 53 s, M12 122 s -> 38 s, PSL(3,4) 6.4 s -> 1.4 s, S5 wr S2 19.4 s -> 4.6 s.
  (The count of IMPRIMITIVE all-zero records is not de-duplicated across threads and may vary.)
* onepiece_v8 (census): --all-zero mode (state property = defect 0; genus enforced by keying the
  half-cusp tables on the side genus), probe before the canonical form. Checks: genus-2 all-zero
  d4..8 = 0,0,0,0,2 (as the complete genus-2 censuses), genus-1 d8 = the 9 one-piece orbits,
  C5/C7 controls; one-piece d9 genus 2 932 s -> 164 s; all-zero d9 genus 2 = 1 imprimitive orbit.
* Genus >= 2 attack running at data/one_piece_20261105/allzero_g2 (run.sh):
  groups2.txt (43 groups: PSL/PGL(2,q) q<=61, PSL(4,2), PSL(3,5), PGL/PGammaL(3,4),
  S9/A9 on 2-sets, product actions of S4/S5/A5/L2(7), diagonal A5xA5, M22?, and the previously
  skipped AGL(4,2), AGL(3,3), ASL(3,3), PGammaL/PSigmaL(2,32)) with --min-genus 2, plus the
  degree-10 all-zero census for genus 2..9 (covers every primitive group of degree 10).

## Nov 8: all-zero orbits that are NOT of C_p type (genus >= 2)

classify_allzero.py classifies every genus >= 2 all-zero record by its generated group and block
quotients. Of the records from the first run (genus >= 1 search): 2350 cyclic, 162 towers over a
cyclic positive-genus quotient, 768 other solvable (all composition factors of prime order, e.g.
order 2058 = 7^3.6 inside AGL(2,7), order 500 inside AGL(2,5)), and 273 records with a
NON-SOLVABLE group.

verify_nonsolvable.py (independent closure, verify_cex.py) reduces these to 245 DISTINCT orbits,
all confirmed defect 0 at every state, all imprimitive: degree 25, genus 2..20, generated groups of
order 300, 600, 1200, 3600, 7200 containing A5; two block systems (rows and columns of the 5x5
grid). Their quotients are degree-5 covers, often with identity entries (3- or 2-point covers):
S5/A5 pieces of genus 0, 1 or 2, and D5/F20/C5 pieces (C5 with two identity entries = the Kummer
cover z -> z^5 branched at two points), i.e. fibre products of two degree-5 covers.
243 of the 245 are non-isotrivial and NOT one-piece everywhere (separating splits only, like
C_p); 2 (genus 2) are one-piece everywhere and isotrivial.
Smallest non-isotrivial one: genus 2, 24 states, group order 7200:
(3 5)(8 10)(13 15)(18 20)(23 25); (1 18 2 16 3 17)(4 20)(5 19)(6 13 7 11 8 12)(9 15)(10 14)(21 23 22)(24 25);
(1 17 8 14 25)(2 18 9 15 21)(3 19 10 11 22)(4 20 6 12 23)(5 16 7 13 24); (11 16 21)(12 17 22)(13 18 23)(14 19 24)(15 20 25)
So "every all-zero orbit is built from C_p pieces" is false even in genus >= 2 (imprimitive
non-solvable examples); the open statement is only the PRIMITIVE genus >= 2 one.
Files: classify_out.txt, NONSOLVABLE_ALLZERO.jsonl.

## Nov 8 (evening): GENUS-3 PRIMITIVE COUNTEREXAMPLE

The genus >= 2 attack (allzero_v3, groups2.txt, --min-genus 2) found 2 orbits in PSL(3,5) acting
on the 31 points of PG(2,5): genus 3, 4 states each, tuple generates all of PSL(3,5) (order
372000), primitive, non-abelian, defect 0 AND one-piece at every state; isotrivial
(a fixed genus-3 curve). Passport 2^12 1^7, 4^6 1^7, 4^6 1^7, 4^6 1^7.
Independently confirmed by verify_cex.py (Python + SymPy): GENUS3_COUNTEREXAMPLES.jsonl.
(2 15)(3 13)(4 16)(5 14)(6 12)(7 19)(8 21)(9 18)(10 20)(11 17)(22 23)(24 26);
(1 3 5 2)(7 24 28 16)(8 12 25 21)(9 30 17 26)(10 18 14 31)(15 19 27 23);
(2 31 11 16)(3 27 8 14)(4 28 10 12)(5 29 7 15)(6 30 9 13)(22 24 25 23);
(1 16 12 13)(2 25 10 20)(3 9 28 22)(5 18 26 29)(6 27 19 11)(8 24 17 31)
Consequences: "genus >= 2, primitive, all-zero => C_p" is FALSE; so are "positive-genus
one-piece everywhere => imprimitive" and "isotrivial => imprimitive" in genus 3.

Census v9: sigma3-orbit normalization of the half-cusp tables (one g3 per sigma3-orbit; sigma3
keeps g1, g2, p and the cluster). Identical results to v8 on d9 g3, d9 g2, d8 g1, d9 g2 all-zero,
d8 g2 all-zero, d7 genus-2 need-good; 8-16x faster (d9 g2 all-zero 343 s -> 25 s).

## Nov 8: work split across two machines

* macbook (MacBook Air M3, 4P+4E cores, 16 GB; `ssh macbook`): groups 33-43 of groups2.txt
  (groups2_remote.txt: L2(7) wr S2 on 64 and 49, the three A5xA5 diagonal groups, M22?,
  PGammaL/PSigmaL(2,32), AGL(4,2), AGL(3,3), ASL(3,3)) with the same allzero_v3 binary (hash
  e52356ba...), nice 5, 8 threads, caffeinate; directory ~/onepiece_remote_20261108. Results are
  copied back to T7 allzero_g2/remote_macbook.jsonl.
* local M1 Max: groups 25-32; stop_after_32.sh stops the local run after group 32.
* Benchmark (PSL(3,4) + M11, identical tuple/defect counts on both machines): macbook 2.45 s
  (8 threads, while its own load was ~120); local 2.0 s idle earlier, 3.8 s while busy. So the
  macbook is worth roughly 80% of the local machine on this code, i.e. about +80% throughput.

## Nov 8 (night): macbook share finished; PSL(3,7) started there

* macbook finished all 11 groups of groups2_remote.txt (files remote_macbook_part1.jsonl and
  remote_macbook.jsonl on T7): no primitive non-abelian all-zero orbit of genus >= 2 in
  L2(7) wr S2 (64, 49 points), the A5xA5 diagonal groups (60), M22 (order 443520 confirmed),
  PGammaL(2,32) = PSigmaL(2,32), AGL(4,2), AGL(3,3), ASL(3,3) (the last is inside AGL(3,3)).
  Imprimitive all-zero records: 725, 24/16/30, 62, 312, 267, 70 (to be classified by structure).
* Local degree-10 all-zero census (v9): genus 2: 8 orbits, genus 3: 0, genus 4: 8; none primitive.
* PSL(3,7) on the 57 points of PG(2,7) (groups3.txt, order 1876896) now running on the macbook,
  every genus >= 1, caffeinate -dims, output psl37.jsonl (synced to T7). Purpose: does the
  projective-plane family PSL(3,q) keep giving primitive one-piece orbits (q = 4: genus 1,
  q = 5: genus 3)?

## Nov 8 (night): PSL(3,7) also gives primitive one-piece orbits

Run on a second machine (groups3.txt): 2 orbits in PSL(3,7) on the 57 points of PG(2,7), genus 1, 4 states,
generating all of PSL(3,7) (order 1876896), primitive, defect 0 and one-piece at every state,
isotrivial, not all-smooth. Entries: three homologies of order 2 (2^24 1^9: a line of 8 points
plus a centre fixed) and one elation of order 7 (7^7 1^8: an axis of 8 points fixed).
Independently confirmed (verify_cex.py): PSL37_COUNTEREXAMPLES.jsonl.
Family so far: q = 4 (genus 1, 28 states), q = 5 (genus 3, 4 states), q = 7 (genus 1, 4 states),
all built from central collineations (homologies/elations) of PG(2,q). Conjectured infinite family.

## Nov 8: degree-10 all-zero census COMPLETE (onepiece_v9, every group of degree 10)

Genus 2..9: 8, 0, 8, 0, 12, 0, 7, 0 all-zero orbits (35 in total), NONE primitive
(T7 allzero_g2/az_d10_g*.jsonl). So no primitive all-zero orbit of genus >= 2 at degree 10,
including A10 and S10 (C_p does not occur at degree 10).

## Nov 9: correction and publication

The q = 4 collineation search printed 5 PSL(3,4) orbits; an independent re-closure keyed by the
catalogue's canonical state key (modulo simultaneous S_21 relabelling) shows only 2 distinct
orbits. Distinct primitive one-piece orbits: 29 = 24 projective-plane (q = 4: 2, 5: 4, 7: 6,
8: 2, 9: 2, 11: 4, 13: 4) + 5 on the 5x5 grid. All 29 are in public catalogue v18
(evidence/onepiece_publication_20261108).

## Nov 9 update

- PSL(3,7) whole-group search COMPLETE: every 4-tuple in the 57-point action,
  1,229,837,638,065 tuples, 22 classes. Exactly 2 all-zero orbits, both primitive genus 1,
  and they are the two collineation-search q=7 genus-1 orbits (same catalogue IDs). Zero
  cyclic, zero imprimitive-subgroup. The q=7 genus-4 collineation orbits generate PGL(3,7)
  (order 5,630,688), so they lie outside this search, consistently. Files:
  psl37_complete_20261109.jsonl, q7_group_check_20261109.jsonl.
- q=11 full collineation search (partial results at job 3500/6290; stopped at 4500, see below):
  9 hits, independently re-closed by verify.py, 9 DISTINCT IDs modulo S_133.
  Genus 6: 3 orbits, 4 states, one-piece everywhere, isotrivial; all three are ALREADY in
  catalogue v18 (IDs 72a5942f571f.., d1414a2a13e7.., e24285a1c95d..), so not new.
  Genus 30: 3 orbits, 12 states; genus 42: 3 orbits, 4 states; both defect 0 everywhere
  but NOT one-piece (non-isotrivial). Partial; not in the public catalogue.
  Files: q11_partial_hits_20261109.jsonl, q11_partial_verify_20261109.jsonl.

## Nov 9: all searches stopped

- q=11 full collineation search stopped at job 4500/6290 (about 15 h), 10 all-zero records; the
  10th (genus 30) came after the 9-hit check above and is not yet verified. Partial
  collin_q11_full.jsonl/.log kept on T7 allzero_g2; no resume file (restart = job 1 or a job range).
- MacBook: nothing running since PSL(3,7) finished Nov 9 00:00. Full account of the MacBook
