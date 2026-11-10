# C209: exact small-degree exclusions from full-cycle twists

The goal remains unmet: no new non-isotrivial genus-three all-smooth component.
The progress is mathematical coverage rather than an increased example count.
See [the proofs and precise populations](MATHEMATICS.md) and CHECKS.json.

- Source-genus-three degree11 is now completely negative for identity-free
  transitive quadruples: two elementary obstructions plus115390 independently
  reconstructed and literal-negative necessary S11 states.
- The ENTIRE degree13 triangle index budget(6,6,6,12) is excluded by a
  fixed-point rotation argument, not only the previously completed involution
  passport. Other degree13 budgets remain eligible.
- A balanced-support theorem forces uniform ordinary profiles when one entry
  is a single cycle and two others have the matching number of cycles.
- Its general support bounds remove all mixed ordinary profiles for degree15
  budget(4,10,10,10) with special types5,4 2,3^2,3 2^2. The last two pure
  cases plus the pure2^4 case now have44800 further literal negatives.

All raw evidence is retained on T7. These are not independent random samples,
new braid components, or a global nonexistence theorem. Geometric proof review
is separate; no Hodge/classifier exclusion was used in these finite negatives.

Reproduce ONLY with a coordinated spare CPU slot and a fresh output directory:

```sh
python3 scripts/run_genus_three_chords.py NEW-OUT PATH-TO-VALIDATED-AUDITOR
python3 scripts/explore_genus_three_saturated.py NEW-OUT COMPLETIONS.tsv AUDITOR
```

The second command accepts the old verified pure3^5 library; it rechecks every
entry, uniqueness and the exact character count before reuse. Both commands
freeze sources, retain pending local-test survivors, and run corruption tests.
Do not reschedule the completed cases as new discovery work.

Resource handoffs retain all ten stream cursors,1793 known identities and
literal words. Five broad discovery workers continued while the auxiliary
slot did exact checks; the final continuation returns to six total workers.
Check the live pointer/viewer, not a launch PID or the historical382 counter.
No publication, TeX change, deletion, storage migration or GPU rollout.
