# Exhaustive smoothness target, degrees 5--8

User-authorized after the original degree-nine restricted passport test.
All product-one tuples in S_n whose generated group acts transitively are
eligible: proper groups, odd entries, all genera, and arbitrary passports.
Identity entries are INCLUDED in the smooth-implies-zero test; they are
reported separately and excluded only from the nontrivial everywhere-smooth
orbit search. No weak-intransitivity filter is used.

Any smooth tuple has an index-additive first or complementary pair. Rotating
the latter pair to the front preserves product one, transitivity and defect,
and stays in its braid component modulo simultaneous conjugation. Enumerate
every index-additive ordered pair up to simultaneous S_n conjugation:
enumerate every product cycle partition, every a in S_n, and set b=a^-1 p.
Retain ind(a)+ind(b)=ind(p). No passport or genus is omitted. Canonicalize
the possibly disconnected pair by rooted equivariant traversal in each
component, sorting the component keys. Enumerate ALL n! third entries for
each pair representative; d=(abc)^-1 is forced. This covers every smooth
tuple and every possible everywhere-smooth component, not every nonsmooth
tuple. It is an exhaustive search of the requested smooth population, NOT
a census of all braid components.

For each transitive occurrence, independently count complementary pair
orbits and compute defect; save any implication counterexample literally.
Pair orbit masses use exact centralizer sizes (root automorphisms times
factorials of repeated components). Weighted left-geodesic counts and
both-geodesic counts give ALL labeled smooth transitive tuples by inclusion-
exclusion: 2*left - both. Counts of reduced occurrences are distinct from
labeled tuples and quotient states. Do not present incomplete weighted
counts as complete counts. Coverage proof awaits human mathematical review.

Each identity-free start receives a literal nonsmooth braid word, membership
in an exactly closed smooth certificate, or an explicitly unresolved status.
Short witnesses use the existing literal braid implementation. Closure and
S_n state keys reuse the independently checked smooth engine; caps/timeouts
are not negative results. Outcome logs reconstruct each seed from the saved
pair frame and lexicographic permutation rank. Flush after every frame.

Run only under the shared campaign lock, with maximum four single-CPU degree
jobs (within the latest seven-slot budget). GPU is not needed for these small
exact jobs. Freeze every Rust source, plan, verifier helpers and compile-time
fixture, hash source and binary. Build/regressions precede launch. Default
one-hour per-degree safety cap; unfinished frames/closures remain unresolved.
Raw output belongs on mounted T7, requiring at least 30 GiB free. Prior weak
campaign and completed restricted smooth mission remain untouched. No TeX edit.
