# Same verification problem, independent implementation

For each geodesic ordered pair (a,b), enumerate c in lexicographic S_n order
and set d=(abc)^-1. The native permutation table is independently generated
by next-lexicographic-permutation; indices are recomputed by cycle traversal.
Its rank is therefore exactly the Python itertools.permutations rank.

All ranks increment the enumeration counter. Under POSITIVE_ONLY, omit any
identity entry, then tuples with ind(ab)+ind(c)+ind(d)<=2n-2. The first pair
is checked index-additive, so this is precisely the nonpositive-genus test.
Keep both omission counts separately and match OMITTED.tsv.

Compute transitivity as connected components of permutation edges for a,b,c.
Since d belongs to their generated group, this is the four-entry group as well.
This union-find implementation is separate from Python's orbit traversal.
For every retained transitive tuple, explicitly recompute the two pair orbit
counts and check n-ind(ab)-r12-r34+1=0 and r34=1, as in the old auditor.
Count second-pair index additivity separately; the labeled smooth-population
mass uses symmetry/inclusion-exclusion 2*left-both, unchanged from the old audit.

Every eligible rank consumes exactly one little-endian rank/token record.
Dictionary tokens must be consecutive; move words must be nonempty and in
0..5. Positive/partial statuses are unsupported and rejected, not counted as
negative. For each word, perform literal Hurwitz moves on all four permutations:
(x,y)->(xyx^-1,x), inverse (x,y)->(y,y^-1xy). Then recompute smoothness using
the definition on the resulting tuple. The endpoint must not be smooth.
There is no witness sampling, representative cache, or inherited production
filter predicate. Full audits require exact stream exhaustion.

The Python controller independently regenerates all geodesic first-pair frame
classes through cycle-splitting paths for every product cycle partition.
It checks disconnected pair canonical keys, stabilizer/orbit masses, complete
frame population, every per-frame producer count, full completion markers,
all omission counts and aggregated population/labeled counts. Native outputs
cannot authorize a census conclusion unless all these checks finish.

Soundness is the same as the original audit: every possible smooth tuple is
represented by an index-additive first or last pair (swapped if necessary),
and each eligible tuple has a literal nonsmooth braid translate. First-pair
conjugacy reduction does not discard any choice of third entry because every
c in S_n is still enumerated. This is a negative-only verifier; it does not
assert that closed-component certification or arbitrary-degree enumeration
has been accelerated. Human review of mathematical arguments remains pending.
