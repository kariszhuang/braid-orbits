# Exact reductions (human mathematical review pending)

Write the tuple as (a,b,c,d), with composition `(xy)(u)=x(y(u))`, product one,
p=ab, and indices A,B,T,V. The constructed pair has ind(p)=A+B.

## Safe index feasibility

Identity-free positive genus requires 1<=T,V<=n-1 and
`A+B+T+V >= 2n`, with even total. Everywhere smoothness requires a smooth
state after moving the entries into each of the three unordered pairings.
At least one side of each pairing must have index sum <=n-1. The constructed
pair already supplies that side for (a,b)|(c,d). The remaining conditions are
`A+T <= n-1 OR B+V <= n-1` and
`A+V <= n-1 OR B+T <= n-1`.
These are necessary, not sufficient. Precomputing allowed V for every T and
drawing only T with a nonempty set removes impossible positive-genus passports;
the generated V is still tested. This changes sampling weights, not the exact
decision on a fixed candidate. In particular T>=n+1-(A+B).

Equivalence to the old sorted-index passport test: sorted indices x0<=...<=x3
have largest minimum paired sum min(x0+x3,x1+x2), among the three matchings.
Thus all three matching minima <=n-1 iff x0+x3<=n-1 OR x1+x2<=n-1.

## Ten short words reduce exactly to five factors

Set F0=b^-1ab, F1=a, F2=b, F3=aba^-1, F4=p^-1bp.
The first-pair product indices after the old ten tested words equal ind(Fj c):

| Old word (move IDs) | j | Compared index sums |
|---|---|---|
|1|0|A+T, B+V|
|4|1|A+T, B+V|
|0,1|2|B+T, A+V|
|0,4|3|B+T, A+V|
|3,1|4|B+T, A+V|
|3,4|2|B+T, A+V|
|2,1|2|B+T, A+V|
|2,4|3|B+T, A+V|
|5,1|4|B+T, A+V|
|5,4|2|B+T, A+V|

IDs0,1,2 are positive adjacent Hurwitz generators and3,4,5 their inverses.
Substitution of d=(abc)^-1 and conjugacy/inversion invariance of index proves
the identities. Five tests in j order0..4 preserve even the first failing
word of the old screen. All product-one S4 triples and larger random controls
check the identities literally. Direct composition cycle traversal avoids
allocating/materializing each product.

## Transitivity reduction

Decompose a,b into A+B transpositions. Their edge graph has the same point
orbits as <a,b>. With r connected components, A+B>=n-r, so r>=n-(A+B).
Every p-cycle is contained in a pair orbit, so r<=#cycles(p)=n-(A+B).
Equality follows: p-cycles are exactly the <a,b> point orbits. Consequently
<a,b,c> is transitive iff <p,c> is transitive. This reduction is used only
for the constructed geodesic pair, not for arbitrary tuples.

## Early closure and proof cache

End Hurwitz moves leave the first-pair product unchanged and merely swap
indices within an end pair, so preserve smoothness. The two middle successors
can be checked raw, before any canonicalization. All six generators still
enter BFS. Indices are carried by adjacent swaps; defect and pair orbits need
not be recalculated during rejection BFS. Closed components are subsequently
saved with all literal statistics and independently verified with all six moves.

Rooted traversal canonicalization still minimizes over every root and uses
the identical three128-bit words. Only lexicographically losing suffix packing
is avoided. Existing state hashes/certificates remain compatible.

For a reached state with root path w and a nonsmooth root witness v,
w^-1 v is a nonsmooth witness from that state, up to simultaneous relabeling.
Each cached word is replay-checked literally before insertion AND use.
Bounded eviction merely repeats work. Capped/timeout results are never cached
as negative; they remain replayable leads. The negative cache is RAM-resident
within each job; no cross-job persistence is claimed.

## Recovery/I/O contract

One open lead handle per job replaces open/close per promoted candidate.
The lead buffer is still flushed before BFS. Resolutions flush periodically
and at completion; if interrupted between flushes, the durable seed remains
an unresolved lead rather than an exclusion. Neither old nor new writer uses
fsync; this contract covers process interruption, not guaranteed power-loss
durability. VERIFIED is written only after complete independent auditing.
