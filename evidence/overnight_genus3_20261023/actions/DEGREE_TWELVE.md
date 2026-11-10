# Complete negative result in source genus three, degree twelve

October23,2026. Product-one transitive quadruples, no identity entry; full
six-generator Hurwitz orbit, states modulo simultaneous S12 conjugation.
This is a finite-degree exclusion, NOT global nonexistence.

## Scalar reduction and transposition cases

For sorted positive indices a<=b<=c<=e, Riemann--Hurwitz gives sum28.
Each of the three pairings must have a side of sum at most11. Equivalently,
a+c<=11 and min(a+e,b+c)<=11. The only budgets are

    (1,7,10,10), (1,8,9,10), (1,9,9,9), (2,8,9,9).

C225's independently source-reviewed Gabai--Kazez argument excludes every
positive-genus full all-smooth tuple containing a transposition, in any
degree. Thus only(2,8,9,9) remains. Its index2 entry is either one3-cycle
or two disjoint transpositions. Braids can place the entries as(a,c,b,e)
of indices(2,9,8,9). Smoothness forces p=ac to be a full12-cycle, since
the complementary index sum17 exceeds11. Relabel p=(0 1 ... 11).

Write Q1(a,c)=(aca^-1,a). Then Q1^(2k) conjugates a,c by p^k.
Following by Q2^-1 puts a_k=p^k a p^-k and b first. Their product must
have index10; otherwise this is a literal nonsmooth endpoint. Applying
Q3^-1 BEFORE those moves instead puts a_k and e first; their product must
have index11. These necessary twists are used below, not asserted sufficient.

## C226: single3-cycle special entry

C209's saturated single-cycle theorem forces both index9 entries to have
cycle type4^3. Since c=a^-1p has three4-cycles, the three support points
of a split the cyclic order of p into three equal arcs. Up to a rotation
there is one marked frame:

    a=(0 4 8),   c=(0 1 2 3)(4 5 6 7)(8 9 10 11).

The orientation of a follows the cyclic order: the reverse orientation
does not split p into three cycles. Residual rotations by4 stabilize the
frame; we do NOT divide candidate counts freely by3 or call them state counts.

For every shift k, the support{k,k+4,k+8} must meet all three e-cycles.
Thus each4-cycle of e contains one element in each residue class modulo4.
Choose three bijections matching residues1,2,3 to the three residue0 points,
and orient each resulting4-cycle: (3!)^3*(3!)^3=46,656 labelled choices.
Set b=p^-1 e^-1; this is forced by product one, not sampled.

There are23,401 choices with ind(b)!=8. Each of the remaining23,255 has
a b-cycle containing two letters of one residue r modulo4. Consequently
the a_r support repeats a b-cycle and is not additive. The literal word
Q1^(2r),Q2^-1 is nonsmooth; r ranges0..3, so at most7 moves suffice.

The Python producer constructs the46,656 choices using residue bijections.
The independent Rust auditor instead enumerates all5,775 unordered partitions
into4-blocks and all216 cycle orientations, the complete4^3 conjugacy class
of size1,247,400; exactly216 support partitions pass the residue condition.
It reconstructs producer ordinals, checks every literal endpoint and obtains
the identical46656*12-byte enumeration SHA256. No cap/unresolved case occurs.

## C227: double-transposition special entry

C209's active-support theorem applies with index2, support4 and two active
cycles. Every cycle of c or e has length between ceil(12/4)=3 and
floor(2*12/4)=6. Their possible types are therefore6,3,3 or5,4,3 or4,4,4.

Enumerate all1,485 double transpositions a. Exactly990 have a*p with three
cycles; exactly180 also satisfy these cycle-length bounds. Simultaneous
conjugation preserving p is precisely a cyclic rotation, giving16 necessary
marked frames. This is a coverage reduction, not a count of braid states.

For each frame, partition the12 letters into three prospective e-cycle
supports, each size3..6. There are42,735 such unordered partitions. For each
k the two edges of a_k must induce a tree on these three supports: neither
edge can be internal, and they cannot connect the same pair of supports.
This condition is exactly the additivity/full-cycle condition for a_k*e,
independent of the orientation of e's cycles. Across16 frames377 support
partitions pass all12 shifts. Enumerating every cyclic orientation gives
104,232 labelled e choices. For each, b=p^-1 e^-1 is again forced.

-57,314 have the wrong b-index and cannot have this passport.
-46,830 of the remaining46,918 have a literal negative Q1^(2k),Q2^-1.
-The final88 ALL have a literal negative among Q2^2, Q1 Q2^-2, Q2^-2.

All88 late witnesses are at most3 moves; no BFS cap or unclosed candidate
is counted negative. The initial Rust output labels them UNRESOLVED because
that stage intentionally has no claim about subsequent moves. The separate
Python closure/witness check and independent full auditor resolve all88.

The independent Python auditor imports no producer functions: it enumerates
the1,485 edge-pair matchings by combinations of66 edges, constructs support
partitions by restricted-growth3-colour words rather than subset masks, and
orients cycles with itertools permutations. It reconstructs EVERY tuple,
replays EVERY negative and checks exact equality, including no duplicates,
of the full104,232-record (a,e) set against the Rust binary enumeration.
No unresolved cases remain. Counts are necessary marked labelled slices,
not distinct braid components; no full A12/S12 generation assumption is used.

## Conclusion and limitations

Combining both index2 cases with the global transposition obstruction excludes
ALL identity-free transitive source-genus3 all-smooth full braid orbits in
degree12, for every transitive generated permutation group. Together with
inherited complete9/10/11 exclusions and the general d>=2g+3 bound, any
target must have d>=13. Degree13 star cases remain open. No global upper
bound or global nonexistence theorem is claimed.

Reproduction sources and immutable result hashes are in D12_CHECKS.json.
The complete negative proof is the coverage argument plus literal moves;
random sampling, geometry classification of candidates and long walks are
not used. The C225 topological input is separately stated and primary-source
reviewed, not replaced by computational evidence.
