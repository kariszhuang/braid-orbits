# Smoothness and defect

Elementary derivation; human review pending. Assumptions: the four permutations
have product one and together act transitively. Fixed points count as cycles.

## Smooth implies zero defect, with a stronger pair conclusion

Put p=ab=(cd)^-1, C=cycles(p), L=orbits(<a,b>), R=orbits(<c,d>).
For the disconnected pair cover, the sum of its component genera is

    gamma_L=(ind(a)+ind(b)+ind(p))/2-n+L >= 0.

Also L<=C because each p-cycle lies within a pair orbit. If ind(a)+ind(b)=ind(p),
then gamma_L=L-C. Consequently L=C and gamma_L=0.

The bipartite incidence graph with pair-orbits as vertices and p-cycles as edges
is connected by transitivity of the whole tuple. Hence C>=L+R-1. Substituting
L=C forces R=1. Therefore defect=C-L-R+1=0. The same argument with the two
pairs exchanged handles the other smooth equality.

Thus the opposite pair to a geodesic pair is transitive. Smoothness anywhere
in a braid component supplies a transitive pair somewhere in that component:
the complement can be moved to the front by braid block exchange (equivalently
the cyclic shift, whose conjugation terms preserve that pair's transitivity).
In particular no certified weak-failure component can contain a smooth state.

There is no counterexample to smooth-implies-zero-defect under these assumptions.
Dropping transitivity is different: on three points (12),(12),id,id has
ind(ab)=ind(id)+ind(id)=0 but defect=-1. Identity/triple strata are separate.

## Passport and construction constraints

Every product index is at most n-1. Since braid moves realize all orderings of
entry conjugacy classes, every partition of the four indices into two pairs
must have at least one pair sum <=n-1 in an everywhere-smooth component.
This is a necessary condition, not a certificate.

In the proposed family the indices are 2,2m,2m,2m, n=3m+f. Its genus is
1-n+(2+6m)/2=2-f. The larger pair sum 4m exceeds n-1, so smoothness requires
the side containing the single 3-cycle to be geodesic. Multiplication by its
two-transposition factorization has additive index precisely when the three
points are in distinct b-cycles: it joins those three cycles. Equal-length
cycles can be relabeled and permuted, giving the two declared normal forms.
Thus exact third-entry enumeration covers all possible all-smooth components
in this passport modulo simultaneous S_n conjugation.

This does not imply all zero-defect components are smooth; zero defect alone
allows both pair genera to be positive or planar pairs with several boundary
cycles in one component. The primary new question is additional everywhere-
smooth, identity-free transitive components, with full A_d/S_d generation
separately marked. The supplied degree-nine tuple generates order 324, not A9.

## An explicit infinite imprimitive genus-two family

For k>=3 let the points be three blocks B0,B1,B2, each with coordinates 1..k.
Let x=(1 2 ... k), y=(1 2), z=(yx)^-1 and tau=x^-1 y x y^-1. Tau is a single
3-cycle. Define a to act as tau on B0 and fix the other blocks. Define b to
send (j,i) to (j+1,i), with block indices modulo three. Define c to send

    (0,i) -> (1,x(i)), (1,i) -> (2,y(i)), (2,i) -> (0,z(i)),

and d=(abc)^-1. The three block maps for c multiply to zyx=1, so c^3=1.
For abc, the return map on B0 is tau*y*z*x=1 by the commutator definition;
therefore (abc)^3=1. Since both c and abc cyclically permute the blocks, they
have no fixed points and all cycles have length three. The passport is

    (3 1^(3k-3)), (3^k), (3^k), (3^k),

and its covering genus is always 2. All four entries are nonidentity and even.

Every braid translate still has one single 3-cycle supported within one block,
and three entries of type 3^k cyclically permuting the blocks. This follows
because every entry preserves the block system and the induced quotient C3 is
abelian, so Hurwitz moves only swap the induced block permutations. The three
points of the small cycle belong to distinct cycles of each large entry: every
large cycle meets a given block just once. Thus their product joins three
cycles, increasing the index by exactly two. Whichever half contains the
small entry is geodesic. Product one makes the two half-products inverse, so
the professor's smooth equality holds at EVERY braid translate. Defect is zero
throughout by the preceding proposition.

The generated group is (S_k wreath C3) intersect A_(3k), of order 3(k!)^3/2.
Indeed e=b^-1*c is the block-preserving element (x,y,z). The conjugates e^j*a*e^-j
are localized consecutive 3-cycles on B0 and generate A_k there; conjugating
by b gives independent A_k actions on all three blocks. The sign vector of e
has even total parity and is nonzero because y is odd. Its cyclic rotations
generate the entire two-dimensional even-parity sign space. Hence the base
group is exactly the even-total-sign subgroup of S_k^3, and b supplies C3.
This is transitive but imprimitive, never full A_(3k) or S_(3k).

Thus the professor's existence question has an infinite affirmative family
if proper transitive generated groups are allowed. This is a mathematical
construction, NOT a claim of literature novelty or a new weak counterexample.
Finite profiles measure quotient sizes separately and may be capped.

The geodesic/noncrossing interpretation is standard background, not a novelty
claim: Philippe Biane, Some properties of crossings and partitions (1997),
https://doi.org/10.1016/S0012-365X(96)00139-2 . The derivation above is explicit
and does not rely on an experimental sample.
