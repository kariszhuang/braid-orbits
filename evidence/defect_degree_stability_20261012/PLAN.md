# Degree/genus stability pilot — pre-acquisition plan

Question: does the genus-coin distribution approximate defect at d>=60,
and do fixed-genus distributions stabilize near d>=2g?

Two distinct estimands: uniform labeled identity-free transitive product-one
tuples conditional on genus (not uniform braid components), and states sampled
inside specified braid components (passport fixed). Full A_d/S_d generation
is not certified. No finite pilot verifies “all g.”

1. Independent ambient S/A starts, degrees30,45,59,60,61,80,100,120;
   three seeded replicates of80000 attempts/cell. Retain ALL genus strata,
   including sparse ones; display cells with >=1000 observations. This is
   exact iid conditional sampling within each observed genus stratum.
2. Conditional S-permutation MCMC, genera1,2,4,8,16,24,32; degree grid
   max(g+2,2g-2),2g,2g+2,60,80,120, deleting invalid/duplicate cells.
   Four separately seeded starts/cell,50000 burn-in proposals,8 blocks of
   1250 observations, stride12. Propose75% a right transposition in one of
   the first three entries (solve fourth);25% one of six braid moves.
   Reject wrong genus, identities or intransitivity. Symmetric proposals
   make uniform conditional labeled tuples stationary ON EACH connected
   proposal class. Large-degree irreducibility/mixing are NOT proved.
   Starts are random relabelings of (a,a^-1,c,c^-1), a full cycle and
   ind(c)=g, with randomized forest cycle types and braid repositioning.
3. Adversarial braid controls: four-full-cycle tuples at60,61,80 (fixed
   defect parity) and the established genus2 three-block family at60,90,120
   (identically zero defect by the prior theorem). Two chains/control,
   same block schedule; these are proper-group or special-passport controls,
   NOT evidence about typical full symmetric groups.

Record joint (C(g1g2),l,r), δ=C-l-r+1, genus, parity and chain/block IDs;
save each literal initial/final tuple for independent replay. Compare TV to
the original coins; also compare fixed-genus degree cells and replicate/half
discrepancies. Predeclare TV<=0.05 as a descriptive close-fit threshold and
cross-degree TV<=0.05 as tentative stability; neither is an equivalence proof
when MCMC diagnostics fail. No iid p-values for MCMC or braid observations.

One auxiliary CPU, no GPU, <=10-minute subprocess budget, RAM histograms,
buffered output, frozen code/binary on mounted T7 with30GiB guard. The existing
4+2 smooth campaign and its frozen files/lock/controller stay unchanged.
The isolated128-point statistic kernel does not change production MAX_D=34;
validate it against literal Python permutations and finite small-degree tests.
