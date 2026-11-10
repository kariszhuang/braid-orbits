# Proof-directed special smooth components

Authorized October12: find additional special everywhere-smooth cases, prove
reductions first, then implement and run efficient searches. No new weak target.

Population: odd cyclic block counts r=3,5,7; block sizes k=3..8 with rk<=25.
The r=3,k=3 population is omitted entirely (known single-3-cycle construction).
In every other r=3 job, identity and single-3-cycle local entries are separated
and excluded from promotion. Five/seven-block jobs retain nonidentity entries.
For each r,k use x=(1..k), every y in S_k, every q with q and q+1 units mod r,
voltage positions0,q,-q and voltages x,y,(yx)^-1. This is an exact finite
construction population, NOT all tuples or all tuples satisfying the theorem.

Classify all y before closure; export up to two representatives per (q,local
cycle type). Those representatives are deliberately selected, not random or
an exhaustive component census. The entire cheap qualifying population is
counted, but only selected representatives have state-count attempts. Different
passports prove different components; equal passports require state-key closure
deduplication. Inherited three-block/single-3-cycle passports are not new patterns.

Separate Rust target reuses existing exact permutation, group order and
six-move quotient closure. Construction is guaranteed smooth by the hypothesis
in MATH.md; all exported tuple/state assertions checked separately. A capped
closure leaves state count unresolved; it is not a closed-component certificate.
Human review of mathematical proof remains pending. Isotriviality is untested.

Default five single-CPU jobs, at most existing seven-slot budget. CPU-only:
classification is tiny and irregular closure has no demonstrated GPU speedup.
Never launch a competing search. Preserve old paused campaign. Freeze sources,
build/regressions, binary and known positive controls; share campaign lock.
T7 storage with30GiB guard; per-job600seconds, per-closure3seconds/10000states.
Independent literal audit of every selected tuple, every saved state and all
six edges of closed components; group-order validation uses SymPy separately.
Run until finite jobs and their audits finish; no ETA to geometric novelty.
