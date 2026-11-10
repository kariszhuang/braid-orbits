# Forest acceleration: pre-specified development and benchmark plan

Objective: faster identity-free, positive-genus everywhere-smooth search without
changing literal smoothness or simultaneous-S_n state equivalence. This is a
biased search, not a degree census. The live degree10 audit and its two frozen
guided lanes must continue unchanged during development.

Implement a separate target retaining the old implementation as an oracle:
fixed first-pair factors, direct composition-index calculation, feasible third
indices, bounded frame reuse with regularly fresh frames, allocation-free third
construction, compatible canonical keys, early nonsmooth successors, carried
entry indices, persistent lead handles and interval resolution flushing.
Always flush a seed before closure; a missing resolution remains recoverable.
Caps/timeouts remain unresolved. Independent Python audit remains mandatory.

Validation: every product-one S4 triple; randomized larger-degree comparisons;
known closed positive-genus controls; old/new complete small closures; literal
negative-word replay; identical canonical keys; all feasible-index conditions
checked against the literal passport rule. Reference proposal mode keeps the
same candidate stream for end-to-end finite comparisons.

Benchmark: fixed deterministic candidate traces, same degrees11--14, paired
alternating old/new order, warm-up and seven repetitions; publish medians,
paired speed ratios and descriptive bootstrap intervals. Measure screen,
pre-closure pipeline, closures and T7 I/O separately. Do not multiply these
ratios or claim faster discoveries. Shared-machine results are not idle-machine
SLAs. Feasible sampling/frame reuse changes proposal weights, reported separately.
Bounded one-CPU tests may coexist with the audit, not a competing search pool.

GPU: test a batchable integer screen, include dispatch/readback and full CPU
shadow validation; activate no GPU discovery lane absent a verified end-to-end
benefit. Do not reuse the weak-search GPU predicate for smoothness.

Deployment: freeze the validated successor and queue it behind the existing
shared campaign lock; never rewrite a running snapshot or interrupt the audit.
Keep six total CPU workers maximum (current4+2 budget), single-threaded child
producers and independent audits. Preserve all old capped candidates for replay.

Design/software acknowledgment: experimental-design procedures materially
informed paired blocking and run order. Timothy Kassis, Vinayak Agarwal, Yuhuan
He, Darshil Patel and Aubrey M. Brueckner (2026), *Scientific Agent Skills: A
Library of Procedural Knowledge for Research Agents*, arXiv:2609.00065 (latest
record checked2026-10-12; revised2026-09-02),
https://doi.org/10.48550/arXiv.2609.00065.
