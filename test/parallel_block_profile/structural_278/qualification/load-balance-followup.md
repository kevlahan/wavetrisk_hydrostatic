# Job 55504: load balance and granularity follow-up

Derived 2026-10-14 from the already verified return; no new run or source edit.
Primitive/scalar/velocity counters below cover the twelve-step validation trajectory;
profile timings cover the nine ordinary steps in the separate coarse runs.
They are work proxies and phase observations, not a causal partition of full time.

| Measure | 160 ranks | 256 ranks |
|---|---:|---:|
| Final catalog blocks | 640 | 640 |
| Average blocks/rank | 4 | 2.5 |
| Ranks with exactly one final block | 34 | 122 |
| Total primitive events | 71,052,000 | 71,052,000 |
| Mean primitive events/rank | 444,075 | 277,546.875 |
| Maximum primitive events/rank | 552,120 | 464,400 |
| Maximum/mean primitive events | 1.2433 | 1.6732 |
| Maximum/mean scalar targets | 1.2359 | 1.6712 |
| Maximum rank ordinary mean RHS arithmetic, ms | 33.782 | 33.512 |

At 256, the 160 Domain-owning ranks average 342,750 primitive events each; the
96 zero-Domain ranks average 168,875. All do numerical work, but equal participation
is not equal work. Zero-Domain ranks average 155.85 ms RHS exchange versus 139.96
ms for Domain owners, and 68.94 versus 43.17 ms inverse wait. These observations
are consistent with late producers and work imbalance; waits also contain actual
communication and cannot isolate producer lateness. Whole-trajectory counters do
not establish per-phase/per-step critical-rank imbalance.

The mean work per rank falls by the ideal 37.5% for a fixed total, but maximum
primitive work falls only 15.9%, and the maximum rank's observed ordinary RHS
arithmetic time barely falls. This is evidence for a load distribution/granularity
problem alongside fixed per-rank overhead, not proof that all ranks have too little
useful work or that communication transport is dominant.

The current catalog splits each Domain into at most four root-child subtrees and
uses next-fit assignment of subtree weights, relaxing the imbalance allowance
until assignment succeeds (parallel_block_mpi.f90:4923, 4990, 5250). The RHS follows
block owners, assigning boundary/scaffold work to a deterministic first block
(parallel_block_rhs_placement.f90:10). This does not balance measured complete-step
costs, retained Domain work or producer dependencies. At 640 ranks this catalog
would average only one unequal indivisible block/rank; further subdivision is an
architectural requirement for flexible balancing, not a demonstrated speedup.

Next investigation must distinguish three mechanisms: inadequate assignment of
existing blocks; indivisible blocks that are too coarse; and a fixed-cost/latency
floor as useful work per task falls. Use actual operation counts, phase timing,
retained Domain work and placement/locality costs. Compare weighted assignment
and finer subdivision separately, preserving the known rejected Stage273/277
history rather than importing those candidates automatically. Do not multiply
1.67 into a claimed solver speedup or infer network dominance from MPI waits.

Commit recommendation: preserve the measured v7 state now on a candidate branch,
without declaring it the accepted production baseline. HEAD is still Stage266 and
prior Stage278 files remain uncommitted/untracked, so record a coherent source
snapshot using sealed manifests. Prefer a tested55503 predecessor commit followed
by the six-file v7 repair/evidence commit. Do not blanket-stage raw archives,
binaries or unrelated work. This follow-up makes no Git commit or promotion.
