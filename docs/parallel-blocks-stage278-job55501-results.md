# Stage278 job 55501: useful repair, performance acceptance still fails

Reviewed 2026-10-13. **Keep Stage266 accepted. Retain kit-v4's qualified repair,
but do not promote it or resume a broader ownership migration yet.** The repair
preserves connected ownership and reduces communication substantially. Its
absolute performance, resource cost and scaling still miss the dual goals.

## Verification

Archive SHA256 `acce384f39454d802b12f2b8d8eb00afa1f5873b0fbad25378147ec313474937`.
All 4,874 inventoried files verify, and every returned kit file matches sealed
kit-v4 exactly. Its five solver-file changes relative to kit-v3 are explicitly
recorded in [metrics](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/review-55501/metrics.json); main and266 source controls are identical
between the two jobs. All four serial builds and 16 runs pass. Recorded topology
and 83/160/256 placement plans are identical to job 55500 on bb13.
No repeat group was used; each timing configuration has one bare run.

Slurm: **12:56:00–13:21:37 America/Toronto, elapsed 25:37**, October 13.
Raw sacct has no explicit timezone; runner UTC 16:56:01.317001–17:21:37.223434
corroborates the local display. Runner elapsed 1535.906454s; queue and collection
are excluded. The runtime ledger records actual build/run schedule and accounting.
No remote action or new job was performed during this review.

Every run completes 12 steps, two reloads and two remaps, with exact accepted-family
checkpoint/numerical gates, established main-relative tolerance and 62 VTK members.
Checked and independent-oracle gates pass. All 96 zero-Domain ranks at 256 execute
real authoritative RHS/RK work: 16,212,000 primitive events, 23,040–438,240 per rank,
and 48 retention transactions each. Primitive and input/RK traffic totals are
unchanged by the repair. Whole-timestep independence still excludes named
boundary/scaffold/RK adapters, inverse/lifecycle publication, separate split Simple
physics and retained geometry/workspaces. This remains the supported restarted
compressible non-mode-split RK4 climate, Simple physics, no NCAR/SSO fixture.

The review rederives reports from raw logs/counters and verifies cluster-generated
binary/checkpoint/VTK hashes. Those binary payloads are omitted from the review
archive, so it is not an independent solver rerun or a recomparison of payloads.
Eight-ulp tolerance is confined to derived floating report sums across Python
versions; scientific tolerances, integer counts, source and raw member hashes
remain unchanged. [Verified report](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/review-55501/verified.json).

## Absolute performance and scaling

Ordinary-step medians, milliseconds, with the same step0/7/8 exclusions:

| Ranks | Main | Stage266 | Repair278 v4 | v4 vs266 | v4 vs main |
|---:|---:|---:|---:|---:|---:|
| 83 | 542.531 | 747.157 | 1110.224 | +48.59% | +104.64% |
| 160 | 356.657 | 553.278 | 807.274 | +45.91% | +126.34% |
| 256 | unsupported | 579.960 | 786.209 | +35.56% | — |

The candidate still needs runtime cuts of 32.70%,31.46%,26.23% at 83/160/256 merely
to match 266, and 51.13%/55.82% at 83/160 to match main. These are break-even
requirements, not forecasts. Complete12-step trajectories also regress versus 266:
14.187→18.521s (+30.55%),12.062→14.885s (+23.40%),13.210→16.398s (+24.13%).
Their improvement versus v3 is smaller than the ordinary-step improvement:
2.24%,3.38%,4.02%. At256 the launcher elapsed is68.636s versus64.175s for 266;
launch-level timing includes costs outside the measured step trajectory.

Candidate 160→256 speedup is **1.02679×**, relative efficiency 64.17%;60% more
ranks yields 2.68% speedup. Whole-trajectory time instead rises 10.16%.
The prior candidate was 1.00341× and 62.71%; matched266 is 0.95399× and 59.62%.
These slightly better ratios cannot compensate for a much slower absolute time.
256 is intermediate for J5 D=160; useful >=640 execution remains unqualified.

## What the repair established

Changes between job 55500 v3 and55501 v4. The normalized column compares the
candidate/266 time ratio in the two jobs; it is a control-drift check, not a
statistical correction or substitute for a same-allocation repeated A/B experiment.

| Ranks | Candidate time | 266 control time | Candidate/266 ratio | Derived messages | Derived payload |
|---:|---:|---:|---:|---:|---:|
| 83 | -2.52% | -0.21% | -2.31% | -8.68% | -29.18% |
| 160 | -6.26% | +0.51% | -6.74% | -11.87% | -36.92% |
| 256 | -8.40% | -1.09% | -7.39% | -17.76% | -39.52% |

Main control changes are−0.67%/+0.66% at 83/160. Together with identical recorded
placement, unchanged controls/fixture and unchanged numerical work, these results
support an improvement from the combined repair, especially at 160/256. One run
per configuration in separate jobs cannot establish confidence bounds or isolate
how much comes from demand pruning, field batching or communicator reuse.
Do not extrapolate the improvement linearly or treat it as recovered parity.

Current derived-stream totals over the whole validation trajectory:

| Ranks | Sends | Payload GiB | Mean payload bytes | Max-rank preparation s | Max-rank RHS evaluation s |
|---:|---:|---:|---:|---:|---:|
| 83 | 26,255,160 | 21.403 | 875.3 | 6.416 | 1.811 |
| 160 | 54,637,680 | 26.770 | 526.1 | 4.060 | 1.603 |
| 256 | 87,006,360 | 32.000 | 394.9 | 4.111 | 1.778 |

At256 the remaining derived stream is 87,006,360 sends,32.000GiB and only395 bytes
per send; it accounts for99.63% of recorded connected stage messages. Total
input+derived+RK traffic is 87,326,416 sends and40.740GiB, plus 8-byte epoch headers
per send. Counts omit other solver MPI and topology metadata. At160 and256 mean
derived sends per rank are341,486 and339,869: almost unchanged despite more ranks.
This is evidence of a substantial per-rank message floor as arithmetic is divided
more finely, consistent with the remaining phase/layer schedule limiting scaling.
It is not direct critical-path attribution.

Narrow max-rank RHS evaluation falls 18–21%, whereas preparation falls only 0.5–1.5%.
The256 validation uses oracles, unlike83/160, and these timers exclude other work;
do not add rank maxima or subtract them from bare step medians. The repair removes
per-route context multiplicity and combines four thermodynamic transactions
into two, but retains layer-by-layer synchronization and conservative cross-level
field demand, including canonical-owner refreshes. More aggressive aggregation
and precise phase liveness remain possible, subject to dependency/epoch proof.

## Memory result

Solver child wait4 measurements from the bare timing runs:

| Ranks | 266 sum of rank peaks GiB | v4 sum of rank peaks GiB | v4/266 | v4 vs v3 |
|---:|---:|---:|---:|---:|
| 83 | 21.46 | 30.27 | 1.41× | -0.25% |
| 160 | 31.57 | 62.72 | 1.99× | -0.86% |
| 256 | 43.05 | 117.61 | 2.73× | -0.97% |

The maximum256-rank process RSS is 595.41MiB versus269.00MiB for 266. The new
component ledger maximum is only 26.22MiB/rank; it excludes geometry/retained
workspaces, descriptors, compilation transients and MPI internal storage.
Summing individual rank peaks is not a simultaneous allocation peak. Major faults
do not alone establish swapping. **The communicator change has not repaired the
memory regression.** It does not prove a leak or identify which excluded component
dominates; little change in observed RSS is evidence against claiming context
multiplicity alone explains the previous increase.

## Decision and next work

1. Retain the direct block-owned numerical implementation and qualified repair as
   an unpromoted candidate. Keep266 accepted. Ownership/retirement has passed;
   absolute performance, memory and useful scaling have not. No new Domain
   migration or partition change should be stacked onto this regression.
2. Reduce the remaining communication rounds and unnecessary consumers. Use exact
   per-phase reader/write lifetimes and explicit producer-ready groups to batch
   compatible fields/layers/peers. Preserve hydrostatic column order, final-write
   alias semantics and exact science. Compile that schedule once per epoch;
   do not recreate the per-value runtime dependency interpreter. The87-million
   small-message residual is the concrete target; payload reduction alone is
   demonstrably insufficient.
3. Resolve the memory gap with a bounded validation-only measurement at prepare,
   execution and release boundaries, accounting for omitted retained allocations
   and compiler peaks. This has a named decision: distinguish persistent redundant
   storage from transient/high-water allocator or MPI effects before choosing
   what to release or redesign. Do not promise an RSS saving from communicator
   lifetime changes without evidence.
4. Keep explicit shared windows as a separate transport experiment. Replacing
   MPI envelopes could help intra-node overhead, but no result here establishes
   that gain, and it would retain the numerical dependencies and inter-node path.
   Require within-NUMA/cross-NUMA/inter-node comparison and correct visibility,
   epoch completion and local first-touch. Do not disguise a phase barrier as a
   shared-memory scaling fix.

A local 1/2/3/4-task comparison remains useful for fixed overhead and memory effects,
but all those counts are below J4 D=40. The returned cluster evidence already
establishes the need for another schedule/memory repair; an unchanged broad
cluster repeat would not change the acceptance decision. Work balance is also
unchanged:96/256 ranks (37.5%) receive22.82% of primitive events;640 blocks provide
no balancing slack at 640 ranks. Subdivision/load balance must be separately
qualified, not mixed into the communication repair. No numerical final efficiency
threshold has been agreed.

[Authoritative scorecard](/Users/kevlahan/.codex/worktrees/26d8/wavetrisk_hydrostatic/docs/parallel-blocks-scorecard.md).
[Raw archived evidence](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/review-55501/review/results/results.json).
