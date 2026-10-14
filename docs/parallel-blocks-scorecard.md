# Block refactor scorecard

Updated 2026-10-14 after independent review of job 55504.
**Stage 266 remains accepted. V7 is the successful, cluster-qualified structural
repair candidate for this fixture; legacy parity, memory parity, whole-timestep
independence and useful scaling remain unfinished.**
Objectives: [authoritative reset](parallel-blocks-objectives.md).

Active development branch: `codex/stage278-structural-repair`, based directly on
accepted Stage266. [Branch roles](parallel-blocks-branches.md) preserve the separate
historical experiments. The qualified predecessor and v7 source are committed;
this records the tested candidate without changing production acceptance.

## Source and campaign decisions

| Item | Current state |
|---|---|
| Accepted solver | Stage 266 `ff740744608731a05d53771c32a3055309f3a586` |
| Legacy control | Main `6a245db7f1632b0205c2bebe95e5f0cd535c2582`, established checkpoint I/O hardening |
| Stage 273/277 | Separate unpromoted partition/cache candidates; excluded from 278 |
| Stage 278 v5 /55502 | Previous ownership candidate; science pass, performance superseded by v7 |
| Attribution 55503 | 12,099 files, 37 passing runs, 38:06; construction diagnosis informed v7 |
| Structural v7 /55504 | 15,672 files verified; five builds, 48 passing runs, 42:33; exact six-file repair delta |
| Fallback/cache proposals | Excluded; fallback lacks standalone case, cross-epoch route-cache probe had zero hits |
| Next campaign | None prepared/submitted during review; choose next architecture repair from remaining dependencies and storage costs |

Implementation: `/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic`. Working source matches the delivered v7 candidate.
Controls and all prior sealed sources/evidence are preserved. Retain v7 as next
development base; this review makes no solver source edit or automatic promotion.

## Dual-goal gates

| Gate | Evidence / remaining work |
|---|---|
| Connected RHS and real RK ownership | Pass in candidate; all 96 zero-Domain ranks do numerical work, 16,212,000 primitive events at 256 |
| Paired retirement | Domain compatibility RHS, scalar capture/transport, velocity remainder export, full temperature-stage publication and displaced interior RK remain retired; v7 adds no numerical retirement |
| Science/lifecycle | All 48 runs pass; exact accepted-family/main policy, 12 states, 2 reloads, 2 remaps, 62 VTK members; checked/oracles retained |
| Whole-timestep independence | Open: boundary/scaffold/RK adapters, inverse/lifecycle publication, split Simple physics, geometry/workspaces |
| Absolute performance vs 266 | Intermediate gate surpassed on this fixture: v7 is 16.32%/8.09%/8.00% faster at 83/160/256; every candidate median below every 266 control median |
| Legacy P<=D performance | Open: v7 remains 17.03%/38.45% slower than main at 83/160; further 14.55%/27.77% candidate time reductions required |
| Intermediate P>D scaling | Unresolved:160→256 slowdown 7.72%, speedup 0.9284, relative efficiency 58.02%; previous 61.48%, 26658.08% |
| Memory | Improved 9.84%/18.97%/24.22% vs previous, but 16.00%/38.67%/42.42% above 266 summed rank peak RSS |
| Useful P>=4 D | Unqualified; J5 requires >=640 matched multi-node ranks and sufficient subdivision/load balance; no invented efficiency threshold |

Scope: restarted compressible non-mode-split RK4 climate/Simple physics, no
NCAR/SSO, established short J5 fixture on nk/bb13. Ownership credit is for the
qualified slice, not an independent whole timestep or all solver configurations.

## Latest matched ordinary-step performance

Milliseconds, mean of four previous/v7 and two main/266 run medians. Ordinary
steps exclude startup/reload/saves and retain remaps; steps are not independent
repetitions. Production saves remain rare(O1000 steps), remaps about every five.

| Ranks | Main | 266 | Previous | V7 | V7 vs previous | V7 vs 266 | V7 vs main |
|---|---:|---:|---:|---:|---:|---:|---:|
| 83 | 540.900 | 756.425 | 790.529 | 632.989 | -19.93% | -16.32% | +17.03% |
| 160 | 365.127 | 550.008 | 598.310 | 505.522 | -15.51% | -8.09% | +38.45% |
| 256 | — | 591.899 | 608.270 | 544.532 | -10.48% | -8.00% | — |

Every adjacent bare A/B pair favors v7. Ordinary means improve 19.04%/14.62%/
10.96%; remap and non-remap groups improve at each count. Gains are supported
within this allocation; long representative workloads and multi-node performance
are not qualified. Remaining main gaps 92.09/140.40 ms cannot be declared solved
by the now much smaller preparation budget.

## Attribution, scaling and memory

Route construction falls 151.51→37.33, 92.41→33.57, 79.95→40.72 ms at 83/160/256;
input plans 33.83→14.28, 27.15→14.70, 31.02→19.15 ms. Named preparation self subtotal
falls 216.94→62.74,141.50→58.56,133.12→72.52 ms. RHS exchange is almost unchanged
at 160/256. This supports structural construction/hash repair rather than a
transport/batching explanation. Individual mechanisms were not ablated.

V7 is faster at 256 than before, but gains more at 160:92.79 vs 63.74 ms saved.
This increases the 160→256 median gap from 9.96 to 39.01 ms. The scaling failure
persists; neither better absolute time nor correct zero-Domain work closes it.
Candidate ordinary-mean scaling also fails(speedup 0.9201). Coarse window growth
includes inverse wait, RHS exchange, adaptation, native wavelet, migration and
catalog/route construction. Wait time includes late producers; causality is not
isolated. Coarse mean perturbation reaches 9.26% at 256; do not partition bare
median growth with these profiled means or omit unprofiled logged-step time.

| Ranks | Previous RSS GiB | V7 RSS GiB | 266 RSS GiB | V7 vs previous | V7 vs 266 |
|---|---:|---:|---:|---:|---:|
| 83 | 27.65 | 24.93 | 21.49 | -9.84% | +16.00% |
| 160 | 54.04 | 43.79 | 31.58 | -18.97% | +38.67% |
| 256 | 80.89 | 61.30 | 43.04 | -24.22% | +42.42% |

RSS is summed individual rank peaks, not concurrent job peak. Lower memory is
measured; MPI/solver allocation attribution and long-run growth are unresolved.
All 1,094 integer work/traffic/narrow-storage rank files exactly match 55503.
Existing v7 source lifetime limits and metadata payload tests remain valid.

## Load balance and granularity evidence — job 55504 follow-up

[Detailed derivation](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/review-55504/load-balance-followup.md),
[rank statistics](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/review-55504/load-balance-followup.json). The same 71,052,000 trajectory
primitive events are spread more unevenly at256: maximum/mean rises1.24→1.67.
Only640 blocks exist; average blocks/rank falls4→2.5, and122 of256 ranks have one
final block. The maximum rank ordinary RHS arithmetic mean barely changes:
33.78→33.51 ms. At256 Domain owners average342,750 primitive events per rank,
zero-Domain owners168,875; the latter wait more in RHS exchange and inverse.
These are explicit load-distribution/granularity signals, not just transport cost.
Counters cover the validation trajectory; coarse timings cover ordinary steps.

The next repair assessment must explicitly separate assignment quality, unequal
indivisible block size and insufficient useful work relative to fixed overhead.
The current next-fit subtree-weight model does not balance measured whole-step
work or retained Domain dependencies. Evaluate finer subdivision and cost-aware
assignment separately with locality/halo/memory costs, without automatically
adopting rejected273/277 candidates. No new timing or source change is implied.

## Next deliverable and blockers

- Retain v7's demonstrated gain and ownership retirements; do not repeat the same
  campaign merely to reconfirm a large benefit.
- Investigate the inherited inverse/lifecycle and RHS execution dependencies,
  remaining storage and rank-count growth using existing evidence before selecting
  the next repair. Map actual producer/readers and lifetimes, then identify the work,
  dependency or allocation removed and a defensible complete-step benefit.
- Avoid another metadata-only plan as the whole answer to main parity. At 160 the
  remaining preparation subtotal 58.56 ms is much smaller than the 140.40 ms bare
  main gap (different statistics; not an additive forecast).
- Require an architecture repair that addresses 256 scaling and memory while
  preserving ordinary absolute performance. Generic batching/shared memory and
  count-only cross-epoch reuse remain unsupported causal shortcuts. Any extra
  instrumentation needs a named unresolved implementation decision.
- Finish whole-timestep ownership, granularity/load balance and matched >=640
  qualification as separate required gates. Stage 266 remains accepted.

## Evidence and history

[55504 full report](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/docs/parallel-blocks-stage278-job55504-results.md),
[verified return](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/review-55504/verified.json), [metrics](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/review-55504/metrics.json),
[prior scorecard](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/review-55504/scorecard-before-review.md),
[runtime ledger](parallel-blocks-campaign-runtimes.json).
[55503 attribution report](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/docs/parallel-blocks-stage278-job55503-results.md)
and [v7 local qualification](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/docs/parallel-blocks-rhs-structural-repair-278-v7.md)
retain prior evidence and limitations. Original returned payloads are not rerun;
raw reports/hashes are independently verified. No remote actions during review.
