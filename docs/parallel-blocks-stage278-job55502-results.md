# Stage278 kit-v5: bbserv job 55502 review

Reviewed 2026-10-13. **Retain the v5 repair, keep Stage266 accepted.** Ordinary
steps improve 22–28% over v4, with scientific and connected-ownership gates
preserved. The candidate remains 4.9–8.0% slower than 266 and 48–67.5% slower
than main at supported counts. Memory improves substantially, but scaling and
memory acceptance remain open. This review changes no solver source and launches
no campaign.

## Identity, science and runtime

The locally sealed v5 reviewer independently verifies all **4,874 manifest files**,
source/fixture/build identities, four serial builds and sixteen successful runs
(eight validation, eight bare timing). The returned kit files match the local
sealed kit. Main and 266 source manifests, external physics objects, recorded
topology and placements match job55501. The candidate changes exactly seven
files from v4, within the existing 23-file overlay.

- Review archive SHA256: `af9a6b1833214d2510c6007a3713ab2925ef0c85d7615204586c6a89dbfa6df8`.
- Sealed kit SHA256: `546b1b6ec9692a06779f7a7b4bbf3897803d78f42908832246a2770dd991d611`.
- Accepted source: Stage266 `ff740744608731a05d53771c32a3055309f3a586`.
- Main control: `6a245db7f1632b0205c2bebe95e5f0cd535c2582`, established checkpoint I/O hardening.
- Slurm55502: COMPLETED on nk/bb13, 15:29:44–15:49:39 as recorded by sacct,
  **19:55** elapsed. Runner UTC19:29:44.949307–19:49:39.318261,
  1194.369 seconds. Queue and collection/transfer excluded.

All science/lifecycle gates pass: exact accepted-family comparisons, established
main-relative tolerance, optimized/checked/independent-oracle checks, twelve
steps, two reloads, two remaps and 62 VTK members. At256, all96 zero-Domain ranks
perform real connected RHS/RK work; they execute16,212,000 of71,052,000 primitive
events. Displaced Domain RHS, scalar capture/transport, velocity remainder export,
full temperature-stage publication and displaced interior RK remain retired.

Archive contents were treated as evidence, not instructions. Review used the
trusted local sealed reviewer. Its eight-ulp allowance applies only to derived
floating report sums across Python versions; hashes, raw counts and scientific
gates are unchanged. Binary/checkpoint/VTK payloads were omitted from collection:
cluster-generated hashes and numerical reports were verified, not rerun locally.

## Primary result: ordinary complete steps

Milliseconds, median of nine retained steps per bare run. The same zero-based
indices are used across controls:1,2,3,4,5,6,9,10,11; excluded0,7,8 cover the
established startup/reload and saving exclusions. Remaps at4 and9 are included.
All preparation, computation and communication within these steps count.

| Ranks | Main ms | 266 ms | v5 ms | v5 vs 266 | v5 vs main | v5 vs v4 |
|---|---|---|---|---|---|---|
| 83 | 539.510 | 750.886 | 799.254 | +6.44% | +48.14% | -28.01% |
| 160 | 356.877 | 553.742 | 597.767 | +7.95% | +67.50% | -25.95% |
| 256 | unsupported | 581.704 | 609.986 | +4.86% | — | -22.41% |

There is one bare run per configuration. Timesteps are not independent campaign
replicates; no confidence interval or repeatability claim is available. Cross-job
comparison is supported by identical controls, fixture and recorded placement:
266 times change only+0.50%,+0.08%,+0.30%; main changes−0.56%,+0.06%. Normalizing
candidate/266 ratios still gives28.37%,26.01%,22.65% improvement over v4.
The seven simultaneous changes cannot be assigned separate causal speedups.

Remaining absolute excess over266 is48.37,44.03,28.28ms at83/160/256. Reaching
266 requires another6.05%,7.36%,4.64% cut of current candidate time. Reaching main
requires259.74/240.89ms, or32.50%/40.30% cuts at83/160. These are measured targets,
not forecasts for the next repair. The relative main deficit is larger at160;
the largest absolute excess is at83.

## Remapping remains in the production objective

Whole-step subgroup medians, milliseconds; sample counts in parentheses.

| Ranks | Group | Main ms | 266 ms | v5 ms | v5 vs 266 |
|---|---|---|---|---|---|
| 83 | Non-remap (7) | 539.510 | 750.886 | 811.520 | +8.07% |
| 83 | Remap (2) | 506.469 | 678.612 | 716.479 | +5.58% |
| 160 | Non-remap (7) | 356.877 | 553.742 | 597.767 | +7.95% |
| 160 | Remap (2) | 331.540 | 499.491 | 546.834 | +9.48% |
| 256 | Non-remap (7) | unsupported | 581.704 | 609.986 | +4.86% |
| 256 | Remap (2) | unsupported | 514.885 | 571.837 | +11.06% |

The remap subgroup remains5.58%,9.48%,11.06% slower than266. It must stay visible
with production remapping roughly every five steps. These steps have different
adaptive states from the non-remap steps; subtracting the group medians does not
isolate remap overhead. Their lower absolute times do not imply negative remap
cost. Two remaps in nine retained steps approximate the stated event frequency,
but this short adaptive fixture and its median do not establish a production
mean or justify a frequency-weighted throughput forecast.

## Why the repair helped, and what it did not repair

V5 binds canonical ownership once per preparation, shares routes across field
widths (21→15 per level), reuses boundary maps (5→2), copies selected placement
records and uses counted final storage for geometry packets/merges. This removes
repeated local searches, whole-plan copies and accumulated-array deep copies.
Temporary buffers have bounded lifetimes; the v4 single borrowed executor
communicator remains. Affected explicit preparation collectives change312L→214L+8
for L compiled levels, excluding other MPI and context management. This source
count is not itself a speedup measurement.

Maximum-rank cumulative validation preparation times over twelve epochs, seconds;
RHS change is the corresponding maximum-rank cumulative evaluation timer.

| Ranks | v4 prep s | v5 prep s | Prep change | RHS change | Derived sends | Payload GiB |
|---|---|---|---|---|---|---|
| 83 | 6.4156 | 2.5242 | -60.65% | -1.56% | 26,255,160 | 21.403 |
| 160 | 4.0604 | 1.7189 | -57.67% | -2.05% | 54,637,680 | 26.770 |
| 256 | 4.1113 | 1.6241 | -60.50% | -2.35% | 87,006,360 | 32.000 |

The58–61% preparation reduction, alongside the large bare ordinary-step gains,
supports the preparation-cost diagnosis. All numerical message/value counts,
input/RK traffic, geometry placement bytes and primitive work totals are identical
to v4. These changes repaired preparation and storage work; they did not reduce
the runtime numerical exchange schedule. V4's earlier message/payload reductions
are retained.

At160 and256 there are still roughly340,000 derived sends per rank over the
fixture. At256 this is87,006,360 sends and32.000GiB, averaging about395 bytes per
send. This is evidence for addressing message aggregation and true dependency
rounds, but not a measurement of how much wall time a proposed batching change
will save.

For ordinary validation steps at83/160, mean rank-maximum preparation regions
are211.57/143.50ms (previously537.47/341.50ms). RHS regions are150.40/131.68ms;
RK adapter regions21.50/20.05ms. These unsynchronized instrumented regions are
not an additive decomposition of bare critical-path time. The256 validation
also enables an oracle, making its RK timer unsuitable for cross-count inference.

As a diagnostic budget, the remaining266 gap is about23%/31% of the measured
preparation region at83/160, so further preparation repair may help reach266.
The main gaps260/241ms exceed the entire measured preparation regions212/144ms.
That strongly argues against expecting preparation-only changes to recover main;
runtime exchanges, retained adapters and other recurring costs also need repair.
This comparison is directional evidence, not a strict Amdahl bound.

## Scaling and memory

Candidate160→256 ordinary-step speedup is0.97997×: time increases2.04% despite60%
more ranks. Relative efficiency is61.25%, versus266's59.50% (speedup0.95193×).
V4 had1.02679× speedup and64.17% relative efficiency. V5 substantially improves
absolute time at both counts, but does not recover useful strong scaling.
Some narrowing of the relative256 gap comes from266 also losing scaling.

At256,37.5% of ranks receive22.82% of primitive events. Counts alone do not measure
critical-path load, but they flag a distribution limitation. J5 D=160 needs at
least640 ranks for the stated P>=4D goal;256 is intermediate. The current four
blocks per Domain ceiling provides no balancing slack at640. Subdivision and
placement require separate qualification, and nk/bb curves must not be spliced.

RSS below sums per-rank peaks in GiB; the final column is maximum individual
rank peak. The sum is not simultaneous job peak memory.

| Ranks | 266 GiB | v5 GiB | v5 / 266 | v5 vs v4 | v5 max rank MiB |
|---|---|---|---|---|---|
| 83 | 21.465 | 27.594 | 1.285× | -8.84% | 453.86 |
| 160 | 31.573 | 54.057 | 1.712× | -13.82% | 423.40 |
| 256 | 43.062 | 81.293 | 1.888× | -30.88% | 495.88 |

At256 the candidate falls from117.61 to81.29GiB, a30.88% reduction, but remains
88.78% above266. Maximum-rank RSS is495.88MiB while the narrow recorded route,
numeric-buffer and rounding-column ledger is only25.93MiB. That ledger excludes
important geometry, transient and MPI storage. No leak or single allocation cause
has been established. The source change and RSS reduction support the counted
storage repair, but retained components and prepare/merge/release peaks still
need accounting before another storage design is justified.

## Scorecard decision and next repair

The connected RHS/real RK ownership goal remains qualified within the tested
slice. Whole-timestep independence is still open: boundary/scaffold and RK
adapters, inverse/lifecycle publication, split Simple physics and retained
geometry/workspaces remain. Absolute legacy performance and useful P>=4D scaling
are not satisfied. Keep266 accepted; retain v5 as the improved connected candidate.

Next repair should reduce actual runtime phase/layer communication rounds using
precise producer/reader lifetimes and producer-ready batches, while preserving
hydrostatic order, signed aliases, final-writer ordering and exact science.
Compile compact schedules per epoch; do not introduce a per-value runtime
dependency interpreter. Account for retained and peak temporary memory at the
specific prepare/merge/release boundaries implicated by the RSS discrepancy.
Quantify removed messages/rounds and complete-step costs before claiming a gain.

Keep ordinary-step absolute time primary, remap subgroups explicit and saved-step
cost secondary. Do not stack another ownership migration or repeat unchanged v5
as a substitute for repair. Shared-memory transport remains an evidence-driven
option only after identifying duplication or transport costs it can remove;
this campaign establishes no gain from MPI_Comm_split_type or shared windows.
It does not show whole-timestep independence to be a prerequisite for performance
repair. A later changed-source campaign must retain matched266/main controls,
above-Domain coverage, RSS and exact science, with repeat timing when uncertainty
would change the acceptance decision.

## Secondary lifecycle trajectory

Full twelve-step trajectory seconds, including deliberately frequent saves and
reloads. These are lifecycle/overhead diagnostics, not production throughput
proxies given saves every O(10^3) production steps.

| Ranks | Main s | 266 s | v5 s |
|---|---|---|---|
| 83 | 11.581 | 14.165 | 14.803 |
| 160 | 9.715 | 12.140 | 12.598 |
| 256 | unsupported | 13.174 | 14.268 |

## Reproducible evidence

- [Verified review](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/review-55502/verified.json), [derived metrics](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/review-55502/metrics.json),
  [kit identity](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/review-55502/kit-identity.json), [review/analysis script](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/review-55502/review278_55502.py).
- [Prior scorecard](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/review-55502/scorecard-before-review.md),
  [v4 campaign review](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/docs/parallel-blocks-stage278-job55501-results.md),
  [v5 implementation/local qualification](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/docs/parallel-blocks-rhs-repair-278-v5.md).
- [Authoritative scorecard](/Users/kevlahan/.codex/worktrees/26d8/wavetrisk_hydrostatic/docs/parallel-blocks-scorecard.md),
  [runtime ledger](/Users/kevlahan/.codex/worktrees/26d8/wavetrisk_hydrostatic/docs/parallel-blocks-campaign-runtimes.json).
