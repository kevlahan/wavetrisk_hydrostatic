# Stage 278 job 55503: bottleneck attribution review

**Verified: all 37 runs pass; Stage 266 remains accepted.** The campaign locates
material costs in route preparation and exchange execution, with inherited inverse/
writeback costs contributing to the 256-rank problem. It does not establish network
transport as the dominant cause or justify blind communication batching.

## Integrity, science and runtime

The user supplied `rk278-review-55503.tar.gz`. The trusted local reviewer verifies
**12,099 files**, the full run schedule, raw profile partitions, timing, numerical
reports, placement, resource records and ownership gates. Returned source and kit
seals exactly match the delivered attribution kit. No instructions or executable
helpers from the returned archive were adopted or executed.

Five serial builds and 37 serial runs completed on nk/bb13 in **38:06** (2,286 s).
Slurm records start `2026-10-13T17:49:28`, end `2026-10-13T18:27:34`; runner UTC
records and the raw accounting timezone distinction are retained. GCC 13.2.0 and
the recorded MVAPICH environment built optimized and checked variants successfully.
All twelve states, exact accepted-family comparisons, established main-relative
policy, 62 VTK members, two reloads and two remaps pass. All 96 zero-Domain ranks
perform production work: 16,212,000 primitive events in the 256-rank validation.

Across corresponding validations, all integer ownership/work/traffic/storage
fields in **1,094 rank files** exactly match job 55502. This is attribution-only
source, not a new numerical repair. The canonical fallback proposal was excluded.
Binaries and checkpoint/VTK payloads are omitted from the returned archive;
cluster-computed hashes and numerical reports were verified, not rerun locally.

## Absolute performance and scaling

Ordinary complete-step milliseconds. Original v5 entries are means of two run
medians; main/266 have one control run per configuration. Retain the same nine
ordinary steps, including remaps; no confidence interval or independent-timestep
sample interpretation is claimed.

| Ranks | Main | 266 | Original v5 | v5 vs 266 | v5 vs main |
|---|---:|---:|---:|---:|---:|
| 83 | 540.186 | 753.402 | 797.111 | +5.80% | +47.56% |
| 160 | 357.787 | 546.473 | 599.258 | +9.66% | +67.49% |
| 256 | — | 581.370 | 611.384 | +5.16% | — |

The remaining v5 cuts for 266 parity are **43.71 / 52.79 / 30.01 ms** at 83/160/256
(5.48% / 8.81% / 4.91% of v5 time). Main still requires **256.92 / 241.47 ms** at
83/160 (32.23% / 40.30%). The 160-rank 266 control is about 1.31% faster than in
55502; original v5 is only about 0.25% slower. The change in the reported gap
therefore is not evidence of a new algorithmic regression.

Original v5 still slows from 599.26 to 611.38 ms: **160→256 speedup 0.9802×,
relative efficiency 61.26%**. Attribution-off is similar (0.9793×). Accepted 266
also slows, 546.47→581.37 ms (0.9400×). Thus the failure is not confined to the
new RHS path. Roughly 78–83% of the main-relative absolute gap is already present
in 266 in this campaign; RHS-only tuning cannot be presumed to close it.

## What the profiles establish

The table reports means of self times on each ordinary step's **same slowest
whole-step rank**, averaged over the two coarse runs. These are profiled means,
not the bare run-median statistic above. Do not subtract them as an exact bare
critical-path prediction or sum independent regional rank maxima.

| Coarse region, ms/ordinary step | 83 | 160 | 256 | 256 minus 160 |
|---|---:|---:|---:|---:|
| RHS route construction | 151.61 | 92.75 | 84.20 | -8.54 |
| RHS input plan construction | 34.36 | 26.98 | 31.67 | +4.69 |
| RHS arithmetic residual | 40.58 | 26.90 | 23.04 | -3.86 |
| RHS exchange residual | 132.39 | 127.17 | 150.93 | +23.77 |
| inverse MPI wait | 47.07 | 38.39 | 51.39 | +13.00 |
| Domain writeback MPI | 11.42 | 10.69 | 16.61 | +5.92 |
| native wavelet work | 44.05 | 28.22 | 34.01 | +5.79 |
| adaptation residual | 26.43 | 20.75 | 29.39 | +8.63 |

The coarse whole-step mean rises **30.66 ms** from 160 to 256. RHS exchange rises
23.77 ms and inverse wait 13.00 ms; reductions elsewhere offset some of their
combined increase. Route construction falls 8.54 ms and RHS work falls 3.86 ms.
The new RHS regions together grow about 12.53 ms; the remaining regions grow
about 18.13 ms. The measured scaling loss therefore spans both the connected
RHS and retained surrounding paths.

At 83 ranks, route construction is the largest individual self region, about
152 ms, compared with 41 ms for RHS arithmetic/workspace work. At 160, the route
region is 93 ms and RHS exchange 127 ms. At 256, RHS exchange reaches 151 ms,
about 24% of the profiled mean step. These are material repair budgets; optimizing
small arithmetic or accounting loops alone is not supported as a full-gap repair.
The route region also initializes dynamics/workspace storage and its MPI context;
it includes metadata MPI and local construction. RHS arithmetic includes copies/
checks. Exchange includes input/RK/metadata routes as well as the derived numerical
stream. These labels do not isolate pure compute or transport, or identify how
much of construction can safely be removed.

## Perturbation and imbalance limits

The attribution binary with profiling off differs from original v5 by
−0.71% / +0.25% / +0.34%. Original-v5 repeat spreads are 1.02% / 1.11% / 0.57%.
There is no evidence of a substantial disabled-hook regression or gain.
Coarse profiling increases observed run-median time by **2.17% / 2.50% / 2.84%**.
The single fine/off pairs differ by **9.13% / 14.17% / 13.44%**. These comparisons
include run variation; they are not exact measurements of profiler overhead.

Fine RHS completion-wait means are 90.74 / 84.77 / 123.70 ms. However, total fine
exchange-region time exceeds coarse by about 58–65 ms/step, close to the 62–68 ms
increase in whole profiled means. The fine breakdown is intrusive and cannot
supply an uninstrumented network-time estimate or a guaranteed saving.

At 256 ranks, coarse rank-group means show Domain ranks doing 26.38 ms of RHS
arithmetic/workspace work and zero-Domain ranks 19.40 ms. Their exchange means
are 144.46 and 159.23 ms respectively. Fine mean waits are 103.54 versus 128.69 ms.
This is consistent with uneven producer work and late arrival contributing to
waits. Aggregated rank times do not prove which producer delayed which consumer;
transport, placement and progress effects remain mixed with imbalance. Both rank
groups appear among slowest whole-step ranks. Do not equate waits with network
latency or conclude that moving all work to Domain ranks would help.

## Memory and epoch evidence

Mean sum of per-rank peak RSS for ordinary unprofiled controls, GiB. These sums
are not simultaneous job peaks.

| Ranks | 266 | Original v5 | Ratio |
|---|---:|---:|---:|
| 83 | 21.466 | 27.645 | 1.288× |
| 160 | 31.582 | 53.964 | 1.709× |
| 256 | 43.096 | 80.824 | 1.875× |

Memory acceptance remains open: about 80.82 GiB versus 43.10 GiB at 256.
No leak or single allocation cause is established. Numerical/storage counters
match v5 exactly; the narrow ledger cannot attribute total RSS.

The ownership audit records **12 RHS plan epochs for 12 steps**, each reused for
four RK stages. The measured source takes its RHS generation from the writeback
plan generation; any generation change or mask change triggers full preparation.
Preparation rebuilds the owner directory, direct/alias identities and per-level
phase routes. This gives a concrete construction/invalidation design to repair.
It does not prove that successive topologies, placements or masks are identical:
unchanged counts are insufficient evidence for cross-epoch cache reuse.

## Next implementation decision

1. **Prioritize structural route preparation for low-rank absolute time.** Separate
   geometry/ownership/alias identities from mask-dependent phase plans. Establish
   exact producer/reader and epoch keys, then avoid repeated construction only
   where dependencies are proven unchanged; otherwise use compact direct rebuilds.
   Keep communicator/storage lifetime bounded to those dependencies and release
   obsolete storage. A small canonical-fallback shortcut alone remains insufficiently
   supported. Routes alone would need roughly 29% / 57% / 36% reductions to cover
   the current 266 gaps at 83/160/256: a screening budget, **not a speedup forecast**.
2. **Investigate execution dependencies and placement across RHS and retained
   inverse/compatibility paths before selecting a transport rewrite.** The measured executor packs/posts, immediately
   waits, and installs each phase. Any removed dependency or readiness change must
   preserve the exact producer/consumer ordering and demonstrate less critical
   waiting; batching by message count is not evidence. The rank-group imbalance
   and inherited inverse/writeback increases must be considered before selecting
   transport changes. Do not repeat the previously rejected asynchronous/local-first
   designs on the strength of wait time alone.
3. **Qualify the next substantive candidate locally before another performance
   campaign.** Require exact science, unchanged ownership retirement, bounded
   memory/context lifetimes and an explicit measured-cost argument. Do not run a
   standalone campaign for the small fallback proposal or another broad profiling
   repeat. A targeted additional diagnostic is justified only if a named causal
   ambiguity blocks that implementation. Whole-timestep independence, shared
   memory and >=640-rank placement remain separate milestones, not assumed shortcuts.

There is enough evidence to begin this structural repair design. There is not
an evidence-based unconditional 10–15% forecast, a proved network-only bottleneck,
or a completed 256-rank repair. Both original goals remain open.

## Evidence

- [Independent verification](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/review-55503/verified.json), [derived metrics](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/review-55503/metrics.json),
  [raw campaign results](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/review-55503/review/results/results.json), [Slurm accounting](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/review-55503/review/accounting.json).
- [Measured preparation/invalidation source](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/review-55503/measured-source/src/parallel_block_rhs_backend.f90:48),
  [route construction](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/review-55503/measured-source/src/parallel_block_rhs_distributed.f90:71),
  [exchange execution](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/review-55503/measured-source/src/parallel_block_rhs_exchange.f90:475),
  [generation source](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/review-55503/measured-source/src/parallel_block_mpi.f90:17732).
- [Previous 55502 review](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/docs/parallel-blocks-stage278-job55502-results.md),
  [scorecard](/Users/kevlahan/.codex/worktrees/26d8/wavetrisk_hydrostatic/docs/parallel-blocks-scorecard.md), [runtime ledger](/Users/kevlahan/.codex/worktrees/26d8/wavetrisk_hydrostatic/docs/parallel-blocks-campaign-runtimes.json).

The original user archive and sealed kit are unchanged. This review submits no
job, transfers no files remotely and changes no solver source.
