# Stage 278 job 55504: structural repair v7 review

Reviewed 2026-10-14. **V7 is a demonstrated absolute-performance and memory
improvement, and exceeds Stage 266 performance at all three tested counts. It
has not repaired 160→256 scaling or reached legacy parity.** Retain v7 as the
qualified next development candidate; Stage 266 remains the accepted solver.

## Verification and scope

The user-supplied [original archive](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/review-55504/rk278-review-55504.tar.gz) has SHA 256
`3cb8181d252a20a23d8d9964af481543c0171ccfff22adbbafcfb1c6dc1b407b`. All 15,672 manifest members verified; all returned kit
members exactly match the sealed structural-kit-v 1 delivery (archive SHA 256
`78d70c48b8f89b2329108e0427c4d54c7240f5ffce9fb7c9f7c91b9af1d804c5`).
The trusted local reviewer independently recomputed reports from raw returned
logs/rows. No returned Python/shell helper or embedded instruction was executed.

Five serial GCC 13.2.0/MVAPICH 4 builds and all 48 serial runs pass on nk/bb13:
eight science/ownership validations, 24 repair/previous bare runs, ten main/266
bare controls and six coarse runs. Slurm completed 09:47:28–10:30:01 on 2026-10-14,
elapsed **42:33 (2,553 seconds)**. Runner timestamps are 13:47:28.786718–
14:30:00.670716 UTC; runner elapsed 2,551.884 s. Slurm times are retained exactly
as supplied without an offset; they align with the client Toronto time. Queue
and collection are excluded. The previous 45–60-minute allowance remains an
estimate, not rewritten to match this actual duration.

Delivered source is the exact locally qualified six-file v7 delta over the
attribution candidate measured by 55503. Stage 266/main controls are unchanged;
Stage 273/277 partition/cache work and the fallback shortcut remain excluded.
The previous control here is the 55503 attribution candidate with profiling off,
not the separately built original-v 5 control used in earlier tables.

All runs preserve twelve states, two reloads/two remaps, 62 VTK members, exact
accepted-family checkpoint fields/topology and the established main-relative
science limits. Checked and independent-oracle runs pass. All 96 zero-Domain
ranks at 256 perform numerical work:16,212,000 primitive events in aggregate.
All 1,094 integer production/cost rank-audit files exactly match 55503, including
work, numerical traffic and narrow storage counters. Floating timings are
excluded from that integer comparison; temporary metadata has separate local
payload tests. Existing Domain retirements remain active; this repair adds no
new numerical ownership or retirement claim.

The review archive omits binaries/checkpoint/VTK payloads. Their recorded hashes
and cluster-computed comparisons were checked; those payload comparisons were
not rerun locally. Tested scope remains restarted compressible, non-mode-split
RK4 climate/Simple physics without NCAR/SSO, on the established short J5 fixture.

## Ordinary-step performance

Milliseconds, **mean of run-level medians**: four bare repetitions for previous
and v7, two for main/266. Exclude established indices 0, 7, 8; retain ordinary remaps
4, 9. Nine ordinary steps/run are not nine independent repetitions. Lower is better.

| Ranks | Main | 266 | Previous | V7 | V7 vs previous | V7 vs 266 | V7 vs main |
|---|---:|---:|---:|---:|---:|---:|---:|
| 83 | 540.900 | 756.425 | 790.529 | 632.989 | -19.93% | -16.32% | +17.03% |
| 160 | 365.127 | 550.008 | 598.310 | 505.522 | -15.51% | -8.09% | +38.45% |
| 256 | — | 591.899 | 608.270 | 544.532 | -10.48% | -8.00% | — |

The improvement is well outside observed run variation: every v7 run median is
below every previous and 266 median at its count, and all 12 adjacent A/B pairs
favor v7. This is evidence across repetitions within one allocation, not a
confidence interval across independent jobs or a general workload guarantee.

| Ranks | Previous run medians, ms | V7 run medians, ms | Adjacent paired savings, ms |
|---|---|---|---:|
| 83 | 779.551, 788.568, 789.018, 804.978 | 630.683, 631.710, 634.891, 634.674 | 148.87–170.30 |
| 160 | 596.107, 600.162, 598.571, 598.401 | 505.677, 505.793, 505.652, 504.968 | 90.43–94.37 |
| 256 | 607.330, 609.223, 607.041, 609.488 | 542.507, 549.925, 540.454, 545.242 | 59.30–66.59 |

Ordinary means independently improve 19.04%, 14.62%, 10.96% at 83/160/256.
Remap-group means improve 16.76%, 14.46%, 9.32%; non-remap means improve 19.61%,
14.66%, 11.39%. The gain therefore is not confined to checkpoint steps, hidden
in a median or purchased with a remap penalty. Short-fixture trajectory elapsed
and child CPU totals include frequent validation saves/reloads; they are not
production throughput proxies.

V7 closes 63.1%/39.8% of the previous absolute main gap at 83/160. Remaining main
gaps are 92.09/140.40 ms, requiring a further 14.55%/27.77% reduction in v7 time.
The intermediate 266 performance target is surpassed on this fixture. Legacy
parity and representative long adaptive workload qualification remain open.

## Why the structural repair helped

Coarse self times on the same slowest whole-profile-window rank, averaged over
ordinary steps; one matched profile per source/count. These are attribution
measurements, separate from the bare-run medians above.

| Ranks | Route construction, ms | Input-plan construction, ms | Retention plan, ms | Transform preparation, ms |
|---|---:|---:|---:|---:|
| 83 | 151.51 → 37.33 | 33.83 → 14.28 | 10.95 → 1.50 | 7.19 → 2.83 |
| 160 | 92.41 → 33.57 | 27.15 → 14.70 | 7.00 → 1.34 | 5.18 → 2.75 |
| 256 | 79.95 → 40.72 | 31.02 → 19.15 | 7.48 → 1.34 | 4.43 → 3.71 |

Route construction falls 75.36%, 63.67%, 49.06%. Summing only the disjoint named
preparation self regions 68–73 and 84 gives 216.94→62.74 ms,141.50→58.56 ms and
133.12→72.52 ms. The corresponding profiled whole windows fall 168.24, 89.23 and
57.22 ms. These reductions support the original structural-construction cost
case. They do not isolate the individual contributions of hash mixing, indexed
metadata pulls, communicator removal or allocator behavior; all are one repair.

RHS exchange self time is 136.62→129.16 ms,128.50→130.66 ms and 143.46→145.05 ms.
It barely changes at 160/256 despite the clear bare speedup. The result supports
repairing construction and local lookup costs rather than assuming message
batching or transport tuning is the cause. New indexed metadata pulls sit in
preparation parents, while prior generic pulls also contributed to exchange;
region 78 is not an identical pure-transport measure across sources. MPI waits
also include late producers and cannot establish transport causality.

Coarse perturbation versus bare mean-run-median is 1.11%, 2.50%, 4.91% for v7;
versus ordinary means it is 1.26%, 4.52%, 9.26%. Logged coarse ordinary means exceed
the profiled windows by 3.39, 9.56, 29.23 ms. Do not omit that difference, infer its
cause without evidence, or subtract profiled means from bare medians. Attribution
is directional cost evidence; bare runs establish the performance gain.

## Scaling remains unresolved

At 160→256, previous median time rises 598.31→608.27 ms (+1.66%); v7 rises
505.52→544.53 ms (**+7.72%**). Speedup changes 0.9836→0.9284 and relative
160→256 efficiency 61.48%→58.02%. Stage 266 itself gives 0.9292 speedup/58.08%
efficiency. All four 256-rank v7 medians exceed all four 160-rank medians.
Ordinary means agree:493.71→536.61 ms, speedup 0.9201 and efficiency 57.50%.

This is an absolute-time improvement at 256, but no scaling improvement. Savings
are 92.79 ms at 160 and 63.74 ms at 256; that 29.05 ms difference accounts for the
increase in the cross-count gap from 9.96 to 39.01 ms. Removing a larger low-count
cost makes relative scaling look worse without making either count slower than
before. It does not excuse the remaining scaling failure.

The v7 coarse profile window grows 506.48→557.08 ms (+50.60 ms). Its disjoint
RHS self regions 68–84 contribute+30.12 ms net; remaining regions+20.48 ms net.
Largest individual increases include inverse MPI wait+19.97 ms, RHS exchange
+14.39 ms, adaptation+12.38 ms, native wavelet+11.63 ms, migration MPI+8.13 ms,
catalog construction+7.51 ms and route construction+7.15 ms, partly offset by
other decreases. Do not sum those increases alone as a full-step budget or treat
MPI waits as network-only. Profile perturbation is larger at 256, so this single
profile pair cannot precisely partition the bare 39.01 ms median penalty.

The inherited inverse/adaptation/publication paths therefore remain relevant
alongside RHS execution and plan construction. J5>=640-rank efficient scaling
on matched multi-node hardware remains untested;256 is only intermediate.

## Memory and lifetimes

Mean summed per-rank peak RSS, GiB. These are sums of individual child wait 4 peaks,
not synchronized node peaks or memory-at-a-specific-step measurements.

| Ranks | Previous | V7 | 266 | V7 vs previous | V7 vs 266 |
|---|---:|---:|---:|---:|---:|
| 83 | 27.65 | 24.93 | 21.49 | -9.84% | +16.00% |
| 160 | 54.04 | 43.79 | 31.58 | -18.97% | +38.67% |
| 256 | 80.89 | 61.30 | 43.04 | -24.22% | +42.42% |

The repair saves 2.72/10.25/19.59 GiB by this metric, but v7 still exceeds 266 by
3.44/12.21/18.26 GiB. Largest-rank mean peaks also decline 457.50→411.61 MiB,
420.06→351.69 MiB and 494.10→377.31 MiB. The directory is preparation-local and
sixteen disposable contexts are removed by source construction; the campaign
shows lower RSS, not a decomposition of MPI allocation versus solver allocation
or proof against long-run growth. The numerical narrow ledger is unchanged and
cannot account for all this RSS improvement or the remaining excess.

## Decision and next work

1. Retain v7 as the measured development base for the next repair. The substantial
   low-count construction issue is now demonstrably reduced and 266 absolute
   parity surpassed. Do not roll back successful work merely because the final
   dual goals remain incomplete. No further campaign is needed to reconfirm this
   large same-fixture gain without new uncertainty or changes.
2. Before choosing the next implementation, use the existing source and rank
   evidence to identify the remaining inverse/lifecycle and RHS producer/reader
   dependencies, retained allocations and rank-count growth. At 160 the entire
   remaining profiled preparation subtotal is 58.56 ms, versus a140.40 ms bare
   main gap; even its elimination is not a defensible claim to close main parity.
   These are different statistics, so this is a scale comparison, not an additive
   speedup forecast. Another metadata-only optimization cannot be the whole plan.
3. Prioritize a concrete ownership/dependency or storage-lifetime change that
   addresses the inherited paths and the 256 penalty together, preserving direct
   fused kernels and paired Domain-path retirement. Specify what actual producer
   lateness, local work or redundant storage it removes; generic communication
   batching, unconditional cross-epoch reuse and MPI shared memory still lack an
   isolated causal case. Add only bounded instrumentation needed to resolve a
   named implementation decision, not another broad diagnostic campaign.
4. Keep absolute time, scaling, memory and whole-timestep independence separate
   on the scorecard. Stage 266 remains accepted pending the broader gates; v7 is
   a successful tested repair, not completion of independent whole timesteps,
   legacy parity or useful >=4 D execution. No source edits, new kit, upload or job
   submission were made during this review.

[Independent verification](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/review-55504/verified.json), [derived metrics](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/review-55504/metrics.json),
[integer audit comparison](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/review-55504/integer-equivalence.json),
[raw Slurm accounting](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/review-55504/review/accounting.json).
