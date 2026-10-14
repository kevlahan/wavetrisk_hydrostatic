Current reviewed result: [kit-v4 job55501](parallel-blocks-stage278-job55501-results.md) passes science/ownership and improves traffic/time, but still fails absolute performance, memory and scaling acceptance. Stage266 remains accepted.

# Stage278 job 55500: ownership passes, performance fails

Reviewed 2026-10-13. **Do not promote this candidate. Stage266 remains accepted.**
The connected RHS/RK migration is scientifically validated on the cluster,
including actual production work on every zero-Domain rank. Its complete-step
performance fails the other goal: ordinary steps are 46–56% slower than 266, and
160-to-256 scaling is essentially flat. Preserve the implementation and its
ownership evidence while repairing the cost regression; do not restore duplicated
Domain arithmetic or count this as a completed efficient timestep migration.

## Verified campaign and source

Archive SHA256 `3be6b82863ee0437ea9fc59cb676488af9bb2bf35758ada34ba328f5820e2220`; **4,874 members verified**.
Returned kit identity, source manifest, inventory and entry gate exactly match
sealed kit-v3 (`e72cd9e60334c709e249f3569dbe20ecd8b2af436cde06c2fdb3a1fdcfbb1091`).
The current 23-file candidate overlay also matches the tested source. The reviewer
matches the sealed helper byte for byte. Its eight-ulp allowance applies only to
derived floating-point report sums across Python versions, not raw hashes,
scientific tolerances, integer counters or source identities.

All four serial builds and all sixteen solver runs pass: eight validations and
eight bare timings. The GCC13 checked-build failure from job 55499 is resolved.
Slurm records **11:28:31–11:54:09 on October 13, elapsed 25:38**, on nk/bb13.
Runner UTC is 15:28:31.603825–15:54:09.310219, elapsed 1537.706436s. Queue and
collection are excluded. Observed execution lies within the 25–40-minute allowance.

Every run completes12 steps, two reloads and two remaps with the required
scientific gates. Accepted-family checkpoint comparisons are exact; main-relative
comparisons satisfy the established tolerance. Each run records 62 VTK members,
with same-rank/family hashes checked. Optimized, checked and independent-oracle
coverage is retained. Bare timings contain no ownership audit or duplicate oracle.
Binary, checkpoint and VTK payloads are not in the returned archive: this review
recomputes from raw logs/counters and verifies cluster-generated hashes/reports;
it does not rerun the solver or independently recompare omitted payload bytes.

## Dual-goal result

**Ownership:** all 96 zero-Domain ranks at 256 execute nonzero authoritative work:
16,212,000 primitive events in aggregate, minimum 23,040 per rank, with 48 real
RHS/RK retention transactions on each rank. This advances beyond Stage277, whose
zero-Domain ranks performed no primitive work. The selected path retires Domain
compatibility RHS, scalar tendency capture/transport, velocity remainder export,
full temperature-stage publication and displaced interior mass/velocity RK.
Runtime counters and the guarded source paths support this qualified slice.

Whole-timestep independence remains unfinished: boundary/scaffold and named RK
reader adapters, inverse/lifecycle canonical publication, separate split Simple
physics and old geometry compilers/workspaces remain. This candidate is still
unpromoted. The qualified scope remains restarted compressible non-mode-split
RK4 climate with Simple physics, no NCAR/SSO, and the delivered diffusion fixture.

**Performance:** ordinary-step medians below exclude the same steps0,7,8. One bare
run per configuration supplies no confidence interval or independent repetitions;
the large consistent penalties suffice to withhold promotion and prioritize repair.

| Ranks | Main ms | Stage266 ms | Candidate278 ms | 278 vs266 | 278 vs main |
|---:|---:|---:|---:|---:|---:|
| 83 | 546.195 | 748.756 | 1138.931 | +52.11% | +108.52% |
| 160 | 354.309 | 550.466 | 861.217 | +56.45% | +143.07% |
| 256 | unsupported | 586.341 | 858.293 | +46.38% | — |

The matched266 gaps to main are 37.09% at 83 and55.36% at 160; the candidate widens
them to 108.52% and143.07%. Candidate runtime must fall 34.26%,36.08%,31.69% at
83/160/256 merely to match 266, and 52.04%/58.86% at 83/160 to match main.
These are measured break-even requirements, not predicted gains.

Complete12-step time is also worse:266→278 is 14.237→18.944s (+33.07%) at 83,
12.026→15.406s (+28.11%) at 160, and13.202→17.085s (+29.41%) at 256.
Full launcher elapsed is 29.001→33.624s,40.081→43.482s,63.970→68.338s respectively.
The regression is therefore present beyond the ordinary-step statistic.

Candidate 160→256 ordinary-step speedup is **1.0034×**, relative efficiency 62.71%:
60% more ranks provide only 0.34% speedup. The complete trajectory gets 10.90% slower.
Stage266 speedup is0.9388×, relative efficiency 58.68%. A slightly better efficiency
ratio from a much slower candidate is not useful performance acceptance.
256 remains an intermediate test for D=160; useful>=640 execution is unqualified.

## Resource and replacement evidence

| Ranks | 266 max-rank MiB | 278 max-rank MiB | 266 sum-peaks GiB | 278 sum-peaks GiB | Sum ratio |
|---:|---:|---:|---:|---:|---:|
| 83 | 373.8 | 489.7 | 21.46 | 30.35 | 1.41× |
| 160 | 266.1 | 472.7 | 31.57 | 63.27 | 2.00× |
| 256 | 269.5 | 598.4 | 43.07 | 118.76 | 2.76× |

These are solver-child wait4 peaks. Summing independently timed rank peaks does
not give a simultaneous allocation peak. Major faults include executable/file
faults and do not alone establish swapping or memory exhaustion.

Whole-trajectory validation counters (input + derived + RK stage streams):

| Ranks | Payload GiB | Send messages | Max-rank preparation s | RHS evaluation s | RK adapter s |
|---:|---:|---:|---:|---:|---:|
| 83 | 37.22 | 28,857,168 | 6.451 | 2.210 | 0.261 |
| 160 | 50.78 | 62,216,196 | 4.122 | 2.024 | 0.220 |
| 256 | 61.65 | 106,115,536 | 4.130 | 2.211 | 1.738 |

Each message carries an additional 8-byte epoch header. Geometry/placement metadata
is additional; the named placement stream sends 0.503/0.588/0.602GiB at 83/160/256.
The 256 validation uses independent oracles, unlike83/160; its timers are not a
matched bare performance decomposition. Counts describe the selected connected
routes, not all solver MPI. Rank maxima may come from different ranks and must
not be added or subtracted from ordinary-step medians to forecast savings.

The counted component memory upper bound is only 26.7 MiB/rank at 256, while observed
candidate max RSS is 598.4 MiB. The component ledger explicitly omits compiler
transients, retained geometry/workspaces, descriptors and MPI internals, so it
cannot attribute the RSS increase by itself. Source inspection finds a separate
`MPI_Comm_dup` for each compiled phase route, across levels and variable channels;
this is a concrete resource/preparation suspect, not a measured causal verdict.

Load distribution also remains uneven: the 96 zero-Domain ranks each have one
block and perform 22.82% of primitive events despite comprising 37.5% of ranks.
Primitive-event share is not a cost-weighted balance measure. Keep subdivision/
partition experiments separate from numerical ownership and communication fixes.

## Next decision and bounded work

Keep 266 accepted; retain 278's science, ownership and retirement proof. Do not
advance the unchanged candidate to a multi-node/640-rank performance campaign.
A repeat of this exact broad campaign is not needed merely to distinguish these
large penalties from a small mixed result. No statistical repeatability claim is
made; a focused repeat is justified only if an identified confound changes the
implementation decision.

The next substantial change should reduce the connected replacement's phase
communication and preparation/resource costs while preserving resident ownership,
fused kernels, ordering and exact science. Evaluate shared communicator lifetime
with collision-safe phase tags and the actual remote dependency footprint first.
The 106.1-million-message256 trajectory and per-route communicator construction
justify those targets; they do not promise that either alone closes the 272 ms
ordinary-step gap. Avoid another per-value dependency interpreter, a full primitive
return, or a mixed partition/numerical cutover. If attribution is needed, bound it
to communicator/route storage and compilation peaks to choose the implementation.

Gate the revised connected candidate locally in optimized/checked/oracle modes,
then repeat the matched 83/160/256 bb13 comparisons with all preparation and adapters
charged. Required recovery is the break-even reduction above, followed by closing
the remaining main gap. Only then consider the separate matched bb readiness and
160/320/640 milestones; no numerical efficiency threshold has been agreed.

Evidence: [independent review](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/review-55500/verified.json),
[derived metrics](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/review-55500/metrics.json), [original archive](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/review-55500/rk278-review-55500.tar.gz),
[scorecard](/Users/kevlahan/.codex/worktrees/26d8/wavetrisk_hydrostatic/docs/parallel-blocks-scorecard.md).


## Design reassessment: repair scaling before completing independence

2026-10-13 follow-up. **Repair communication, memory and scaling now.** Whole-timestep
independence may eventually eliminate compatibility materialization, duplicate
representations and Domain-owner bottlenecks. It has not been demonstrated to make
this distributed schedule efficient, nor is it necessary to change the schedule.
A narrowly scoped ownership move can accompany a repair when it eliminates a
specific measured crossing; broad migration first is not justified by this evidence.

The new traffic decomposition is decisive about scope, though not about elapsed
time attribution. At256, derived-field routes carry105,795,480 of106,115,536 recorded
stage messages (99.70%) and52.91 of61.65GiB payload (85.82%). Input and RK routes
together account for0.30% of messages and14.18% of payload. Even deleting all those
input/RK messages—an unrealistic upper bound on their removability—leaves the
large derived-field stream. Remaining inverse/physics/legacy MPI is outside this
ledger; message shares are not time shares or speedup forecasts.

| Ranks | Derived-field messages | Share of counted stage messages | Mean payload bytes |
|---:|---:|---:|---:|
| 83 | 28,751,280 | 99.63% | 1129 |
| 160 | 61,993,200 | 99.64% | 735 |
| 256 | 105,795,480 | 99.70% | 537 |

Concrete design barriers in the tested source:

- `parallel_block_rhs_distributed.f90:124` runs the complete phase sequence for
  each of30 physical layers. Each phase completes `MPI_Waitall` before returning
  (`parallel_block_rhs_exchange.f90:305`); nonblocking calls therefore do not
  provide useful computation overlap in this routine. This exposes many small
  exchanges and dependency waits as work per rank shrinks. It is not a global
  barrier at each phase, and not every dependency can be overlapped or merged.
- `parallel_block_rhs_routes.f90:32` requests every slot/component in the retained
  program footprint for a phase. It does not restrict that phase's requests to
  values actually read before their next overwrite. Later identity pruning and
  producer-slot deduplication help, but do not provide phase-specific liveness.
  This is a correctness-conservative distributed workspace refresh, not a proof
  of a minimal numerical halo. Exact aliases, signs and writer order must survive
  any pruning; the recoverable volume is not yet measured.
- `parallel_block_rhs_exchange.f90:126` duplicates the communicator for each route,
  performs global metadata collectives and allocates rank-sized route arrays.
  Plans rebuild across12 recorded epochs in12 steps. Communicator reuse alone
  does not reduce the runtime message count. MPI internal memory and compilation
  peaks remain unmeasured contributors to the RSS increase, not proven causes.
- The96 zero-Domain ranks supply37.5% of rank capacity but22.82% of primitive events,
  each with one block. This repeats the imbalance seen in Stage272. Event counts
  are not cost weights; the640-block ceiling also leaves at most one block per
  rank at640. Finer subdivision and runtime-aware assignment remain separate
  qualification work, not an automatic cure for message overhead.

The earlier Stage251–255 failure established that optimizing a costly execution
structure piecemeal can leave most of the regression. Stage278 uses direct kernels
rather than the old per-value owner interpreter, but its coarse numerical kernels
are surrounded by a fine-grained distributed schedule. Correctness and ownership
alone do not show that schedule will scale. The source provides evidence of
structural obstacles; it does not establish that block ownership itself is doomed.

Revise the immediate priority from communicator reuse alone to **phase-specific
live dependency sets and fewer communication rounds**, alongside bounded resource
lifetime. Evaluate field/phase aggregation and layer batching only where dependency
and exact ordering permit: hydrostatic column recurrence means a30-layer batching
speedup cannot be assumed. Maintain per-column order and measure any additional
workspace. Use focused attribution for unaccounted MPI/compiler/geometry memory.
Then test load balance/granularity separately. Preserve all qualified ownership
and scientific gates throughout; do not reinstate the rejected interpreter or
whole-state primitive returns.

Recovering266 parity is necessary but insufficient, since266 itself scales poorly
above D. Require the revised candidate's matched complete-step improvement,
controlled memory, legacy-gap closure and separately qualified useful P>D scaling.
Use a bounded inter-node check after the local execution design recovers; do not
wait for all timestep ownership work to finish before testing the relevant MPI
architecture. No final efficiency threshold is invented here.

[Source-bound stream decomposition](../.stage278-validation/design-reassessment-278.json)
and [historical failure evidence](/Users/kevlahan/.codex/worktrees/26d8/wavetrisk_hydrostatic/docs/parallel-blocks-measured-restart.md:376).
