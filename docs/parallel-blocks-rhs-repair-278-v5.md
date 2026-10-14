# Stage278 preparation and communication repair — kit-v5

Implemented after job 55501. Stage266 remains accepted; the repaired candidate passes local science and compiler qualification. Performance acceptance remains open locally and on the cluster; cluster memory and above-Domain qualification are also pending. No Slurm job was submitted for this repair.

## Why the smaller-task gap matters

The absolute ordinary-step excess over 266 in job 55501 falls from 363.067 ms at 83 ranks to 253.997 ms at 160 and 206.249 ms at 256. This is not sufficient evidence for a latency-only explanation. Against main, the relative penalty actually increases from 104.64% at 83 to 126.34% at 160. Between 160 and 256, candidate time improves 2.61% while 266 worsens 4.82%; some narrowing of the gap is the control losing scaling.

Existing raw validation timers give a more specific target. On ordinary steps, mean per-step rank maxima for backend preparation/RHS evaluation are 537.472/152.810 ms at 83 and 341.496/134.519 ms at 160. Preparation is rebuilt in all 12 timestep epochs. The validation ordinary medians exceed bare controls by 1.02%/4.79%; these unsynchronized region maxima are not an additive decomposition of the bare critical path. The 256 validation has an oracle and is 74.42% slower than bare, so its total timing cannot be treated as production.

Inspection found repeated global hashing and delivery of interior definitions/requests for every phase and field width; repeated boundary-map construction; a whole-plan deep copy before filtering; and repeated growing-array concatenation during placement. Larger local problems at smaller P amplify these memory/compilation costs. This is a source-supported mechanism, with the local measurements below testing the repair; it does not attribute every millisecond of the cluster gap.

For scale, the 83/160 gaps over 266 are about 68%/74% of the measured preparation region. Main gaps are larger than that entire region. Those ratios are diagnostic budgets, not a forecast: preparation repair is a plausible large contribution toward 266 parity, but it cannot establish main parity or good high-P scaling by itself.

## Changed design

- Canonical ownership is published once per executor preparation. Retained ghosts and physical boundary aliases bind to that directory once. Per-phase final-writer selection uses indexed local canonical slots; only nonlocal writer overrides and consumer requests pass through MPI metadata exchange. Accepted write ordering remains authoritative. The generic hash compiler remains an independent component-test reference.
- One route handles one/two field widths: 21→15 persistent phase routes per level. Runtime offsets use the actual width, with separate state stamps for every transaction. Duplicate numeric work or fields are not reintroduced.
- Build node/edge boundary maps once per level: 5→2 constructions. Build component key layouts once per preparation. For this route/map portion, explicit preparation collectives change from 312L to 214L+8 for L compiled refinement levels, excluding communicator creation/free, geometry distribution and other unchanged solver preparation. Runtime numerical message counts are unchanged by this repair.
- Placement copies selected records directly instead of deep-copying a full Domain program. Pressure/producer arrays use a count/fill allocation, preserving patch order and ordinals without repeated append/copy growth. Geometry packet packing and receiving also use counted final storage; received programs are unpacked directly into final entries. Merge arrays are allocated once and filled in original order with identical address/operand offsets. Packet/identity buffers are released before merge, and sender selection temporaries expire before transfer. These changes remove the accumulated-program deep copies and growing merge arrays; one temporary program copy per merge input remains.
- One borrowed context per executor remains; routes release before its communicator. The owner directory and home/alias tables exist only during preparation. No cross-epoch geometry cache, shared-memory window, partition change or additional numerical migration is introduced.

## Local qualification and quantitative results

Fresh checked and optimized builds; full 12-step checked, optimized and independent-oracle replays; exact accepted-family checkpoint/state/62-member VTK comparisons; two reloads and two remaps. The four-rank fixture exercises 40 Domains, so it does not replace above-Domain cluster qualification. All 23 external physics objects match the accepted build manifest.

The approved bbserv compiler check separately replays all 73 checked compile/link commands with GCC 13.2 and MVAPICH 4, using the final source manifest and unchanged external cluster physics objects. It does not run the solver or submit a campaign.

Component tests cover canonical-owner/hash-route equivalence, signed aliases, cycles, final-writer overrides, rank without canonical ownership, missing/ambiguous owners, variable widths, invalid epochs/tags/extents, repeated release and stress. Failed attempts remain archived: an unused test declaration and a local MPI sandbox socket rejection were corrected without weakening gates.

Optimized validation totals (maximum rank cumulative timer over 12 epochs, seconds; separate runs, not isolated causal attribution):

| Preparation stage | Backend preparation | RHS evaluation |
|---|---:|---:|
| Reviewed v4 | 18.798152 | 15.176867 |
| Width/map/placement repair | 12.211702 | 15.171398 |
| Owner binding, before linear geometry | 9.098407 | 15.281565 |
| Final v5 with linear geometry assembly | 8.415920 | 14.344197 |

The final preparation timer falls by 55.23%. Exact derived traffic remains 697,200 sends and 3,301,517,760 payload bytes. Input/RK traffic and numerical primitive counts also remain unchanged. Persistent counted route storage summed across ranks falls by 2.29%. This narrow ledger excludes compilation peaks, retained geometry, descriptors and MPI; it does not prove an RSS reduction.

Bare local J4 ordinary-step medians in seconds; negative percentage changes mean faster. At 4 ranks, numbers average two run medians in forward/reverse order. At 1/2/3, v4 and final-v5 have one fresh paired run each. A fresh 266 three-task control follows the final candidate to check the earlier parity exception. Main controls at 1/2/3 and 266 controls at 1/2 come from the preceding matrix and are not contemporaneous repetitions. All runs are serial on the same local machine, four usable cores, no MPI binding. Steps within a run are not independent repetitions.

| Tasks | Main | 266 | v4 | v5 | v5 vs v4 | v5 vs 266 | v5 vs main |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 5.677164 | 14.963110 | 12.143285 | 10.050743 | -17.23% | -32.83% | +77.04% |
| 2 | 3.840410 | 8.665677 | 10.181685 | 6.898400 | -32.25% | -20.39% | +79.63% |
| 3 | 3.144740 | 6.063612 | 11.457831 | 6.649054 | -41.97% | +9.66% | +111.43% |
| 4 | 2.669437 | 5.038540 | 5.126727 | 4.383616 | -14.49% | -13.00% | +64.22% |

The final local matrix still exceeds the 266 ordinary-step control at 3 tasks. The fresh control confirms a 9.66% ordinary-step deficit; the full-trajectory deficit is 22.45%. Local 266 parity is therefore not established uniformly. The single-run comparisons, external controls at 1/2 and unbound local scheduler limit causal interpretation; the exception is retained, not averaged away.

The earlier owner-binding matrix exposed a three-task penalty of 14.69% against 266. Its full 20 runs are preserved alongside the final 15 runs (the 14-run result is also preserved before adding the fresh three-task 266 control). The geometry change removes growing-array/deep-copy work by construction; separate whole-run samples do not establish an isolated causal speedup for that subchange.

Complete 12-step trajectory time, summing ordinary and exceptional solver steps (seconds). These include adaptation/remap steps; launcher elapsed times, including startup and output, are recorded separately in each run. The same external-control and repetition limits apply.

| Tasks | Main | 266 | v4 | v5 | v5 vs 266 |
|---:|---:|---:|---:|---:|---:|
| 1 | 86.163 | 202.206 | 166.417 | 138.768 | -31.37% |
| 2 | 53.221 | 111.877 | 137.450 | 93.778 | -16.18% |
| 3 | 42.397 | 77.761 | 165.030 | 95.216 | +22.45% |
| 4 | 36.258 | 66.599 | 70.906 | 59.226 | -11.07% |

This is a local architecture/performance screen, not a prediction for J5 bb13. Main uses the established numerical tolerance; 266/v4/v5 retain exact same-rank science and VTK gates. Per-run ordinary, complete-trajectory and launch times are in the evidence.

## Dual goals and next acceptance decision

Connected block-owned RHS and real RK consumption remain active. Displaced Domain RHS arithmetic, scalar capture/transport, velocity remainder export, full temperature-stage publication and interior Domain RK stay retired. Whole-timestep independence remains incomplete: boundary/scaffold and inverse/lifecycle adapters, split Simple physics and retained geometry/workspaces remain named work.

The next matched bb13 campaign must test this source against 266 and main at 83/160, and 266 at 256, with the existing checked/oracle/zero-Domain science gates. Include bare ordinary and complete-trajectory time, rank RSS, replacement costs, topology and placement. Compare the returned source and controls to 55501; do not promote based only on a better candidate-relative result. Required v4 cuts were 32.70%/31.46%/26.23% to 266 and 51.13%/55.82% to main.

The high-P floor is still open: v4 sends 87,006,360 derived messages at 256 with roughly 395-byte mean payload. This repair reduces preparation rounds; it does not remove those runtime exchanges. If preparation falls but full-step parity is still missed, the next targeted repair is producer-ready batching across compatible layers/phases and removal of unnecessary intermediate refreshes while preserving hydrostatic order and signed alias snapshots. Do not restore Domain arithmetic to manufacture a speedup.

Overall memory remains unqualified until matched RSS measurements; a bounded preparation directory is not a proof that the 2.73× v4/266 memory ratio disappears. MPI_Comm_split_type/shared windows remain a separate option requiring evidence about local transport cost and NUMA placement; they do not remove the unnecessary work repaired here or establish inter-node scaling.

J5 requires at least 640 ranks for the fourfold goal. 256 remains intermediate; subdivision and balance need separate qualification. Stage266 remains the accepted baseline until the full scorecard permits promotion.

## Evidence

- [GCC 13/MVAPICH compiler check](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/repair-v5/compiler-check/results.json)
- [Source-bound qualification](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/repair-v5/final/qualification.json)
- [Bare matrix](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/repair-v5/final-bare/results.json)
- [Cost comparison](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/repair-v5/final/cost-comparison.json)
- [55501 timer derivation](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/repair-v5/rk278-gap-analysis.json)
- [Changes relative to v4](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/repair-v5/final/solver-changes.patch)
- [Previous campaign review](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/docs/parallel-blocks-stage278-job55501-results.md)

## Delivery

The sealed kit is on the local Desktop; it has not been transferred to bbserv or submitted. All 12 campaign helper tests, packaged imports, shell syntax checks and fresh archive extraction passed. The authorized bbserv activity was the separate source-overlay compiler check.

- [Sealed kit](/Users/kevlahan/Desktop/rk278-connected-rhs-kit-v5.tar.gz)
- [Delivery identity](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/repair-v5/delivery.json)

## Production timing clarification (2026-10-13)

The user prioritizes ordinary complete-step performance: production checkpoint/data
saves occur only every O(10^3) steps, while remaps occur about every five steps.
Use the established ordinary metric as primary; it excludes save/startup/reload
steps and retains ordinary remaps. Separate remap/non-remap costs in the returned
analysis. Treat the short fixture whole-trajectory totals above as secondary
lifecycle/overhead evidence rather than production throughput. This updates the
review priority without changing the sealed kit or its source and fixture hashes.
