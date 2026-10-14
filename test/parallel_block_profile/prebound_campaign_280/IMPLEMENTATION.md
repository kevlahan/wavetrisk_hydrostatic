# Stage280 — direct transform route binding and compact output storage

Implemented2026-10-14 on the existing repair worktree. The exact Stage279/v3
source from job55507 is the base; only parallel_block_rhs_exchange.f90 and
parallel_block_transform_backend.f90 change. Stage266 remains accepted; v7 and
Stage279 remain exact controls. This is a repair of the same provisional ownership
slice, without another numerical migration, partition change or unsafe epoch reuse.

## Evidence and implementation contract

Job55507 demonstrated clear83/160 regressions: +7.76%/+12.75% versus v7,
with new transform preparation128.83/117.23ms in separate coarse ordinary means.
Audited home resolution and route/adapter binding dominate preparation. Whole-
process summed rank peak RSS also grows4.71%/5.97%/6.35% at83/160/256.

The epoch's canonical ownership directory already proves exactly one owner and
workspace address per physical node. Ordinary scalar/vector, coefficient, identity
and publication phases declare only those canonical producers; they cannot have
competing definitions. Stage279 nevertheless ran final-writer discovery for every
phase. Stage280 passes the already-validated producer address directly to the
exchange compiler, with zero competing definitions. The compiler rejects a mixed
provider/definition contract, invalid ranks/addresses, inconsistent distributed
binding modes and incoming sources outside the sender's workspace.

For the qualified J4/J6 and J5/J7 ranges,16 routes take this direct path. Each
avoids three metadata deliveries(two collectives each) and two writer-validation
reductions:128 redundant preparation collectives per topology/mask build.
Both outer-vector routes retain their full ordered writer proof, including
noncanonical in-place aliases. Both fixed-input routes also retain the prior
compiler. Plan signature consensus and source-request exchange remain; no claim
that all preparation synchronization disappears. The numerical exchange engine,
packing snapshot, signs, duplicate-demand handling and peer state stamps are unchanged.

Physical-node lookup now mixes both key axes before table reduction, replacing
the unmixed linear key polynomial. The mapping is still compiled once for the
current epoch. No isolated speedup is attributed to this hash change.

Output views now bind to their required value: a local output reads its canonical
workspace row; remote outputs share one import row per unique physical node.
Only unique remote imports become publication exchange destinations. This retires
local publication snapshots and duplicate remote output rows, including unnecessary
send/receive scratch demand. Scalar/vector coefficient and input arrays remain
separate because their simultaneous versions are required. Existing raw-RK
boundary/fixed-scaffold adapters, alias saves and final-transform workspace remain.

All route buffers remain bounded; sequential routes share maximum-sized scratch.
One communicator per epoch is borrowed by the routes and freed after release.
Each runtime exchange drains its requests. The physical directory and temporary
compiler bindings expire after preparation. Topology, ownership and mask epochs
are still checked; all12 trajectory states still rebuild and reuse their routes
across three provisional transitions. No count-only cross-step cache was introduced.

## Source-bound qualification

Full local checked, optimized and independent-oracle trajectories pass. All36
provisional transitions preserve exact forward inputs/coefficient results, inverse
outputs and real subsequent RHS inputs. All11 full solver replays pass12 states,
2 reloads,2 remaps and62 VTK members, with exact accepted-family science. The
old Domain provisional paths stay retired; global compiled kernel censuses match
Stage279. The current source is bound by255 path hashes and a two-file delta.

48 serial MPI cases pass under optimized and checked builds at1/4 ranks. They
cover original competing writers, prebound providers, signed aliases, cycles,
identity pruning, empty/all-empty ranks, large rendezvous traffic, stale masks/
states, duplicate phases, insufficient scratch, invalid provider rank/slot,
mixed definitions and inconsistent distributed binding mode. Complete serial
GCC13.2/MVAPICH J5 checked and optimized builds pass in an isolated bbserv copy.
Only the two corrected source files were uploaded for that compiler check; no
solver or Slurm campaign was run remotely. All compiler/runtime checks remain on.

## Quantitative local replacement cost

Same J4/J6 four-rank trajectory; this is below40 Domains and cannot qualify256
or>=640-task execution. Matched separate bare processes use A/B/C/C/B/A order
(v7/Stage279/repair). Two repetitions per kind; nine ordinary steps exclude0,7,8,
retain remaps4,9 and are not independent repetitions.

| Variant | Run1 median (s) | Run2 median (s) | Mean run median (s) | Ordinary mean (s) |
|---|---:|---:|---:|---:|
| v7 | 3.867520 | 3.866214 | 3.866867 | 3.960152 |
| stage279 | 3.841657 | 3.939667 | 3.890662 | 3.995518 |
| repair | 3.739361 | 3.623598 | 3.681480 | 3.740309 |

The repair's mean run median is5.38% below
Stage279 and4.79% below v7. Ordinary means improve
6.39%/5.55% respectively.
Both repair medians are below all four control medians. This is useful local
performance evidence, not a cluster forecast or proof of legacy parity.

Summed rank maxima of named numerical and route allocations fall
601.553→423.263MiB
(29.64%). Numerical buffers alone fall
580.343→409.443MiB; route payload
21.210→13.821MiB. This is extra-transform
allocation versus Stage279, not total process RSS or memory parity versus v7.
Kernel/geometry metadata, construction peaks, retained final workspace and MPI
are outside that named count; cluster RSS remains required.

Remote numerical payload is exactly unchanged:294,320,610 sent and received
values across the validation trajectory, and3600 total rank-route invocations.
The mechanism removes preparation and local staging work; it does not establish
transport dominance or fewer numerical communication rounds.

One separate local coarse pair shows preparation165.24→157.13ms and exchange
212.19→178.82ms, while RHS exchange increases513.54→585.74ms. Whole profiled
ordinary mean changes3.81676→3.80259s. The smaller coarse benefit and shifted
RHS waiting prevent attributing the entire bare gain to an individual mechanism.
Do not dismiss that waiting growth as unrelated or add these self-times to bare
median differences. Single coarse profiles are diagnostic, not replicated speed tests.

## Scorecard and remaining gate

Provisional ownership/science/retirement are preserved. Quantified replacement
storage is lower and matched local absolute performance improves. This supports
a bounded next cluster comparison; it does not promote the candidate. There is
still no evidence of new83/160/256 gains, improved256 scaling, main parity,
whole-timestep independence, total-memory parity or efficient>=640 execution.
Root-child granularity and first-block scaffold placement are unchanged.

The measured job55507 v7-regression targets remain49.20/64.37ms at83/160 and
main gaps143.07/210.14ms until new cluster results demonstrate otherwise. Removing
128 proven redundant preparation collectives and substantial local output storage
is a concrete architecture change, not a promised latency-derived speedup. Do
not extrapolate four-rank gains or assume the entire117ms cluster plan budget can
be removed. Next qualification must retain exact v7, Stage279,266 and main
controls, above-Domain production/oracle coverage, bare ordinary-time repetitions,
RSS and separate scaling evaluation. No cluster kit has been submitted.
