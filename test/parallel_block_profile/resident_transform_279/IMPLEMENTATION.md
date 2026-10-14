Current delivery: v3, validation-only correction. See [checker fix](VALIDATION_FIX.md).
Solver source remains exact v2.

Current delivery is v2. See [GCC13 correction](PORTABILITY_FIX.md).
The implementation and timing results below describe v1; current source
qualification is in LOCAL_VALIDATION.json.

# Resident provisional transform repair — implementation and qualification

2026-10-14. Development base: job55504/v7, commit
`de558ed03798f721991f38eeecc00f00eafc8370`. Stage266 remains accepted.
The repair is independently selectable with `WAVETRISK_RESIDENT_TRANSFORM=1`.
Default behavior remains v7. The campaign sets this only for candidate/checked.

## Connected ownership and paired retirement

The three provisional RK4 transitions now execute resident canonical input →
scalar/vector forward transform → compression → inverse → real next-RHS input
on numerical owners. The production branch bypasses `WT_after_step` completely
for these transitions. It also asserts that the Domain production-writeback
counter did not change. The final transform and its adaptation/remap/physics
consumers retain the qualified path.

| Displaced provisional work | Replacement |
|---|---|
| Domain-indexed vector forward staging and coefficient import | Fixed compact vector kernels on numerical owners; scalar forward runs there too |
| Compressed coefficient publication to Domain interiors | Resident coefficients feed the inverse directly |
| Geometry-owned scalar/outer/inner inverse arithmetic | Fixed all-column kernels placed with their actual outputs; in-place outer dependencies grouped during compilation |
| Initial canonical inverse gathers and final interior installation | Prebound local canonical reads and local result commits |
| Canonical Domain solution publication | No provisional canonical mirror; explicitly retained uncatalogued compatibility and block-boundary views are updated |

No credit is claimed for already removed in-phase Stage203 gathers, Stage266's
deferred first installation, or full-SOL network returns (already zero).
The existing raw-RK compatibility producer and fixed timestep-sol scaffold inputs
remain named dependencies. Canonical values, dependency aliases and consumer
views preserve their distinct versions. Whole-timestep independence is unfinished.

The independent oracle first computes the candidate result without modifying
production storage, then executes the qualified reference path, compares the
forward footprint/coefficient results and final canonical/Domain/block-boundary
values bit for bit, and restores candidate SOL for actual subsequent consumers.
This duplicate work is validation-only. Separate optimized/checked runs establish
the production cutover and retirement; oracle timings never establish speedup.

## Preparation, locality and lifetime

Kernels retain direct all-column arithmetic and original operation order; the
runtime has no owner lookup, Domain geometry lookup or dependency-depth interpreter.
Physical ownership resolves once per preparation and feeds direct final-writer
route binding. Plans invalidate on topology/ownership and RHS mask generations;
matching counts are never accepted as an epoch identity.

Canonical inputs use existing local block addresses. A first implementation
unnecessarily compiled a general resident-provider exchange and retained its
extra staging arrays. Those plans and buffers have been removed.

All transform routes borrow one communicator per backend epoch. Every operation
waits for all requests before returning. The context is freed only after all
borrowing routes are released, and the backend is released at restart clearing.
One numerical exchange scratch allocation serves sequential routes and is sized
by their maximum requirement; the initial implementation allocated their sum.
The separately bound raw-RK boundary adapter uses bounded requested-reader buffers.
There is no new shared-memory window or per-field/per-block communicator.

The first production version's local bare medians were5.217/5.700s versus
4.837/4.650s for v7. It was not sealed. Direct local reads and pooled exchange
scratch address concrete redundant construction, copies and storage. Summed
per-rank maxima of named transform numerical/route allocations fall from
1450.339MiB to601.553MiB (58.52%) on the same four-rank J4 trajectory.
This compares two implementations of the new transform, not total memory versus
v7. It excludes kernel/geometry metadata, descriptors, construction peaks, the
retained final-transform workspace and MPI. Whole-process RSS remains a separate
cluster gate; old/new final-transform workspace overlap has not disappeared.

Preparation timers distinguish geometry compilation, placement, boundary metadata,
fixed inputs, home resolution and route binding. Coarse regions85–90 distinguish
new plan/input/forward/inverse/output/exchange time. These are attribution aids,
not claims that MPI transport dominates. Exchange includes packing and late
producer waits; route invocations are not global synchronization-round counts.

## Scope and remaining architectural work

The original J5/D160 workload remains decisive. The campaign checks83/160 absolute
time against main and266, and candidate160→256 scaling independently. The fixed
640-block catalogue and first-block scaffold assignment remain unchanged. Per-
block compiled transform operation vectors are captured alongside existing RHS
work to inform a separate granularity repair. These descriptor counts include
inactive/structural operations and are not calibrated FLOPs or time weights.

The measured v7 main gaps remain92.09ms at83 and140.40ms at160 until a new cluster
result demonstrates otherwise. The broad old transform envelopes208.11/142.25/
152.56ms include necessary work and callers outside this slice. They are not
removable-cost or speedup forecasts; this provisional slice cannot defensibly
promise full legacy parity. Efficient J5 >=640 execution, finer ownership, final
lifecycle closure, long production workloads and other solver configurations
remain required work. J6/J8 and explicit node-local shared memory do not replace
those goals.

## Campaign

The sealed kit, once its entry gate passes, contains exact unchanged main,266
and v7 source controls plus this candidate. Five serial builds and49 serial runs
on exclusive nk/bb13: nine validations,34 bare measurements and six coarse runs.
Four v7/candidate repetitions per count use ABBABAAB order; main/266 controls
bracket them. An audited production256 run is separate from oracle256.

Science gates retain12 states,2 reloads,2 remaps and62 VTK members, exact accepted-
family results and the established main-relative tolerance. Ownership/cost checks
require36 provisional transactions per rank, bounded context/epoch lifetimes,
balanced traffic and useful work on zero-Domain ranks. Ordinary timing excludes
startup/reload/saves at indices0,7,8 and retains remaps4,9. Saves everyO(1000)
production steps remain secondary. Review run medians, means, spread, remap groups,
absolute legacy gaps, rank RSS and scaling separately; no automatic promotion.

Final source-bound test results and the local timing table are recorded in the
qualification report delivered with the kit. Local timings are not cluster
predictions: the unchanged control's coarse median varied4.66→6.90s while the
machine was heavily loaded. The exclusive campaign is the decisive cost test.

## Final local qualification

All nine full four-rank J4/J6 replays pass, including optimized production,
checked production and independent oracle. Each retains12 states,2 reloads,
2 remaps and62 VTK members. The portable MPI harness passes24 cases across
optimized/checked, default/shared buffers, aliasing, empty ranks, stale epochs,
insufficient scratch and large rendezvous payloads;32 campaign-helper tests pass.
Local J4 has40 Domains, so above-Domain production validation remains a cluster gate.

| Bare run | Variant | Ordinary median (s) | Ordinary mean (s) |
|---|---|---:|---:|
| A1 | v7 | 7.240100 | 7.075902 |
| B1 | candidate | 7.218158 | 7.389441 |
| B2 | candidate | 6.673557 | 6.991783 |
| A2 | v7 | 6.364683 | 6.682450 |

The mean of run medians is6.802392s for v7 and6.945858s for candidate
(+2.11%); ordinary means are6.879176s and7.190612s (+4.53%). These point
estimates are slower and the run ranges overlap. Substantial uncontrolled host
load prevents a defensible speedup claim. The candidate coarse RHS exchange
region also grows; producer lateness introduced by the new pipeline could
contribute, so this must not be dismissed as unrelated noise or pure transport.
The bounded replacement cost case and scientific qualification permit the
exclusive cluster comparison; they do not constitute performance acceptance.
Do not stack another cutover or promote this candidate before reviewing it.

The sealed49-run kit is prepared for user submission. No remote upload or
submission has occurred. Source hashes bind the candidate, controls and tests;
full raw local evidence is retained under `.stage279-validation/local-evidence`.
