# Proposed Stage 191: block lifetime and communication scheduling

Status: components A and C are implemented experimentally; B remains proposed.
Cluster job 54561 validates a 5.96% ordinary-step improvement attributed to A;
C alone shows no measurable benefit. See `parallel-blocks-cluster-191-54561.md`. Stage 190 remains
the validated numerical and performance baseline, commit
`cdd2eb64` (the following `694e251a` adds investigation only).
“Legacy” means the pinned Domain implementation
`b82647b36eb825aed0c0aae4e92adcf4ab46d608` with the common checkpoint I/O fix.

## Objective and scope

Reduce costs introduced by maintaining Domain and block representations together.
Campaign 54553 measured a 2.56x ordinary-step ratio to legacy. Block-only profiles
show approximately 0.178 s/ordinary step in construction/migration, including
0.044 s in serialization/copy proofs, and 0.078 s in native mass MPI wait.
These are measured block costs, not guaranteed recoverable savings or a causal
partition of the difference from legacy.

The original proposal has three components; the current candidate implements A and C only:

| Component | Concrete change | Success evidence |
| --- | --- | --- |
| A: production installation | Keep structural checks; gate round-trip proofs; transfer owned arrays instead of deep-copying temporary blocks | Proof bytes vanish in production; moved/copied bytes recorded; lower installation time and/or peak RSS |
| B: block lifetime | Reuse storage for surviving subtrees and retain valid plans; rebuild changed dependencies only | Fewer allocations and plan builds despite adaptive grid changes; lower construction/setup cost |
| C: mass communication | Start fine-level tendency exchange before velocity-source work; finish at its first coarse consumer | Same results and payloads; less exposed mass wait and lower ordinary-step time |

A alone is a useful preparation and control experiment, not the substantial
improvement sought. B and C are the main implementation targets. Avoid another
restriction-arithmetic-only stage while these structural costs remain.

## A. Production installation without redundant copy proofs

Introduce an explicit reconstruction-validation option propagated through
`build_source_blocks`, `migrate_blocks` and `install_local_blocks`. Full checking
runs retain pack/unpack/pack and byte-equivalence tests. Production keeps
mandatory allocation/extent, identity, ownership, route and MPI transfer checks
without the redundant serialization proofs.

Do not interpret `full_validation=.false.` as sufficient today. Source proofs,
received serialization verification, and local copy proofs are executed outside
that guard. Gate all three consistently, with rank agreement checked once at
configuration initialization rather than a new collective per timestep.

For production installation, transfer allocatable payload ownership from
retained-source/received blocks into the installed store where lifetime permits.
Preserve the copy-based path for explicit comparison. Move all owned components
as one audited operation, including geometry, fields, boundaries, ghosts and
auxiliary state; a partial transfer is not acceptable.

Source-lifetime correction after implementation inspection: the post-installation
source-patch-count reader is in full ghost-source validation, not in
`build_block_writeback_plan` as this proposal originally stated. Full validation
therefore retains the copy path and source arrays. Production can consume the
source allocations after routing/install metadata has been established.

Measure actual moved bytes, copied bytes, proof bytes and peak staging storage.
Temporary source/receive staging is already released at the end of migration;
this change reduces construction copying and transient peaks, and must not be
advertised as eliminating a permanent leak or halving the full memory footprint.

## B. Retain block storage and plans through adaptation

The current path reconstructs all blocks after Domain-based adaptation. Replace
that whole-store replacement with an explicit lifecycle for surviving, changed,
new, deleted and migrated subtrees.

Use a stable geometric subtree key, such as root Domain plus dyadic child path.
Current local patch indices and catalog IDs are not stable across compaction or
redistribution and cannot be reused as persistent cache identities.

Separate validity of:

- Storage capacity and field layout.
- Geometry/connectivity and stencil addresses.
- Neighbor routes and final ownership.
- Current field values, masks and boundary publication state.

Start with storage-capacity reuse, which can preserve allocations even when a
plan must be rebuilt. Then retain compiled stencils/routes only when their full
dependency signature is unchanged. A neighbor's topology or ownership change
can invalidate a surviving block's halo plan. Mask changes may require rebuilding
an execution list even if its geometry remains unchanged.

Every surviving block must receive fresh fields after physics/adaptation.
Reusing geometry is not permission to reuse old mass, temperature, wavelets,
physics residuals, boundary values or cached divergence. Restart and ownership
changes must invalidate the appropriate state explicitly.

A global “skip rebuilding when the grid is unchanged” shortcut is insufficient:
the accepted fixture changes its adaptive grid almost every step. Reuse must
work for unchanged subtrees within a changing grid. If retaining an execution
plan is too coupled to the current global generation, first accept verified
capacity reuse and report exactly what remains rebuilt.

Required counters: surviving/new/deleted/migrated blocks; reused/rebuilt storage;
reused/rebuilt stencils and routes; moved/copied/refreshed bytes; allocations;
setup CPU/wall; live and peak bytes. Avoid a new persistent duplicate store that
trades lower allocation counts for increased memory pressure.

## C. Restore mass-tendency exchange overlap

Split native mass communication into begin/finish operations. After finishing
the fine-level divergence, begin its AT_NODE exchange; perform independent
velocity-source work; finish and install received values before the first
coarse restriction that consumes that tendency. Retain synchronous exchanges
where no independent work exists.

The previous legacy schedule already overlaps this portion of the work. The
current native routine posts and waits within a single call. Changing that
schedule is more specific than assuming all MPI time is network overhead.

Deferred operations require persistent request and send/receive-buffer lifetimes,
explicit pending state, and guarded buffer reuse. The current routine's local
request array must not escape its scope. Preserve local-copy and remote-install
ordering, message matching and the numerical phase of each dependency. Never
reuse a buffer for the next layer/stage or invalidate its plan until completion.

Measure pack/post, independent work, exposed wait and install separately. Record
messages and bytes to verify that the first experiment changes scheduling rather
than payload semantics. Forced rank delays in small MPI tests should exercise
both early and late arrival. Require no pending requests at RK-stage completion,
adaptation, checkpoint or restart boundaries.

## Direct comparison and ablations

Retain three headline numerical implementations throughout development:

1. Frozen Stage 190, without any Stage 191 patch.
2. Stage 191 candidate.
3. Frozen legacy Domain with the same existing checkpoint I/O hardening.

Build all with the same compiler and external physics objects. Optimized builds
use O2 and FP_CONTRACT=off; native CPU targeting is compiled on the node where
they run. Use a separate checked candidate and a reconstruction-validation run.
Record complete source, patch, executable and physics-object identities.

For attribution, run A-only and C-only candidates against Stage 190
and legacy, and the A+C combination. B needs its own later experiment. Each component needs its own on/off
control with default behavior recorded in the manifest. Prefer explicit build
variants for code changes; any runtime selector must be logged and rank-consistent.
All comparisons use the same checkpoint input and physics state. Do not compare
new candidate times against historical Stage 190 times from a different job.

Add a separately instrumented legacy executable using comparable phase timers:
tendency/shared producers, boundary handling, scalar restriction, inverse,
adaptation, physics, RK and output. Verify its checkpoint against the authoritative
legacy executable. Keep profiles outside performance medians. The existing
legacy instrumentation helper assumes an unpatched archived build; adapt its
source checks for the common I/O delta rather than bypassing them.

Report baseline/candidate/legacy ordinary-step and elapsed medians; candidate
reduction versus Stage 190; both block/legacy ratios; memory; construction counts;
and per-rank exposed waits. A reduction in one phase is insufficient if its
replacement setup, copying or memory cost makes the complete step slower.

## Local acceptance first

Use the existing evolved J4/J6, four-rank fixture and FP_CONTRACT=off. Complete
the optimized/checked/oracle checkpoint comparisons and the existing fixed
legacy numerical screen. Exercise refinement, coarsening, ownership changes,
checkpoint/reload, post-restart remap and appropriate RK3 smoke coverage.

Component-specific validation must cover consumed/moved storage, retained versus
invalidated dependency plans, new/deleted subtrees, changed neighbor ownership,
and outstanding MPI operations. Include poisoned stale state so an incorrectly
reused field or halo cannot silently pass. Keep the previous independent oracle.

Local testing establishes correctness and allocation/work-count changes. Do not
use whole-solver speed ratios when local compression or swapping is active.

## Cluster acceptance on nk

Keep the accepted 83-task J5/J7 RK4 configuration on bb13, one exclusive node,
one thread per rank, identical core placement, and native compilation on bb13.
Keep time_end=0.315 and dt_write=0.01 for the short checkpoint gate. Preserve the
seed checkpoint 3, one newly written/reloaded checkpoint 4, and post-restart remap.

After warmups, use at least five measured runs of each headline implementation
in a balanced rotating order within the same allocation. Retain per-rank RSS,
CPU affinity and node VM counters. Profile Stage 190, Stage 191 and legacy
separately. Reject incomplete runs, missing metrics, different printed
step/time/dt/Jmax/dof trajectories, or numerical failures before timing analysis.

For each accepted structural component, add a longer adaptive window with one
or two checkpoint/restarts. Choose its duration/save interval from an untimed
pilot, freeze them for all variants, and verify observed checkpoint coverage.
This avoids repeating the earlier hard-coded expected-checkpoint mismatch and
checks that reuse does not accumulate stale state or memory growth.

Block candidate versus Stage 190 and checked/optimized comparisons remain exact
for the existing semantic fields/topology/thresholds. Legacy retains the current
fixed numerical tolerances. No relaxed tolerance, missing restart test or oracle
omission is a performance optimization.

## Work after Stage 191

The larger architectural target is moving a complete shared producer/consumer
chain onto final block owners and removing its Domain round trip. That should
follow the paired legacy profile and the above lifetime work. Select the chain
by measured transfer/packing cost and dependency completeness; avoid a partial
port that keeps both full producers or adds another permanent state copy.

Likewise, replace leaf-count-only block weighting with measured costs only after
separating fixed Domain work from movable block work. Blindly rebalancing block
weights cannot move the retained Domain-side producers and may increase migration.
