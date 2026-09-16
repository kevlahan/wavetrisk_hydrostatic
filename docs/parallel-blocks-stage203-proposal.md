# Stage203 proposal: eliminate repeated representation passes

Status: original design rationale. The implemented scope and retained paths are documented in
[Stage203 implementation](parallel-blocks-stage203.md); cluster job 54695 now shows a real gain below the acceptance threshold. See
[the review](parallel-blocks-cluster-203-54695.md).
Baseline: commit `7c406f4d`, built with `STAGE202=on`, the same earlier runtime
selectors, and `FP_CONTRACT=off`. The default `STAGE202=off` build is NOT the
validated fast baseline. Prerequisite diagnostics are in `bd120da0`.

## Evidence and size of the opportunity

Stage202 job 54693 reduced ordinary-step time from 1.055 to 0.925 s, versus
0.5415 s for legacy. The remaining gap is 0.3835 s/step. Both reversed pairs
recovered more than 0.12 s; all 17 runs and 37 correctness/measurement gates pass.
See [the review](parallel-blocks-cluster-202-54693.md).

The observed savings accompany fewer instructions: wavelet instructions fell
14.1–14.9%, and other-dynamics instructions fell 17.5–18.0%. Wavelet CPI was
nearly unchanged; other-dynamics CPI increased slightly. This supports removing
more representation work. It does not establish how much of each surviving
operation is removable; user-space MPI progress is included in these counters.

| Phase | Stage202 rank-mean time | Earlier legacy rank-mean time | Indicative excess |
|---|---:|---:|---:|
| Dynamics wavelets | 0.2445 s | 0.0691 s | 0.1754 s |
| Other dynamics | 0.4853 s | 0.3750 s | 0.1102 s |
| Adaptation/reconstruction | 0.1497 s | 0.0600 s | 0.0897 s |

These legacy phase controls come from Stage201. They rank opportunities; they
are neither a same-job residual attribution nor guaranteed removable seconds.
Their rank-mean statistic is distinct from the 0.3835 s printed-step median gap.
Stage203 should include one fresh legacy coarse control to resolve that limitation.

## 1. One execution representation for the complete inverse transform

### Current, verified source path

`activate_block_native_inverse_transform` in `src/parallel_block_mpi.f90`
executes scalar prolongation/reconstruction and inner-vector reconstruction on
block fields. Outer-vector reconstruction executes in `parallel_block_inverse_mod`
on its Domain-indexed `value` workspace. Between kernels, the synchronization
routines gather block fields into that workspace, apply signed aliases, scatter
back to block boundaries/interiors, and refresh block ghosts.

Stage202 removed the self-route wire-buffer round trip, but still copies between
block storage and the workspace. At each inverse level the code performs two
scalar synchronizations, one outer-vector synchronization and one inner-vector
synchronization. There is also an initial scalar synchronization. Each invokes a
whole scalar/vector ghost-payload refresh; the current API has no changed-grid-level
argument. For the J5/J7 interval this is nine such calls per inverse transaction,
excluding preparation and final compatibility work. This is a source call count,
not nine measured identical costs or a count per whole time step.

Preparation additionally gathers both field families and both components. Final
publication exposes compatibility state to Domain consumers. Production does
not call a second complete legacy inverse calculation; the phase oracle does
independent checking only when enabled.

### Proposed modification

Keep the existing inverse workspace as the execution representation for the whole
transaction, and compile scalar and inner-vector stencils into it, alongside the
existing outer-vector tape. Import authoritative block inputs at entry; install
produced block results at exit. Within the transaction, all three kernels consume
the same prepared addresses and dependency values. Remove the intermediate
block-to-workspace gathers, workspace-to-block scatters and block-ghost refreshes
whose only consumers were those inverse kernels.

This is a bounded prototype using an existing workspace, not a proposal to allocate
another complete solution. Retain compact block storage as the authoritative
representation outside the transaction. Preserve the current arithmetic expression
and per-output accumulation order, including contraction off. Compile topology
addresses once per generation; read changing masks and coefficients at their
proper epoch. Reuse Stage193 scalar-column and Stage194 inner-geometry arithmetic
rather than substituting a different interpolation algorithm.

This relocates scalar/inner-vector execution from block owners to the owners of
the existing inverse workspace. That is a real ownership/load-distribution change,
not a free optimization. Record per-rank stencil work and coarse phase maxima;
reject the path if a worse distribution cancels the eliminated passes. A later
block-owner implementation would be a separate design, not an invisible change
to the ownership convention in this prototype.

### Required semantics before removing any pass

- Current gathers restore fixed scaffold values and alias-overwritten interiors.
  Skipping them without replacement is incorrect. Distinguish produced owner
  values from alias destinations; retain only the required restoration/snapshot
  subset, with explicit phase epochs. Prove source/destination overlap handling
  against the existing ordered alias implementation.
- Remote alias sources are snapshotted before local alias writes. Preserve signed
  destinations, local dependency chains, receive installation order and pentagon
  corrections. An alias update is not necessarily an unordered copy operation.
- Preserve `consumer_alias_value`: downstream physics deliberately sees different
  scalar/vector boundary epochs from the latest inverse-kernel dependencies.
- Keep final Domain publication at the real RK/physics/adaptation consumers.
  Do not assume every writeback or scaffold node is redundant.
- In checked mode, compare after each scalar/outer/inner phase as well as at the
  final checkpoint. Required intermediate publications remain available to the
  oracle even when the production path no longer needs them.

The intended saving is whole local pack/copy/install/traversal passes plus their
secondary waiting, not an assumed benefit proportional to fewer transmitted bytes.
The 0.1754 s phase excess is an opportunity envelope; it is not the expected saving
from this redesign alone.

## 2. Patch-wide scalar field operations across the remaining pipeline

### Current, verified source path

The scalar store already separates 17 field slots from 33 shared geometry slots.
Stage202 accelerates the virtual 50-slot interface and avoids a full record in
direct flux. Several production paths still reconstruct virtual addresses and
translate slots for each sample:

- `finalize_block_scalar_divergence_capture`: compact received physics triplets
  are installed with one range-write per sample.
- `exchange_block_scalar_restriction_ghosts`: local sources are packed into a
  scratch buffer then installed; dynamic ghost fields use sample-wise reads/writes.
- `pack_production_restriction_ghost` and `install_production_restriction_ghost`:
  indexed scalar adapters translate selected fields for every sample. Stage202's
  range batching does not eliminate these indexed adapters.

The physics producer already writes patch-wide fields through
`scalar_write_field_records`. Shared geometry and producer transport are already
deduplicated. Reimplementing those completed improvements would not address the
residual gap.

### Proposed modification

Extend the existing patch-field interface through receive installation, direct
flux, restriction, divergence and ghost handling. Compile logical field selections
to compact field-row indices once; pass patch/column views plus separate geometry
views into kernels. Bulk-install received physics fields. For same-rank ghosts,
copy the required compact field rows directly, eliminating the temporary wire
representation; preserve its snapshot if overlapping storage requires it.
Remote routes pack/unpack the same selected rows in the existing wire order first,
so this prototype can test local work removal without changing the network protocol.

Retain the expanded record API for the independent oracle and compatibility callers.
Keep physical/inactive geometry classes, signed zeros, soil/surface semantics,
mask freshness and geometry ownership. Only remove zero/poison initialization
where a write-before-read proof covers every consumer; no blanket deletion of
initialization or validity checks. Avoid noncontiguous array temporaries that
would silently reintroduce the copies being removed.

Do not independently promise the previous 0.0895 s other-dynamics improvement
again. That gain is already in the Stage202 baseline. This change targets the
surviving generic-adapter and local staging loops, whose current exact cost has
not yet been isolated.

## Scope and acceptance

Implement the inverse transaction and scalar field path as independent prototype
selectors, retaining Stage202 fallbacks. Unit/phase tests can isolate each change;
combined performance must pass the substantial-gain screen. Do not add an unrelated
adaptation redesign to Stage203: its approximately 0.090 s excess did not improve
in Stage202 and needs its own concrete consumer/lifetime analysis.

Before cluster submission, record a source-level before/after ledger of removed
passes, field values copied and adapter calls. In a diagnostic build only, count
those operations without per-node timers. No claim of removing a pass is accepted
if a new array temporary or workspace copy merely replaces it. Use local J4/J6
restart/remap and phase-oracle tests for correctness, not cluster speed predictions.

On bb13, compare freshly built Stage202-on, Stage203 and legacy with identical
physics objects, contraction off, fixed 83-rank 42/41 placement and the existing
checkpoint window. Keep the checked oracle-off/on gates and two reversed bare
pairs. Reuse the minimal Stage202 coarse/core-PMU observers, adding a legacy coarse
control rather than restarting the scaling, stack-sampling or memory-event surveys.
Report instruction counts, phase means/maxima, observer calibration, placement,
checkpoint/workload comparisons and copied-value accounting separately.

**Minimum success: at least 0.10 s/ordinary step recovered in each reversed pair,**
with all exact block/oracle/restart gates and the fixed legacy numerical screen
passing. This would move 0.925 to at most approximately 0.825 s/step, removing
26% of the present gap and bringing blocks/legacy to approximately 1.52x if legacy
stays at 0.5415 s. An ambitious design goal is 0.15 s/step (39% of the gap), not a
forecast. If the prototype cannot meet the minimum, preserve the measurements and
reconsider the representation choice rather than promoting a marginal gain.

Use the same single review tarball and scp collection process as Stage202.
