# Stage239 — accepted ordered package execution on block owners

The validated velocity-plan migration is committed as `1a63854f` on
`codex/performance-investigation`. This owner-execution change is uncommitted.
The authoritative worktree is
`/Users/kevlahan/.codex/worktrees/7173/wavetrisk_hydrostatic`.
All nine qualification cases and `physics_239/review_physics.py` passed.
The source-bound acceptance report is `parallel-blocks-stage239-local-results.json`.

## Bounded ownership change

`PHYSICS_OWNER=on/check` requires `PHYSICS_VELOCITY_PLAN=on/check` and defaults off.
Catalogue-covered physics columns execute the complete external single-column
package on their final block owners. Integer subtree maps associate local patches
and nodes with catalogue IDs. Fixed scaffold columns remain local.

The geometry source gathers each column in the original order using the validated
velocity-version plan, supplies a bounded single-column packet, and waits for the
owner's result before publishing velocity, temperature and soil updates. Sources
are serviced serially; all other ranks service requests, including ranks with no
Domains. The protocol validates destination ownership, source rank, writer order,
source versions, message extents and reply identity.

Each velocity operand carries its local source-node identity, expected preceding
writer ordinal, and catalogue ID when covered. Index zero explicitly retains
compatibility aliases/scaffold. This maps covered sources to owners without
silently replacing stale entry aliases with newer canonical values. It does not
yet move the version authority or numerical gather to owner storage.

This is a conservative execution migration. The numerical inputs and returned
prognostic/soil fields remain Domain-owned. Packets retain the source column
orientation; they are not direct writes into canonical owner field storage. No persistent Domain arrays are
removed, no parallel dependency scheduler is introduced, and no performance
improvement is claimed. The pending atmospheric wavelet transaction and soil
transform stay in their existing lifetimes.

## Package-state audit

`parallel-blocks-stage239-dependencies.json` binds the audit to source files.
Every packet supplies latitude/longitude, time, timestep, mask, first-call flag,
atmospheric inputs and the complete surface/soil column. Location scratch is
reloaded before execution; the package wrapper reloads soil for non-first calls.
The original per-column first-call semantics are preserved.

Package allocation is private and rank-local. A previously idle owner receiving
its first non-first-call packet initializes through a disposable first-call copy,
then reloads the requested inputs for the real call. Orbit/surface parameters are
initialized as in the original path. `icount` and `zDay_last` do not feed numerical
outputs in this package. The empty `physics_checkpoint_restart` routine is not
mistaken for an allocator. No third-party package source is modified.

The check variant repeats the same column package call locally and compares the
entire returned single-precision packet bit-for-bit. Full trajectory checks remain
independent of that transport oracle and compare against Stage238.

## Qualification protocol

One solver job runs at a time on four usable cores; numerical-library threads are
one. The 41-rank case is oversubscribed solely for zero-Domain correctness.
The sequential suite covers a four-rank oracle, full 4/41-rank trajectories,
feature-off control, wrong destination, reordered request, invalid source version,
corrupted output, and a two-step destructive package-soil poison test. The latter
also forces the late first-use initialization path on non-first calls, verifying
that its disposable initialization cannot replace the requested soil state.

Full cases must match twelve steps, two remaps, checkpoint reloads 2 and 3, and
62 VTK members at each rank count. Rank 40 must execute remote physics while
retaining zero Domain columns. Owner and scaffold execution totals must equal the
original source column count. All evidence is bound to frozen source, package
source, binaries, inputs and logs by the final reviewer.

## Accepted evidence

- The checked four-rank run matched all 35,812 returned column packets exactly,
  including 962 remote executions, while retaining the hydrostatic, velocity and
  RK oracles. Its completed state matched Stage238.
- Both full trajectories matched all twelve steps, two remaps, checkpoint reloads
  2 and 3, and 62 VTK members exactly. Each executed 460,608 owner columns. Remote
  execution counts were 14,796 at four ranks and 210,104 at 41 ranks.
- Rank 40 retained zero Domain columns while executing 5,824 remote physics
  columns across the 41-rank trajectory; it executed owner work on every step.
- The feature-off control matched its reference state exactly.
- Wrong destination, reordered request, invalid source version and corrupted
  output each failed at the intended guard/oracle.
- The two-step soil poison and forced late-initialization case completed and
  matched both reference states exactly, without retaining a Domain-sized backup.
- The final reviewer verified live/frozen sources, external package sources,
  binaries, inputs, logs, fault injections and helper identities. Frozen Stage238
  source and binaries remain unchanged. Main remains
  `6a245db7f1632b0205c2bebe95e5f0cd535c2582`.

All solver tests ran sequentially and have completed. No processes remain paused.
No remote job or push was performed. `git diff --check` passes. The qualified
fixture has no scaffold physics columns, so the explicit scaffold fallback is not
separately exercised by these cases. No performance or scaling claim is made.

## Next boundary

Move versioned gathering and returned atmospheric/soil state onto owners, retaining
explicit compatibility adapters for entry aliases and fixed scaffold. Then replace
the serial source service with a dependency-aware schedule; retain the existing
column ordering until independently qualified. Do not infer a speedup from this
correctness-first transport path.
