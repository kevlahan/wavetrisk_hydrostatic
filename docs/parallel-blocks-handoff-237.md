# Stage237 handoff — accepted Simple-physics scratch migration

The validated owner tendency and retained RK migration is checkpointed in commit
`6e528e92` on `codex/performance-investigation`. The subsequent lifecycle audit and
first bounded migration are complete and validated, with Stage237 changes left
uncommitted for review. The authoritative worktree is
`/Users/kevlahan/.codex/worktrees/7173/wavetrisk_hydrostatic`.

## Change and scope

The source-anchored audit in `parallel-blocks-stage237-dependencies.json` covers
nine lifecycle families: RK publication, Simple physics, package initialization,
kinetic-energy ordering, other physics, soil wavelets, adaptation/rebuild, restart,
and output/checkpoint.

`PHYSICS_COLUMN_LOCAL=on` computes Simple-physics hydrostatic pressure and
geopotential in the active column, preserving the original arithmetic and traversal
order. It removes this consumer's dependence on Domain hydrostatic scratch and
creates no Domain-shaped numerical snapshot. Integer visitation metadata guards
against duplicate consumers. `check` additionally evaluates independent legacy
hydrostatic values. The flag defaults to `off`.

Physics execution, prognostic and soil state, geometry, and kinetic-energy stencils
remain Domain-owned. Persistent scratch allocations remain needed by other
consumers and validation paths. This is not complete physics owner migration.

## Acceptance evidence

`parallel-blocks-stage237-local-results.json` was produced successfully by
`test/parallel_block_profile/physics_237/review_physics.py` after all cases finished.
The reviewer binds evidence to live/frozen source, binaries, inputs and helpers.

- Full production runs at 4 and 41 ranks match the Stage236 same-rank reference
  exactly: 12 steps, two remaps, restart/checkpoint and 62 VTK members.
- The four-rank checked oracle passed 1,074,360 layer comparisons, comparing five
  double-precision values per layer, plus surface pressure and RK alias checks.
- Destructively poisoned legacy physics scratch stayed poisoned without a backup
  or restoration, while the completed step matched the reference exactly.
- Missing epoch, duplicate consumer and corrupted pressure injections failed at
  their intended guards; the feature-off control passed and matched the reference.
- Rank 40 in the 41-rank case has no Domains while participating in retained block
  RK work. Physics still follows Domain ownership.

The frozen Stage236 sources and all four binaries remain unchanged. Main remains
`6a245db7f1632b0205c2bebe95e5f0cd535c2582`. No remote jobs or pushes were performed.

## Four-core local execution constraint

This machine has four usable cores. Early four-rank checks overlapped; subsequent
cases were serialized after the user's correction. Numerical-library threads were
limited to one. All cases have completed and no processes remain paused. Future
local solver jobs must run one at a time. The oversubscribed 41-rank test exists
only to exercise zero-Domain correctness; no timing or scaling claim is made.

## Next migration boundary

Before moving Simple physics to final block owners, explicitly represent the
versions and order of neighboring velocity dependencies: kinetic energy reads
velocities that earlier columns may already have updated. Preserve the pending
atmospheric wavelet image and soil epochs through that transition. Adaptation,
restart and output retain separate lifetimes to migrate afterward. Keep new flags
opt-in until those boundaries receive independent numerical qualification.
