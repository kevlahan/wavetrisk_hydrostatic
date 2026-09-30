# Stage238 — accepted explicit Simple-physics velocity operands and versions

All qualification cases and the source-bound `physics_238/review_physics.py`
passed. The acceptance report is `parallel-blocks-stage238-local-results.json`.
The prior scratch migration is committed as `efdd9b9b` on
`codex/performance-investigation`; this new work remains uncommitted for review.
The authoritative worktree is
`/Users/kevlahan/.codex/worktrees/7173/wavetrisk_hydrostatic`.

## Bounded change

`PHYSICS_VELOCITY_PLAN=on/check` requires `PHYSICS_COLUMN_LOCAL=on/check` and
defaults to off. It records the original column traversal and six velocity reads
per column as integer metadata. Each read has a local source identity and expected
writer ordinal: zero denotes the value at entry to the physics step; a positive
ordinal denotes a preceding non-pole column's write. The plan is rebuilt each
physics step and released before soil transforms and grid adaptation.

The planned edge indices supply a six-value local operand array for the pure
`physics_kinetic_energy` kernel. Its arithmetic order matches the original
`diagnostics_mod::kinetic_energy`. The three current-column values also supply
U, V and W to the unchanged physics package. The `check` variant compares every
layer's double-precision kinetic energy bit-for-bit with the independent original.
Guards reject an absent epoch, changed consumer order, mismatched source versions,
duplicate planned columns, and incomplete consumption.

| Slot | Local source | Component |
| --- | --- | --- |
| 1 | current column | UP |
| 2 | current column | DG |
| 3 | current column | RT |
| 4 | south | UP |
| 5 | southwest | DG |
| 6 | west | RT |

Non-pole columns overwrite all three current edges after the package call, at all
atmospheric levels. Poles read the same stencil but do not overwrite velocity.
Temperature and soil updates retain their existing locations and lifetimes.
The source-anchored audit is `parallel-blocks-stage238-dependencies.json`.

## Ownership and ordering limits

This change extracts explicit operands and makes their version contract executable.
Physics execution and numerical source fields remain on Domain owners. Metadata
uses local storage identities, not canonical block-owner identities. An alias
without a preceding local writer must retain its entry value even when another
copy has been updated elsewhere. A future exchange must request the correct
version rather than automatically choosing the newest owner value.

The full original order is retained, preserving both read-after-write and
read-before-write dependencies. This is not yet a parallel dependency scheduler.
The pending atmospheric wavelet transaction and the soil transform are unchanged;
no global Domain storage is removed. No performance improvement is claimed.

## Local qualification protocol

The four-core machine runs one solver job at a time, with numerical-library
threads set to one. `physics_238/qualify.py` enforces sequential execution of an
oracle case, full production cases at 4 and 41 ranks, a feature-off control, and
four isolated fault injections. The 41-rank case is deliberately oversubscribed
to exercise zero-Domain correctness; it is not a scaling measurement.

Each full case must match Stage237 at the same rank count for 12 steps, two remaps,
restart/checkpoint and 62 VTK members. The short oracle and control must match the
first reference state. Missing epoch, wrong column order, stale source version,
and corrupted velocity operands must fail with their specific diagnostics.
The reviewer checks source and binary identities and preserves the frozen
Stage237 evidence. No remote job or push is part of this work.

## Accepted evidence

- The four-rank oracle passed 1,074,360 exact kinetic-energy comparisons alongside
  the existing hydrostatic and RK alias oracles. Its 35,812 columns used 134,781
  entry-version edge reads and 80,091 reads from preceding updates (counts are
  per column, before the atmospheric-layer multiplier).
- Both full production cases passed all twelve exact step comparisons, two
  remaps, checkpoint reloads 2 and 3, and 62 exact VTK members each. The four-rank
  trajectory consumed 460,608 physics columns and 1,033,282 updated edge reads.
- The 41-rank case recorded zero physics columns on rank 40 at every step while
  that rank owned block RK work.
- The feature-off control matched the first reference state exactly.
- Missing epoch, wrong consumer order, stale source version, and corrupted
  velocity operands each failed at their intended diagnostic.
- The final reviewer verified live/frozen sources, binaries, input/log hashes,
  helper identities and the source-anchored audit. Frozen Stage237 source and
  binaries remain unchanged. Main remains
  `6a245db7f1632b0205c2bebe95e5f0cd535c2582`.

All solver cases ran sequentially and have completed. No paused cases remain.
The new flag defaults off; this is correctness qualification, with no performance
or scaling claim. `git diff --check` passes.

## Next boundary

Map these local source/version pairs to canonical block-owner operands while
preserving stale entry aliases and required writer ordering. Audit the external
physics package's mutable state and first-call behavior before changing execution
ownership. Keep the atmospheric pending decision image and soil epochs separate.
