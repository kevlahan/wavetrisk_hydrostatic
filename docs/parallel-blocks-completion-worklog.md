# Remaining migration — active implementation

User authorization: complete the remaining implementation before performance work
on bbserv. Do not stop at another partial handoff and do not launch bbserv jobs.
Use the authoritative worktree `/Users/kevlahan/.codex/worktrees/7173/wavetrisk_hydrostatic`.
The local machine has four usable cores: one solver job at a time, library threads
one; 41 ranks is an oversubscribed correctness check only.

## Accepted checkpoint

Stage239 owner package execution and its nine qualification cases are committed
as `58416b11`. The earlier retained RK migration is complete. Main is unchanged.

## Active change

Stage240 adds `PHYSICS_OWNER_BATCH`: import column records to their final owners,
recompute pressure/velocity inputs there, execute a velocity-dependency DAG, send
only required updated predecessor edges, and publish completed state in bulk.
The original serial path is the reference. The checker replays source-column
order independently and compares complete packets and double-precision state.
The production publication callback reads no numerical Domain operands.

Build currently targets `.stage240-validation/candidate`, with helper
`test/parallel_block_profile/physics_240/build_physics.py`. All ten qualification cases and final source-bound review now pass, including
exact full 4/41 trajectories and destructive poison. This checkpoint is ready to
commit before beginning native input gathering. Preserve all older frozen
evidence and do not silently refresh a frozen tree if source changes.

## Still required before claiming completion

- Qualify the owner-resident physics DAG and bulk publication, including exact
  4/41 trajectories, old-source poison, independent replay and invalid dependencies.
- Replace remaining physics entry export from Domain numerical fields with native
  owner storage and explicitly versioned alias adapters; retain updated owner
  atmospheric and soil fields across the lifecycle transition.
- Migrate the soil wavelet transform and adaptation/reconstruction numerical
  consumers. Preserve the separate pre-physics atmospheric wavelet decision image
  and post-physics physical scaling values. The inverse reconstructs fine values
  from coarse physical values plus wavelets; do not assume physics atmospheric
  updates are simply discarded.
- Replace persistent restart/output numerical mirrors with bounded owner import
  and serialization adapters. Geometry-only Domain metadata may remain.
- Remove obsolete production numerical storage and compatibility round trips only
  after all producers/consumers have moved; preserve independent oracles/off paths.
- Run complete local source-bound qualification of the final integrated route.
  Prepare a precise bbserv performance candidate after correctness is complete;
  no performance claim follows from these local oversubscribed runs.

This file is an implementation continuity record, not an acceptance claim.

## Next owner-input implementation details

Use the existing inverse operand keys for the seven column/velocity storage
identities. Positive catalogue keys must fetch physical scalar/vector, means and
soil directly from final-owner Block_Data; key-zero entry aliases must use the
saved final-inverse consumer epoch, not the tendency epoch. Preserve each alias
storage identity and expected writer ordinal. The six positive-version velocity
inputs still advance only through the DAG edges. Geometry and fixed reference
inputs may remain source metadata. Remove production pack_physics_vars reads
once native entry gathering is wired; retain the old packing in independent replay.

A native alias accessor needs explicit epoch validation and a node-slot to alias
index map; inverse_scalar_operand and inverse_vector_operand are prohibited here
because their contracts select initial/provisional tendency epochs. Current
capture_native_inverse_consumer_aliases stores the required distinct version.
Do not read current value(:,:,slot,SOL) as a substitute for that saved epoch.

After native input qualification, install physical/soil updates into owner storage
while preserving the retained atmospheric wavelet image. Amend the pending
consumer lifecycle contract explicitly; do not bypass its mutation audit. Soil
wavelets and post-adaptation inverse still need vertical-range generalization.
