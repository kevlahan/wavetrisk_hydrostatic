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

## Current Stage241 candidate (after committed Stage240)

Stage240 is committed as `f1bd9886`; all ten cases and source-bound review passed.
Stage241 is now implemented but unqualified. `PHYSICS_NATIVE_INPUT` sends native
column/edge keys to execution owners, fetches canonical atmospheric/mean/soil
values directly from final-owner Block_Data, and uses guarded saved final-inverse
atmospheric aliases. Production queue packing skips the legacy atmospheric pack.
Zero-key means and soil still use an explicit Domain alias adapter. Output remains
the Stage240 write-only Domain publication. This is still not complete migration.

Build is running as session 73427, driver `/tmp/build241-driver.log`, frozen tree
`.stage241-validation/candidate`. Oracle binary has compiled. No solver is active.
New helpers are under `test/parallel_block_profile/physics_241`; qualifier reuses
Stage240 injection helpers and Stage239 package-soil injection. Reviewer and audit
are drafted, and will need source/evidence checks after qualification. Next run:
`python3 test/parallel_block_profile/physics_241/qualify.py` after all builds finish.
No live source changes while the frozen candidate is qualifying.

Temporary draft generators `/tmp/edit241.py` and `/tmp/prepare241_storage.py` have
already been applied. Do not rerun blindly. Storage/exchange/alias drafts are in
`/tmp/native241_{storage,exchange,alias,fill}.txt`; live source includes subsequent
fixes. Always treat current source as authoritative.

Stage241 update: all four builds passed. Qualification runs sequentially in session
94318. The oracle-4 case passed exactly (35,812 column checks); production-4 is
running. Four additional cases are queued after inherited tests: native-route,
native-epoch, native-value, entry-poison. New injection helper and include parsed
successfully. Entry poison deliberately backs up old Domain values, poisons all
atmospheric solution plus canonical mean/soil values during native gathering, and
restores them after execute_physics_batch, before normal publication/adaptation.
This proves entry-read independence only; retain the separate Stage240 publication
poison with no backup/restoration. Reviewer includes the four extra cases and
helper/include hashes. Final review has not run.

Further lifecycle inspection: main.f90 Domain_Consumer_Audit checks transaction
phase and writeback counts only, not numerical payload checksums. Its comments
currently forbid compact transaction changes. Native physical writes need an
explicit named lifecycle contract preserving atmospheric wavelets, not silent
assumptions about numeric audits. Physical output installation must carry the
is_pole flag because pole columns do not write velocity; current 26-field owner
header does not include that flag.

For removal of the remaining physics alias mean/soil reads, capture those bounded
alias fields at the initial inverse-consumer epoch (before RK; they remain unchanged
within RK), and retain them until final physics consumption. This moves the explicit
upstream alias adapter out of physics; it does not itself finish the adaptation /
initial-input lifecycle migration. The saved atmospheric consumer alias remains
its separate final-inverse version. Do not conflate these epochs.

Stage241 full 4/41 trajectories and control have now passed exactly. Each full
run has 3,122,100 canonical native requests and 102,156 explicit alias requests.
Guard/poison tests are still running in session94318; no final review yet.

Next-change drafts (not applied): `/tmp/parallel_block_inverse_242.f90` captures
bounded alias mean/soil references at initial RK entry and offers them as optional
physics_consumer_alias results. `/tmp/native242_state.txt` writes completed
physical temperature/velocity/soil into compiled final-owner routes while leaving
wavelets unchanged. These need integration and tests, not blind copying. Native
output installation must carry is_pole; consider a 27th header field and validate
primary native column key against the execution owner. Preserve keys through result
installation instead of freeing entry_keys immediately after input fetch.

Stage241 is now fully qualified: all fourteen cases passed, final reviewer passed,
and source/main/prior evidence hashes matched. The entry-poison test exercised
20 primary alias columns (ranks0..3:6,4,6,4), so the alias reference adapter has
actual coverage. Qualification session94318 and review11999 both completed0.
No solver jobs remain active. Commit this checkpoint, then continue Stage242.
