# Remaining migration — current implementation record

User authorization: **complete the remaining implementation before performance
work on bbserv**. Do not stop at another partial handoff or launch bbserv jobs.
The authoritative worktree is `/Users/kevlahan/.codex/worktrees/7173/wavetrisk_hydrostatic`,
branch `codex/performance-investigation`. The tool default cwd26d8 is not this tree.
The local machine has four usable cores: one solver job at a time, library threads
one, builds `-j1`. Forty-one MPI ranks is oversubscribed correctness coverage only.
No subagents. Main must remain `6a245db7f1632b0205c2bebe95e5f0cd535c2582`.
Mutations/tests in7173 require the approved escalation mode. Never stage all files;
many unrelated historical artifacts and generated package objects are untracked.

## Accepted checkpoints

- Stage236 retained RK/tendencies: `6e528e92`.
- Stage237 local physics pressure scratch: `efdd9b9b`.
- Stage238 explicit velocity versions: `1a63854f`.
- Stage239 ordered owner package calls: `58416b11` (nine accepted cases).
- Stage240 owner physics DAG/bulk publication: `f1bd9886` (ten accepted cases).
- Stage241 native canonical input gathering/saved atmospheric aliases: `0a69c3ba`
  (fourteen accepted cases, including destructive entry poison).

Frozen candidate trees, binaries and trajectories are under `.stageNNN-validation`.
Never silently refresh a frozen tree after source edits. Preserve failed attempts.
Stage240 first failed replay consumed pressure scratch twice; resetting surface
pressure/geopotential before independent replay fixed it. Failed attempts are kept.
All accepted full runs match12 step states, two remaps and62 VTK members at4/41 ranks.

## Stage242 — accepted; ready to commit

Five solver files changed: Makefile, src/main.f90 comments, src/parallel_block.f90,
src/parallel_block_inverse.f90, src/physics/physics_simple.f90.
`PHYSICS_NATIVE_STATE=on/check` requires native input and retained RK.

- Alias mean/soil references are captured at initial RK entry; those fields remain
  unchanged within RK. Atmospheric aliases use their distinct saved final-inverse
  consumer epoch. Production physics now reads no Domain prognostic/mean/soil
  input; geometry and fixed reference inputs remain geometry metadata.
- After complete DAG execution, canonical temperature/velocity/soil results are
  retained on final owners. Primary alias results are retained in native alias
  storage. Header27 carries is_pole, preserving non-writing pole velocity.
- Mass, means and atmospheric decision wavelets remain unchanged. Main's consumer
  audit preserves pending phase and writeback counts; its comment now explains
  native physical updates. The audit does not hash numerical payloads.
- The check build reads stored values back, checking exact scalar/vector/soil,
  mass/mean/pole invariants; independent serial replay still checks scientific
  equivalence of complete packets and results.
- Final Domain publication still serves legacy soil/adaptation. Native state is
  rebuilt from Domain after adaptation; the complete lifecycle is NOT migrated.

All four builds passed. Frozen tree `.stage242-validation/candidate`, build driver
`/tmp/build242-driver.log`, completed build session92112. Qualification session85276 completed0. All sixteen cases and final reviewer
(session5439) passed. No solver/build jobs remain active. Commit this checkpoint,
then continue the soil/range migration.

Passed: oracle-4 (35,812 native storage readbacks), full4/41, control, route,
dependency, version, output, publication poison, package-soil poison, native-route,
native-epoch, native-value, entry-poison, retained-owner and retained-alias. Each full run retains
460,368 canonical and240 alias columns. The alias test now logs `Inject retained
alias corruption` before mutation; helper and reviewer require the marker so its
shared mismatch diagnostic cannot yield a false positive.

Helpers: `test/parallel_block_profile/physics_242`; audit, handoff and final
source-bound results JSON are complete.
The entry poison now poisons ALL old atmospheric/mean/soil inputs while gathering,
then restores test buffers before legacy consumers. It exercises20 primary alias
columns. Distinguish this from publication poison, which has NO backup/restoration.
Helpers/include parse successfully. Keep source frozen until final review. Update
handoff/worklog and commit focused paths after acceptance, then continue remaining
work; do not end the user's task after this checkpoint.

## Remaining implementation

1. Native soil wavelet transform/inverse, preserving atmospheric decision wavelets.
2. Native numerical adaptation/reconstruction and topology migration. Geometry-only
   Domain metadata may remain, but numerical blocks must not rebuild from mirrors.
3. Replace the initial alias import and restart/output mirrors with explicitly
   versioned bounded adapters. Output may mutate aliases/diagnostic scratch: audit
   actual behavior rather than assuming serialization is side-effect free.
4. Remove obsolete persistent production numerical storage only after all producers
   and consumers have moved. Preserve off paths and independent oracles.
5. Complete source-bound integrated qualification, then prepare the bbserv candidate.
   No remote jobs or local performance claims yet.

## Scientific/lifecycle constraints

Atmospheric wav_coeff is the PRE-physics decision image, while physical scaling
is POST-physics. The legacy inverse reconstructs from coarse physical values plus
wavelets; do not assume physics updates are discarded. Soil WT_after_step computes
coefficients without scalar restriction, compresses, then inverses. Keep that order.
Native inverse consumer aliases are deliberately older than fresh kernel aliases:
post-prolongation scalar aliases, and pre-finest-inner vector aliases. Preserve
those versions for downstream consumers, including a soil inverse.

The initial/provisional accessors inverse_scalar_operand/inverse_vector_operand
must never substitute for the final consumer alias. Stage241/242 physics uses
physics_consumer_alias with readiness and registration guards. The new result
setter updates only physical alias fields; finish closes entry access afterward.

## Next soil/range work — findings and drafts

`/tmp/native243_transfer.txt` is an UNAPPLIED draft adding optional field_first /
field_count to transfer_local_block_inverse_routes, keeping atmospheric defaults.
Its offsets/counts need validation; do not treat the draft as tested code.
Old generators `/tmp/edit241.py`, `/tmp/edit242.py`, `/tmp/prepare241_storage.py`
were already applied. DO NOT rerun them. Live source is authoritative.

Native inverse has68 zlevels references. Stage203 scalar/inner kernels use node
slots and runtime layer loops, so topology is mostly vertical-independent. Add an
explicit selected range, updating transfer counts, address(5), buffers and epoch
contracts together. Do not mutate the global zlevels parameter. Atmospheric
physics/tendency accessors must reject soil-only epochs; soil publication must
not mark the atmospheric physics epoch ready.

prepare_native_inverse currently requires rk_boundary_ready and reads wavelet
aliases from its Domain wavelet argument. A soil entry mode needs native physical
and wavelet seeding, not a fake RK stage. Preserve/export post-physics soil aliases
before a range-plan rebuild reallocates consumer arrays. For generic soil support,
retain all soil scalar/vector aliases, not only S_TEMP. Do not silently assume all
other soil fields are zero (including signed-zero behavior).

produce_block_scalar_wavelets (parallel_block_mpi.f90 around16353) and
compute_block_scalar_wavelet (around16590) already address arbitrary block field
levels, but producer loops skip outside1:zlevels. Inactive soil geometry records
can be zero: explicitly select the physical geometry class for geometry, while
reading numerical values at their actual soil level. The same issue applies to
vector operators. zero_resident_wavelets currently zeros entire arrays: a soil
path must zero ONLY soil coefficients to protect atmospheric decision wavelets.

Stage242 reference arrays currently allocate/capture at initial native RK entry
even when PHYSICS_NATIVE_STATE is off (numerically inert). Guard that cost when
preparing final performance baselines; do not change frozen242 source mid-review.
