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

## Stage242 — accepted d0bf509e

Native physics output retention and alias mean/soil inputs passed sixteen cases.
Frozen evidence: `.stage242-validation/candidate`. Full4/41 trajectories match
12 states, two remaps and62 VTK members.

## Stage243 — accepted

Native soil selected-range transform/inverse passes all ten source-bound cases.
Frozen evidence: `.stage243b-validation/candidate`; results and kernel audit are
in docs/parallel-blocks-stage243-{results,kernels}.json. Full4/41 and full Domain
poison trajectory match12 states, two remaps and62 VTK members. Exact owner
readbacks pass4/41; stale epoch, invalid range and alias corruption are rejected.
The 15 soil kernels normalize exactly to legacy arithmetic. Preserve geometric
alias traversal: canonical-only inverse traversal missed alias writes.
Failed attempts remain under .stage243-validation and reject-poison. The latter
had a test backup capacity error; corrected independent allocation backups and
slice restoration pass the complete trajectory. No active solver jobs remain.

## Remaining implementation

1. Native numerical adaptation/reconstruction and topology migration, without
   rebuilding prognostic blocks from Domain numerical arrays.
2. Initial aliases and restart/output bounded adapters; audit output mutations.
3. Remove obsolete persistent production numerical mirrors; retain off/oracle.
4. Integrated qualification and bbserv candidate preparation; no remote runs yet.

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

## Stage243 — in progress, not accepted

Live source changes add selected vertical ranges to native bulk inverse storage,
explicit post-physics soil alias export/entry, native soil forward/compression and
Stage203 inverse. `SOIL_NATIVE=on/check` requires native physics state and Stage203.
The check path executes independent legacy soil transform then exact comparison.
No adaptation changes yet. Remaining lifecycle scope above still applies.

- `.stage243-validation/range-core`: pre-driver range-generalization checked build
  and `range-core-4` one-step atmospheric regression PASS (all Stage242 checks).
- `soil-attempt1`/`soil-attempt2`: preserved failed Fortran syntax builds.
- `soil-attempt3`: checked native build passed; first exact soil run
  `soil-attempt3-4` in progress (session40739, log /tmp/run243-soil3.log).
- Build helper `test/parallel_block_profile/soil_243/build_soil.py` is new/untracked.
  It preserves frozen sources; do not refresh earlier trees. Other helper stages
  remain frozen. Build sessions31472,14613,34994,13559 finished (last succeeded).
- Temporary generators /tmp/edit243_core.py, /tmp/add243_soil.py,
  /tmp/fix243_setters.py have ALREADY BEEN APPLIED. Do not rerun.
  Live source is authoritative. Core script originally had invalid Python '&'
  dead code, fixed before applying. Setters are explicit calls because single-line
  Fortran IF rejects pointer-function assignment. Original legacy arithmetic has
  been retained with native storage getters/setters.
- No new checkpoint commit yet. Do not claim the soil path verified until the
  exact oracle and full trajectory qualification actually pass.

Stage243 update: attempt3 exact oracle failed on key-zero scalar alias nodes
(level5, differences up to2.8e-11). Do NOT loosen tolerance. The canonical203 tape
omits geometric alias writes performed by the legacy soil inverse. Attempt4 uses
complete native getter/setter traversal of legacy inverse geometry, preserving
those writes. Attempt4 checked build and one-step exact soil entry/output PASS.
The current candidate adds explicit output-check markers, guards reference arrays
when native physics state is off, and keeps only temperature soil references when
soil migration is off. Frozen `.stage243-validation/candidate` builds native,
production,control sequentially (session85305; /tmp/build243-final.log). Native and
production built; control finishing at last check. Run helper soil_243/qualify.py
only after all builds finish. It runs 8 cases sequentially, full12step4/41 cases
plus exact oracle4/control4 and4 adversarial cases. One-step helper explicitly
compares step states; full helper also checks remaps/checkpoint/62VTK payloads.
Temporary /tmp/inverse243_soil.py ALREADY APPLIED, never rerun. Source frozen during
qualification. Need accepted result docs/reviewer and commit, then continue full
remaining lifecycle work. No bbserv runs and no final partial handoff.

Stage243 final qualification currently active: session86168, /tmp/qualify243.log.
Candidate build completed all3 variants. `oracle-4` PASS with exact soil entry and
output checks; full production-4 active at last inspection. Remaining qualifier
cases production41/control4/epoch/range/alias/poison follow serially. Poison test
now runs FULL12step4rank trajectory, checks2remaps/62VTK/checkpoint through existing
validate(), and retains that result in its manifest (not just one-step mass).
Reviewer soil_243/review_soil.py prepared but MUST run only when all8 cases finish.
No solver source edits until this source-bound qualification completes.

Next adaptation design notes (NOT IMPLEMENTED; verify before acting):
- `install_block_native_adaptation_mask_seed` already uses native coefficients,
  then clears all block state. Capture a versioned native numerical transition
  image BEFORE that clear; topology currently extends stable Domain node ids.
- `adapt` still calls Domain compress_wavelets and analytic update_case.
  Climate update_case initializes reference means/topography on new patches and
  refreshes their boundaries; it does not evolve prognostic sol. Treat this as
  a bounded reference initialization adapter, not a prognostic source.
- New source blocks should be populated from the native transition image (old
  physical values and PRE-physics atmospheric wavelets, post-soil wavelets),
  zeroing genuinely new nodes like extend_Domain does, then migrated before
  native compression/reconstruction on the new topology. Do not silently use
  Domain sol/wav fields as fallback. Existing `build_source_blocks` has many
  scalar/vector/ghost/boundary reads and comparisons; use explicit native field
  provider callbacks in all of them. Preserve Domain off/oracle path.
- Before soil plan rebuild, archive atmospheric consumer aliases and coefficient
  aliases; the soil range rebuild otherwise destroys those saved arrays.
- Full native inverse traversal used by soil now exactly preserves geometric
  alias writes. Generalize that traversal for post-adaptation reconstruction;
  Stage203 canonical tape alone failed the exact soil alias oracle. Keep bounded
  numerical transition storage temporary; persistent canonical state remains on
  final-owner blocks. Further remove initial alias import and output/restart
  mirrors only after their consumers move.
- Makefile orders inverse before build module, so build can depend explicitly on
  a native provider if no module cycle is introduced. Main currently invokes
  Domain inverse, closes adaptation, builds/migrates blocks, then retains state;
  native ordering must build/migrate from native old state first, reconstruct,
  then close/retain. Recheck lifecycle flags around clear/migration carefully.

Source arithmetic audit soil_243/audit_kernels.py PASS: all15 copied native kernels
normalize IDENTICALLY to legacy math when getter/setter addressing is reversed.
Result docs/parallel-blocks-stage243-kernels.json. Include this in final reviewer.
Qualification oracle4 and fullproduction4 PASS. Fullproduction41 active now.
Important next-step correction: extend_Domain new atmospheric scalar defaults are
NOT always zero. If split_mean_perturbation=false, both scalar defaults are
`a_vert_mass(k)+b_vert_mass(k)*p_0/grav_accel` in compressible mode; vectors and soil
are zero. Preserve these exact defaults in any native topology extension.

Final243 first candidate qualification stopped at control4: physics_consumer_alias
unconditionally required reference readiness even when mean/soil outputs were NOT
requested. New off-storage guard exposed that latent coupling. Fixed readiness
check to apply only when optional mean or soil is present. Preserve entire old
.stage243-validation tree, including passing full4/41 and failed control. Fresh
final tree is .stage243b-validation/candidate; qualifier/reviewer ROOT updated.
All8 cases must rerun against freshly bound binaries. Re-run kernel audit too.

Fresh243b serial build active session77966 (/tmp/build243b-final.log); no solver
currently active. Previous qualifier86168 ended nonzero at old control4. All old
.stage243-validation evidence remains untouched. New qualifier now runs control4
FIRST, then exact oracle4/full4/full41/four fault cases. Fresh range-core and
soil math are unchanged except conditional optional-reference readiness guard.
/tmp/native244_provider.txt is an UNAPPLIED draft for a native source-block
serialization provider (range/scalar/subtree access), requiring transition-state
APIs. Existing builder has28 sol/wav reads:4 whole-field subtree copies,20 range
reads,4 scalar sample reads. All must switch under explicit native-transition
selection; reference means remain a separate bounded initialization boundary.
`native_velocity_plan_generation()` already exposes current block writeback
plan generation for native inverse re-entry. `retain_post_grid_change_block_reconstruction`
also unconditionally reimports Domain SOL/WAV boundary families; native adaptation
must replace that import, not merely switch initial source-block construction.

Fresh243b build completed all3 variants. Qualification active session40381,
/tmp/qualify243b.log. Corrected control4 PASS, exact soil oracle4 PASS; full4 active
at last check, full41/faults follow. Prior build77966 ended0. Do not re-poll old
sessions. No Stage243 acceptance/commit yet; reviewer waits for all8 cases.
Additional UNAPPLIED draft /tmp/native244_transition.txt provides sparse native
atmospheric/soil transition capture + field getter/defaults/epoch guards. It needs
module declarations/imports, explicit ADAPT_NATIVE guards/hooks, full-range native
inverse entry and native builder selection before use. Draft new-node defaults
must enforce compressible scope. Do not apply blindly. Source providers in
/tmp/native244_provider.txt require these APIs. No Stage244 live source changes.

Further adaptation audit (not implemented): mask_active(block_seed_installed=true)
STILL calls update_bdry1(wav_coeff), then mask_tol_boundary_vars reads Domain
coefficients for aliases. validate_and_install_domain_patch_decisions also reads
Domain coefficients for uncatalogued scaffold patches. Native transition must
cover BOTH dependencies, not just source block packing. Preserve two coefficient
versions: scaffold decisions precede boundary refresh; boundary mask decisions
follow it. A full-range native transition workspace/alias plan must be prepared
while old blocks still exist, so boundary alias exchange can run after mask-seed
clears block state. Exchange uses native workspace values, not stale block routes.
Then update transition coefficient records for new source construction. Refine
happens after mask_active, so old integer alias plan remains valid during this
boundary exchange. Force new plan generation after topology/build/migration.
Do NOT simply gather old block coefficients after constructing this transition:
canonical compressed images and consumer aliases have distinct versions. Soil
workspace holds its exact final consumer image; retain that image for transition.
Existing draft transition/provider files alone do not implement these requirements.

Qualification expanded to TEN cases: two added injected readback modes
`owner-readback` and `owner-readback-41` (variant native) call check_gather203 on
canonical soil coefficients immediately after install, and on final soil physical
values after scatter. This proves actual owner storage, not merely workspace or
compatibility publication. No solver-source edits needed; injection helpers bind
original and injected hashes. Reviewer now requires10 cases.
IMPORTANT: running qualifier40381 loaded the OLD8-case list before this helper
edit. After it finishes, run the two extra reject_invalid.py modes sequentially
MANUALLY into .stage243b-validation/reject-owner-readback and
reject-owner-readback-41, then run review_soil.py. The saved qualifier.py now has
all10 for future reproductions. Keep main/source frozen until review succeeds.

Fresh243b full4/full41 PASS; epoch and range fault rejection PASS. Alias fault and
full poisoned run follow. Draft source-builder rewrite generated by
/tmp/draft244_build.py into /tmp/parallel_block_build244.f90, with origin SHA in
/tmp/parallel_block_build244-origin.json. It replaces all28 prognostic/coeff reads
(4 subtree copies +20 ranges +4 scalar samples) through the drafted native provider,
adds explicit optional native_transition selection and rejects legacy fallback
while a native transition is open. NOT APPLIED/COMPILED; reference mean/TKE reads
are intentionally still visible and need their own bounded adapter work.

Fresh243b qualifier40381 ENDED NONZERO: original poisoned full run completed but
step states diverged at step10; first9 states and BOTH VTK outputs were exact.
This run is preserved in reject-poison (do not delete). Test backup incorrectly
assumed scalar/vector/wavelet allocated capacities share one nn. Float_Array elts
is allocatable; restoring entire arrays can resize them. Corrected poison helper
now backs up EACH variable/layer/Domain allocation separately, verifies capacities
are unchanged during transform, and restores slices to prevent reallocation.
This is a test-harness diagnosis to verify, NOT yet proof the solver is clean.
Fresh corrected run goes to reject-poison-exact-capacity. Reviewer/qualifier updated.
Current active session53529 runs /tmp/finish243b.py (/tmp/finish243b.log): positive
owner-readback4, owner-readback41, then corrected full poison, then final reviewer.
No solver source edits; all fresh candidate hashes remain frozen. Fresh243b normal
full4/full41, control, exact oracle, epoch/range/alias faults already PASS.
Do NOT accept/commit until these remaining cases and reviewer actually pass.


## Stage244 — active, not accepted

Stage243 commits574fcbed/b636e2f5 passed ten cases; tracked manifest now binds
immutable detailed files by hash instead of duplicating millions of counters.

Native adaptation implementation added under ADAPT_NATIVE=on/check. Captures
post-physics canonical values and consumer alias versions before soil; appends
soil native result; compiles full-range temporary transition workspace. Scaffold
mask uses pre-refresh coefficients, boundary mask refreshes native coefficients.
All28 build-source prognostic/WAV reads use selected provider; old-node gaps fail
closed and new-node defaults match compressible unsplit extend_Domain. Analytic
reference initialization still uses Domain means/topography. Post-build fullrange
compression/inverse reuse qualified geometric soil kernels. Native retained
boundary state replaces Domain boundary imports. Numerical transition storage
is released after reconstruction/lifecycle close. Compatibility publication still
serves next-step alias import/output/remap consumers: migration NOT complete.

attempt1 .stage244-validation compiled and ALL4 exact adaptation comparisons
passed, then failed old ghost pre-install equality (ghosts legitimately stale).
attempt2 adds optional verify_previous; native retainer skips only pre-exchange
comparison, preserving poison/install/readback tests. Lifecycle closes after
native reconstruction. Frozen attempt2 compiled. Active run47205 exact one-step
oracle writes /tmp/run244b.log and .stage244-validation/oracle-attempt2.
Source edited AFTER attempt2 freeze: lnorms now selects native solution provider
while transition open, needed for dynamic thresholds; Make dependency added.
Thus fresh candidate required before final qualification. Off/default paths remain.
No stage244 commit or acceptance yet; prior frozen candidates untouched.

Stage244 attempt2 exact one-step PASSED including all4 native/legacy adaptation
comparisons and post-install ghost validation. /tmp/run244b.log completed.
First final candidate (.stage244-validation/candidate) native and production
compiled, but control failed unused private capture routines under -Werror.
Fixed by #ifdef ADAPT_NATIVE around the two capture routines; no arithmetic change.
Fresh candidate2 building control/native/production sequentially, session89916,
/tmp/build244candidate2.log. Control/native binaries built so far; production
pending. Helpers qualification/reviewer point to candidate2 (verify before launch).
All Stage244 source now frozen awaiting9 cases: control4, oracle4, fullprod4/41,
staleepoch, stalegeneration, full Domain poison, ownerreadback4/41. Helpers parse;
15-kernel audit PASS with new Stage244 artifact. Candidate2 reviewer not yet run.

Next lifecycle audit saved in /tmp/native245-plan.txt. Includes pre-RK boundary,
remap pre/post imports and polar numerical update, initial alias capture, cpt_dt
(which is currently entirely Domain numeric), mass diagnostics, output alias
mutation and checkpoint coarserestriction + immediate restart. Native full-range
workspace after244 reconstruct can feed these; temporary alias pack bridges next
atmospheric plan rebuild. Need explicit bootstrap fullrange adapter after restart.
No stage245 source edits applied. Do not claim implementation complete.

Stage244 candidate2 all THREE builds PASSED (session89916 completed). Qualification
session71093 ENDED NONZERO: control4/oracle4 pass; fullproduction4 first TWO states
exact, then fails closed atstep3: missing old transition operand rank1/localDomain6/
node2061. Failed outputs preserved. No full41/fault/poison cases launched yet.
Diagnostic injection running session74251 via /tmp/diagnose244.py, outputs
.stage244-validation/missing-node-diagnostic and /tmp/diagnose244.log. It augments
missing-node error with native key, family/variable/layer, mask/level/old/newlength,
and legacy value/defaultmass; frozen/live solver source UNCHANGED. Await diagnosis.
Likely canonical stored node not registered as an operator dependency: map_patch
sets nodes%key but does NOT node_slot every interior, whereas sourcebuilder copies
whole patches. collect_boundaries already registers all old compact boundary nodes.
If diagnostic keypositive, extend ADAPT_NATIVE metadata closure to all canonical
old patch nodes before plan construction, and gather exact native owners; NEVER
fallback to Domain/default for a missing old operand. If keyzero, investigate
boundary lifetime separately. First2 states verified exact with analyze.py.

Unapplied Stage245 draft /tmp/draft245_inverse.py generates
/tmp/parallel_block_inverse245.f90 + originSHA; API text /tmp/native245-lifecycle.f90.inc.
Adds fullrange bootstrap/alias bridge, lifecycle getters, optional vertical range
in alias transfers. NOT APPLIED/COMPILED; needs main/remap/diagnostics integration
and extensive review per /tmp/native245-plan.txt. Generator continuation parser was
corrected, but future source244 edits mean regenerate/check origin before reuse.
Do not resume245 until244 accepted. No other solver/build jobs active beyond74251.

Diagnostic74251 completed: missingnode key[0,0,0,0], familySOL/scalarmass/layer-10,
oldlength2573 ->2629, mask0/level6. Followup diagnostic2 failed compile due test
insertion overlongimport; corrected diagnostic3 session38445 completed: boundary
index79/start2061/side0/neigh0; ALL12 native slots absent. Thus retired temporary
boundary storage, NOT canonical. No arbitrary default fallback was introduced.

Added ADAPT_NATIVE integer closure registration of COMPLETE canonical patches and
all allocated boundary storage units (including retired side0). Native aliases
retain exact existing values through current input adapter; later sourcebuilder
can copy full boundary units. Current compact metadata collects only referenced
boundaries, insufficient for later offset closure; broader explicitaliasstorage
is necessary. Initial alias adapter removal remains Stage245 work.

Fresh candidate3 production build PASSED; active FULL4 run session53057,
/tmp/run244candidate3.log, .stage244-validation/production3-4. First3 states exact,
including formerly failing transition. Wait complete12steps/output/checkpoint gate.
Then build native/control into SAME candidate3 (sources frozen; helper supports
adding variants), run remaining qualification. Helpers now pointcandidate3 and
control3-4/oracle3-4/production3-{4,41}; qualify skips existing production3-4 only
if completed status/12steps/binaryhash pass. Reviewer requires all9cases. Stage244
NOT accepted/committed. No additional solver/build jobs currently active.

New audit in /tmp/native245-plan.txt records a remaining RK vector coefficient
carrier: native vectorflatstage -> DomainWAV -> import_domain_field_family ->native
compression -> Domainwriteback ->inversealiasseed. Must remove this too before
mirror cleanup. See precise callsite plan; NOT implemented. Don't claim complete.


Candidate3 FULL4 passed12 exactstates,2remaps,checkpoint numerical gate,62exactVTK.
Added native/control builds to same frozen tree (session76005 completedPASS).
Active qualifier85522 /tmp/qualify244candidate3.log: control3-4 PASS; oracle3-4
running, then skipalreadyacceptedproduction3-4, full41 and five injectedchecks.
Do not run another solver/build concurrently. Reviewer still notrun/accepted.

Stage245 additional UNAPPLIED generator /tmp/draft245_integration.py outputs
/tmp/main245.f90, /tmp/remap245.f90, /tmp/parallel_block_mpi245.f90 with originSHAs.
Inverse generator regenerated againstcandidate3 source, now has outputaliasadapter
API including keyzero plus canonical alias-destination slots. Nonecompiled/applied.
Known review needed: optional remap oracle does extra boundaryrefresh909 before
writeback and910 afterward; draft postremap full comparison needs exactphase
alignment. No Domain numericalmirror removal or RKcoeffcarrier changes in245draft.

Full41 has finished12steps through0.3308day; helper validation pending at latest
inspection (qualifier85522 still active). Added optional --oracles to Stage244
run_adaptation helper so we can ALSO test expanded alias coverage with the complete
independent RK/physics oracles, not only adaptation/soil checks. Reason: registration
occurs before original_inverse_nodes, which participates in RK low-scaffold update
and comparison coverage. This is a specific new closure concern, not broad retesting.
After85522 finishes, add oraclevariant to candidate3 and run rk-oracle3-4 with
--oracles --checks before acceptance; update reviewer to requirethis tenthcase.
Don't build concurrently with the running solver. Current saved qualifier still
has original9cases; extraoracle manualfollowup notyetexecuted.

Candidate3 FULL41 PASS:12exactstates,2remaps,62VTK/checkpoint validation (qualifier
printedPASS). Epoch/generation faulttests PASS. Full Domainpoison running under
85522; ownerreadback4/41 follow. Final reviewer updated to require TENcases and
fourbuilds, including rk-oracle3-4 (independent fullRK+physics oracle;16exact RKalias
markers). Saved qualifier now includes this extraoracle forfuture reproduction,
but RUNNING85522 loaded old9-case list: afteritfinishes manually buildvariantoracle
and run rk-oracle3-4 with --oracles --checks, thenreview. Mustn't accept earlier.

Stage244 accepted: all ten cases and final source-bound reviewer PASS, including full independent RK oracle. See handoff244 and results manifest. Continuing lifecycle and storage cleanup; no bbserv run.

Stage245 lifecycle PASS: frozen candidate2 production full4 trajectory (12 exact states, two remaps, 62 exact VTK members and checkpoint gate), checked native one-step4 PASS. Source-bound reviewer passed. Initial candidate compile failed missing BLOCK_PAYLOAD_SOL import; preserved. Native lifecycle retains full-range aliases across steps, adapts to atmosphere-only RK input, and supplies timestep/mass diagnostics and remap poles/boundaries. Remaining: coefficient transport and persistent mirrors. Final 41-rank qualification will cover integrated cleanup.

Stage246 direct coefficient transport passed full4: 12 exact states, two remaps, 62 exact VTK members, checkpoint gate. Vector producer staging transfers directly to owners; compression no longer imports Domain WAV, inverse aliases start from native coefficient epoch. Legacy coefficient publication is oracle-only. Source-bound manifest passed. Independent oracle and empty-geometry-rank checks deferred to final integrated storage qualification.

Stage247 development (not accepted): candidate1 removed active vector staging incorrectly, failed first RK transform. Candidate2 keeps necessary native vector producer workspace, removes obsolete scalar/snapshot mirrors, and supplies scaffold soil from native saved aliases; bounds-checked run reached final RK but revealed physics-ready state was coupled to compatibility publication. Candidate3 separates native consumer activation, adds native lifecycle norm input and selector guards. Failed candidates preserved. Current live cleanup not yet qualified.

Candidate3 passed RK, physics, soil, adaptation then found remaining Domain pole mass read in minimum-mass remap trigger (coarse patch branch was already native). Candidate4 routes pole mass through same native lifecycle accessor. No tolerance changes.

Candidate4 completed first exact step with storage absent; next initial tendency still called legacy update_bdry1067 despite native lifecycle boundary979. Symbolicated crash confirms that call. Candidate5 skips it on native lifecycle, preserving independent oracle. Debug build reuses unchanged candidate4 objects and rebuilds affected dependencies; final production/oracle builds will be fresh.

Candidate5 full checked trajectory PASS (12 exact states, two remaps, 62 exact VTK members, checkpoint gate, runtime storage absence). Optimized build passed; strict separate oracle build rejected unused native_output_allocations outside MIRROR_FREE. Candidate6 guards that diagnostic declaration; arithmetic unchanged. Final integrated qualification uses fresh candidate6 binaries.

Stage247 accepted: frozen candidate6, final6 cases. Full checked4 and optimized4/41 each passed 12 exact states, two remaps, 62 exact VTK members and the existing checkpoint numerical gate. Independent RK/physics oracle passed (16 exact alias markers); off-path control and lifecycle/coefficient/mirror fault rejection passed. Runtime absence checks cover rank40 with no Domain geometry but active compute ownership. Fifteen arithmetic kernels and final source-bound eight-case review PASS. Native lifecycle/direct coefficients and bounded output adapters complete; persistent SOL/WAV numerical views, obsolete RK snapshots and persistent scalar staging removed for qualified scope. Active native vector arithmetic, alias/stencil and reference buffers remain. Main unchanged. Preparing sealed bbserv bundle; no remote run.
