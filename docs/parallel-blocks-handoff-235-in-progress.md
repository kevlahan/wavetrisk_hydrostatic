# In-progress continuation — full tendency migration remains the active task

User's latest request is “Continue to complete full tendency migration.” Do not
stop at this handoff. No commits or branch changes. Working branch remains
codex/performance-investigation, HEAD 4ea640128a71a857f9b8c6b9411dc34fc0aeb641;
main must remain 6a245db7f1632b0205c2bebe95e5f0cd535c2582. Stages231 onward
uncommitted; hundreds of unrelated historical untracked files; never stage all.
No subagents authorized. No remote cluster jobs. This continuation adds real
owner execution; it does not demonstrate a performance gain.

Read Stage234 handoff for earlier work. Latest root source has further changes:

- `parallel_block_resident_operands`: payload7, first5 unchanged, raw mass/temp
  rows6/7 fetched directly, not reconstructed by subtracting means.
- `parallel_block_mass`: pure `mass_restriction_value` extracted exactly;
  old executor calls it, same counters/writes.
- `parallel_block_program`: distributed immutable numerical operation graph.
  Geometry compiler refs5 = owner/origin/local-operation/signed-row/depth.
  Owners retain results; MPI routes operations, geometry, deduplicated input
  reads by depth. Supports resident inputs, gradient/velocity flux, divergence,
  restriction, scalar physics+advective edge flux, mass divergence, zero,
  circulation/curl, velocity physics/direct/restricted source. `program_value`
  reads a current owner result; external callback codes >=100 supply explicit
  other-producer operands. Missing definitions propagate as absent refs (not
  numerical zeros), live read compilation rejects them. Unused compatibility
  boundary writes are poisoned, and no missing operand reads Domain numerical
  scratch.
- `parallel_block_diffusion`: exact topology/order compiler for divu, scalar
  Laplacians, Laplacian(divu), scalar physics + advective flux; also owner mass
  restriction/divergence when BLOCK_MASS_OWNER enabled. Geometry-only definition
  maps, actual numerical execution and storage on final owners. Signed remote
  snapshots then all-level ordered local boundary copies then remote installs.
  Some scratch slots are explicitly initialized to zero and never written in
  ANY level. Topology coverage seeds zeros only in those slots; written slots
  with missing producers remain undefined. Restriction dependencies split by
  selected edge; scalar flux ops split by edge. Native mass compatibility
  delivery remains for RK/remaining consumers. APIs expose mass/diffusion
  producer refs/results for new velocity-source consumers.
- `ops`: scalar physics and advective fluxes use owner results when active.
  Independent old callback/arithmetic compares actual consumed edges. Important:
  native temperature path omits some rho loads; oracle must not evaluate
  unconsumed temperature edges, or checked build traps on signaling NaNs.
- `multi_level`: dispatches diffusion/mass routes; disables old execution in
  production, retains it for check. Resets diffusion/mass active flags at end.
- New `parallel_block_primitive_registry`: integer identities only for latest
  PV/vorticity writes. key7 = kind,owner,origin,relative-packet-or-operation,row,
  level,depth. kind1 regular primitive,2 pentagon,3 new velocity graph,4 explicit
  initialized zero. Tracks all primitive write destinations per layout; seeds
  zeros only outside that coverage. Written-but-not-produced slots remain
  undefined. Reset per physical layer at highest spatial level.
- `parallel_block_primitives` registers regular/SW PV/vorticity owner addresses;
  exports primitive_owner_value (key5 owner/origin/packet/row/level, checks
  packet/row/layer). `parallel_block_pentagons` applies metadata overwrites,
  exports pentagon_owner_value with equivalent checks; per-level layer stamp.
- New `parallel_block_definition_exchange`: metadata7 signed edge exchange;
  optional coverage mode marks receive/copy destinations as written regardless
  of source definition. Avoid passing noncontiguous derived-component arrays
  with allocatable members: GNU Fortran produced temporary/double-free when
  earlier code passed map%flux to a procedure. Current code accesses explicit
  standalone arrays or map fields through helpers.
- New `parallel_block_velocity_owner`: per-level distributed programs for curl,
  circulation, pentagon corrections, diffusion physics, and EXACT ordered
  direct/restricted source action tape. Sources consume actual owner mass and
  primitive results. Dependency masks for Qperp obtained from its existing pure
  arithmetic, not all unused entries of the rectangular stencil. Previous
  spatial-level source/curl results are explicit external operands. Independent
  rot/source checks in check mode. Production no longer returns intermediate
  velocity source to geometry native workspace after latest gradient hookup.
- `parallel_block_gradients` latest ROOT source now fetches source via owner
  addresses/results; legacy native source capture only when velocity owner off.
  Existing direct retained gradient block storage remains.

Make selectors: BLOCK_COLUMNS, BLOCK_PRIMITIVES, BLOCK_DIFFUSION,
BLOCK_MASS_OWNER, BLOCK_VELOCITY_OWNER all off default. Each later selector
requires previous on/check. BLOCK_DIFFUSION currently requires TEST_CASE=climate
because the migrated physics is that callback. Extra dependencies were added;
new owner source dependency pulls utils/comm earlier, so explicit comm/comm_mpi
and utils deps are present. Use fresh build dirs when types/options change.

## Completed local evidence (frozen trees, NOT latest root source)

All tests use .stage230-validation/20260929/fixture and same-rank Stage230 refs.
Helpers currently reused from primitives_234/run_primitives.py. Dedicated new
helpers/review still need writing; no new campaign archive exists.

`.stage235-validation/candidate` is frozen diffusion/scalar owner implementation
WITHOUT mass owner or newer primitive registry/source changes.
- bin-235-check3, probe8-4: one evolved checked step PASS.
- bin-235-native-check, diffusion-native-4: no-oracle path with runtime checks,
  absent scratch poisoned: one step PASS.
- diffusion-full-4: checked diffusion full12steps/2remaps/restart PASS, exact
  checkpoint fields, step states, 62 VTK members vs Stage230 production-4.
  results.json and primitives-results.json exist.
Earlier probes1..7 failed during implementation; keep evidence, do not count
as qualifications. Principal findings: unused boundary definitions; per-edge
liveness; genuinely initialized-never-written zeros. probe7 poison reached a
live scalar divergence until explicit initialization was represented.

`.stage235-validation/mass-candidate` is frozen diffusion + mass owner, WITHOUT
new registry/source/gradient hookup.
- bin-check, mass-probe1-4: one-step independent mass checks PASS.
- bin-native-check, mass-native-4: old mass execution skipped, runtime checks PASS.
- bin-production, mass-full-41: full12steps/2remaps/restart PASS, exact checkpoint,
  step states,62VTK members vs Stage230 production-41. Includes rank40 with no
  Domains. results.json and primitives-results.json exist.
These full tests both completed and exec sessions were polled successfully.

## Currently working / next actions

`.stage235-validation/source-candidate` is a mutable development tree, overlaid
from root. `bin-check` initial source connection built, but source-probe1-4
failed because primitive registry had no explicit initialized-zero definitions.
`bin-check2` included latest direct gradient hookup but preceded zero-definition
fix; built but NOT run (same known issue).
Newest root changes declare all primitive write destinations across all levels
and seed known untouched PV/vort slots with explicit zero. Overlaid into
source-candidate. A fresh build-check3/bin-check3 checked build is RUNNING:
command make -C .stage235-validation/source-candidate -j1 TEST_CASE=climate
PARAM=param_J4 DEBUG=check FP_CONTRACT=off BLOCK_COLUMNS=check
BLOCK_PRIMITIVES=check BLOCK_DIFFUSION=check BLOCK_MASS_OWNER=check
BLOCK_VELOCITY_OWNER=check BUILD_DIR=build-check3 BIN_DIR=bin-check3
log `.stage235-validation/build-source-check3.log`; exec session56255.
No other builds/tests should still be running.

Finish this build, run source-probe2-4 with bin-check3, resolve any remaining
producer/order issues WITHOUT Domain fallback. Then no-oracle checked and
optimized full 4/41 adaptive/restart tests, runtime oracles, fault injection,
feature-off control, source-bound review. Preserve already frozen prior trees.
Local MPI needs exec require_escalated (approved throughout this session).

The user asked for FULL migration, not just another partial stage. Remaining
work after validating source path: direct temperature/mass retention/capture
cleanup; stage-one/restart/post-physics/final-transform alias epochs; reference
means/surface input ownership; remove no-longer-needed Domain pathways/storage.
Current columns arm only stage>1. Native inverse `complete_inverse_operand_epoch`
currently marks provisional only; stage1 cannot simply flip that flag: final
inverse has distinct consumer_alias_value, and adaptation/physics changes state
before next stage1. Existing grid-change publication selects saved consumer
aliases after native sync. Preserve this lifetime exactly; no whole Domain-shaped
solution snapshot and no stale aliases. Need inspect time_integr/multi_level and
parallel_block_mpi grid-change/import lifecycle for correct owner producer bridge.
Do not claim full migration, Domain removal, or performance improvement yet.


## Continuation update: initial stage and compatibility cleanup implemented

Root source now includes BLOCK_TENDENCY_INITIAL=on, initial physical alias import
(inverse prepare consumer_only), require_tendency_operands and first-stage arming.
Stage-one binding occurs after update_bdry(q), before column startup. Canonical
interiors cannot be read through initial compatibility operand getters.

New owner_tendency_without_reference predicate (columns module) is true only for
complete owner route with all compile/runtime oracles disabled. It suppresses
PV/vort/KE/Bernoulli installs in primitives/pentagons, pressure/geopotential/Exner
installs in diagnostics, numerical Bernoulli/Exner restriction, and diffusion
result delivery/install. Integer thermodynamic producer propagation remains.
Final gradient return MUST remain: native_velocity_rk in time_integr consumes it.
It is explicitly documented as RK return, not a tendency operand.

.stage235-validation/qualified-candidate freezes latest root sources and Makefile.
Builds native(check runtime, all selectors on), production, oracle(check selectors),
and control(all off) done. build_tendency helper records source/binary hashes.
Do not edit this frozen tree while tests are active.

New helpers test/parallel_block_profile/tendency_235:
- build_tendency.py: frozen variants and manifest build-235.json.
- run_tendency.py: all-stage validation + explicit boundary descriptions.
- reject_invalid.py: isolated missing producer, stale epoch, corrupt source;
  scratch variant poisons removed pressure/geopotential/PV/vort/KE/B/Exner/divu/
  scalarLap/lapdiv each layer and checks they remain poisoned.
- review_tendency.py: source-bound acceptance, final durable report.
Draft documentation docs/parallel-blocks-stage235-tendency-migration.md.

All-stage tendency KERNEL migration is implemented for climate, but full Domain
storage removal is NOT claimed. Topology/geometry, initial aliases and means,
scalar transport and native RK closure remain explicit boundaries. Default off.
No bbserv jobs or commit; main unchanged.

Current sessions (poll, then launch required remaining tests):
- 68337 old source-full-41, ~10/12 steps as of update.
- 5673 initial-probe1-4 (checked owner initial stage, older source before cleanup).
- 48007 qualified-production-4 full12 steps, final candidate.
- 25934 scratch-4 final native injected object, one step.
- 93417 qualified-all-oracles-4 final oracle, one step.
- 85862 sequential reject-missing / reject-epoch / reject-source final injections.
- 59017 control-4 final control one step.
- 4128 final manifest adds native+production to existing oracle/control manifest.
Source-native-4 PASSED (60190 consumed). Source-probe2-4 earlier passed.
Native/prod final build sessions83784/96314 consumed. Oracle/control98894 consumed.

Still required: qualified-production-41 full with final all-stage candidate (wait
old 41 run finishes to avoid excessive resource competition); source-bound review
report; resolve any numerical failures. The review no longer requires separate
qualified-oracle-4, since all-oracles one-step covers it. Do not stop at handoff.


### Latest RK return correction (qualification still running)
- initial-probe1-4 PASSED; source-full-41 PASSED full12/2remaps/62VTK.
- control-4 and reject missing/epoch/source PASSED on pre-fix source.
- qualified-all-oracles-4 FAILED native mass RK exact comparison.
- Isolated /tmp/probe235_mass.py -> mass-boundary-probe identified undefined
  native_mass output slots NaN versus reference +0, mask8/12. All defined
  numerical mass operations had already passed independent checks.
- Root diffusion install now initializes the RK return tendency buffer to +0
  (original zero_float), then fills defined outputs. Graph missing refs stay
  missing; flux scratch stays poisoned. This is output initialization, not a
  missing-operand fallback.
- Old qualified-candidate renamed pre-rk-fix-candidate. NEW qualified-candidate
  is freezing/rebuilding root with this fix: session52429, log build-rk-fix.log,
  variants oracle/native/production/control. Wait before final tests.
- Production-4 interrupted deliberately for heavy paging; preserved at
  production-4-interrupted-memory-pressure. No acceptance claim from it.
- Still running old pre-fix scratch-4 session25934 (nearly done), and old
  qualified-production-41 session59520 (~6/12). Preserve under pre-rk-fix names
  after completion, then rerun final code. Old fault/cntrl outputs need similarly
  renamed before reusing required final review names.
- Run final tests mostly serially; concurrent 4-rank generic owner graphs cause
  heavy macOS paging. Do not repeat that scheduling mistake.
- Remove/rewrite any stale report claims; source-bound review has not passed.


### Tests scheduling/status after RK fix build
- New qualified-candidate ALL FOUR builds finished (52429 consumed), source-bound
  manifest build-235.json exists and matches root (only new correction is native
  mass RK return initialized +0; numerical graph refs unchanged).
- Final corrected qualified-all-oracles-4 RUNNING session29395. It PASSED the
  previous failing first-stage RK handoff, now progressing stage2; wait for full
  one-step completion before acceptance.
- Pre-fix final41 full all-stage trajectory PASSED12/2remaps/62VTK; its directory
  renamed pre-rk-fix-qualified-production-41 (59520 consumed).
- Pre-fix scratch PASSED (25934 consumed), renamed pre-rk-fix-scratch-4.
- Other old results renamed pre-rk-fix-qualified-all-oracles-4 (failed),
  pre-rk-fix-control-4, pre-rk-fix-reject-{missing,epoch,source}.
- Old source-full41 and initial-probe1-4 both passed and consumed.
- No other MPI runs should remain except final corrected all-oracles session29395.
- NEXT run final suite MOSTLY SERIALLY: production41 full, production4 full,
  scratch4, reject missing/epoch/source, control4. Same helpers and fixture.
  Full reference paths .stage230-validation/20260929/production-{4,41}.
- reject_invalid scratch now additionally backs up one physical layer of q,
  poisons its mass/temp/velocity arrays AFTER explicit alias binding, asserts
  no compatibility write, restores for downstream RK. This test-only backup is
  injected into an isolated object, not the production source. Not yet run.
- Draft docs updated; final report cannot be written until review_tendency passes.


### Current acceptance suite (latest authoritative status)
Corrected final all-runtime-oracles4 PASSED; session29395 consumed. Root source
and final frozen qualified-candidate match. Build manifest now records both link
identities and resolved source contents. No further production source changes.

ONE SERIAL SUITE is running in exec session27184 (MPI escalation approved):
1. final qualified-production-41 full12 (currently3 completed steps).
2. final qualified-production-4 full12.
3. final scratch-4: removed scratch AND primary q fields poisoned after binding,
   restored only for downstream RK; asserts poisoned fields remain untouched.
4. final reject-missing, reject-epoch, reject-source.
5. final control-4.
6. review_tendency.py generates docs/parallel-blocks-stage235-local-results.json.
The Python driver stops at any failure; poll27184 and read current log. Do NOT
start additional large MPI tests concurrently (paging made earlier tests slow).
All old/pre-fix source/full results are preserved separately and do not count
as final-code qualification. First all-runtime-oracle failure was fixed and
corrected run passed. Keep main unchanged. After suite passes, update final docs
and handoff, verify git diff --check/main/source hashes, and report outcome with
explicit remaining input/RK/diagnostic Domain boundaries and no speed claim.

Final corrected qualified-production-41 PASSED12/2remaps/exactfields/states/62VTK.
Serial suite27184 is now running final qualified-production-4 full (1/12 steps).
No production source changes since the RK +0 fix. Only test helper improvements:
source_contents resolves symlinks; epoch injection now calls owner_diffusion_value
with k+1 without changing producer guard; scratch also poisons primary q fields.
All helper syntax compiled. Keep running suite; it will proceed automatically to
scratch, missing/epoch/source, control, and source-bound review. Four-rank full
coverage remains useful because it exercises ordered local copies across many
Domains per rank after grid changes (41 ranks mostly exercises remote copies).

### Lookup correction after final four-rank timeout
The final four-rank trajectory reached 11 complete steps but hit the 2400s
watchdog in the final step; preserved as qualified-production-4-timeout. No
acceptance claim. All MPI processes exited. Session27184 ended.
Remaining pre-hash-fix scratch/missing/epoch/source/control checks run serially
in session50785. Do not rename their qualified-candidate tree while running.
A bounded five-second sample of rank1 in scratch-4 found 2108/3596 main-thread
samples in producer locate, including1963 via diffusion external values. Raw
samples lookup-sample-r{0,1}.txt are preserved in .stage235-validation.
Hash used translated contiguous IDs, causing long linear-probe clusters. Root
parallel_block_program now mixes IDs with two bounded 32-bit rounds using int64
products (no overflow). Numerical operations/ordering untouched. New four-variant
frozen build hash-candidate running session93771; log build-hash-fix.log. After
old suite finishes, preserve old qualified results with pre-hash-fix prefix, swap
new hash-candidate to qualified-candidate, rerun full qualification serially.

### Final hash-corrected acceptance suite
Old scratch stopped deliberately after collecting bounded CPU samples; preserved
at pre-hash-fix-scratch-interrupted (not pass evidence). Old suite50785 exited.
New hash-corrected four builds PASSED (93771 consumed). Hash-corrected scratch
poison PASSED in47.584s (71766 consumed), including primary q arrays and removed
Domain scratch kept inaccessible. Old tree/results preserved with pre-hash-fix
prefix. hash-candidate renamed qualified-candidate; production root source matches.
ONE SERIAL acceptance suite29625 now runs final production4 full12, production41
full12, all-oracles4 one step, missing/epoch/source faults, control4, then review.
The final scratch test is already done and must not be rerun. No other solver
runs remain. Preserve timeout evidence; final report requires new completed runs.

Final hash-corrected production4 and production41 BOTH PASSED all12steps,
2remaps, checkpoint fields, step states and62VTK members exactly. Rank40 owns
blocks with zeroDomains and executed12 initial input epochs plus48 owner programs.
Suite29625 now finishes all-runtime-oracles4, then missing/epoch/source and
control checks and source-bound review. Final scratch already PASSED on this
source; do not rerun. Final handoff drafted in parallel-blocks-handoff-235.md.
review_tendency additionally binds each fault test original source hash to the
final frozen source. Helpers compile; root/frozen hashes and diff check passed.

### COMPLETE — final authoritative handoff
All final acceptance tests PASSED; suite29625 consumed. Source-bound report
written to docs/parallel-blocks-stage235-local-results.json. See
docs/parallel-blocks-handoff-235.md for final status; earlier session IDs and
in-progress qualification notes in this rolling file are historical.
Full tendency KERNEL migration validated; Domain storage removal is not complete,
selectors remain off, no performance claim, main unchanged, no commits.
