# WAVETRISK parallel-blocks optimization handoff — 11 September 2026

## Start here

Continue in the existing `/Users/kevlahan/wavetrisk_hydrostatic` checkout on
`parallel-blocks`. This handoff is committed together with Stages 182/183.
The commit immediately before it is `78ba3574` (Stage 181 allocation census).
Resolve this checkpoint with `git log`; do not assume a later HEAD still names
this stage. Do not reset, re-extract an old ZIP, or replace the working tree.

The next task is **acceptance and measurement of Stage 183**, then a bounded
column-stencil/accessor optimization. Do not repeat the completed rounding
investigation or completed local runs without a concrete new reason. No MPI job
is left running by this handoff. No push or cluster transfer is included in
this commit. The user normally transfers source, builds/runs on the cluster,
and supplies logs. The user has requested creation of the new project chat.

Read in this order:

1. This file: current state, constraints, next work and pitfalls.
2. `parallel-blocks-shared-scalar-geometry-183.md` and its `-results.json`:
   implementation, source identities, oracle settings and validation evidence.
3. `parallel-blocks-mixed-precision-diagnosis-182.md` and its `-results.json`:
   completed numerical diagnosis and its precise limits.
4. Only as needed: `parallel-blocks-allocation-census-181.md`,
   `parallel-blocks-local-protocol-180.md`,
   `parallel-blocks-restart-corrections-179d.md`, and
   `parallel-blocks-local-profiling-179c.md`.

Older stage documents are historical snapshots, not current authorization or
current acceptance status. Some say “uncommitted” or “next” at their original
writing date. This handoff supersedes those status statements, not their data.

## Project objective and decisions already made

- Originally parallel-blocks was about 4x slower than legacy on 83 tasks;
  later user observations still indicated about a 2.5x per-step deficit.
  These are historical observations, not a controlled current Stage 183 ratio.
- Remove dual representations **from computation**. Production must not run a
  legacy Domain producer and then repeat its work in blocks, or repeatedly
  write native values to Domain solely to invoke a callback and re-import them.
- Domain may remain for metadata, geometry ownership, initialization, I/O and
  other interfaces where it does not materially burden the timestep. Do not
  interpret this as permission to retain expensive compatibility everywhere.
- Extra reference computation belongs only to explicitly selected testing
  oracles. Eventually oracle machinery can be retired, after independent
  confidence is established; do not remove it prematurely to make tests pass.
- The user authorized **native geometry-owner velocity execution first** in
  Stage 175. Native mass also runs on geometry owners. Final-owner placement
  remains deferred, not completed or proven necessary for speed.
- Reaching within roughly 20% of legacy speed is an aspiration, not an
  established attainable bound. There is no evidence that the remaining gap is
  irreducible block overhead. Quantify necessary work versus removable costs.
- Prefer substantial complete producer-to-consumer changes with observable
  work removal, not many tiny traffic-only stages with no measured benefit.
- Use the column structure: physical vertical layers share horizontal geometry.
  Geometry sharing must still respect soil/surface, masks, poles and boundary
  roles. Geometry is not equivalent to every field's validity/active state.
- Authoritative `main` must remain unchanged. The user permits numerical-neutral
  instrumentation or causal experiments in clearly identified isolated legacy
  copies; never deploy those experiments as production or commit them to main.

## Source architecture and what is already replaced

| Area | Current implementation / entry point |
| --- | --- |
| Block state and kernels | `src/parallel_block.f90`; catalog construction in `parallel_block_build.f90` |
| RK orchestration | `src/time_integr.f90`, particularly `prepare_block_multistage_tendency`, `RK_sub_step`, native wavelet calls |
| Shared producer, geometry-owner native programs | `src/multi_level.f90`, `trend_ml`, native mass/velocity plan and execution routines; primitive operations in `src/ops.f90` |
| Transport, ownership, scalar execution and oracle orchestration | `src/parallel_block_mpi.f90` |
| Native inverse workspaces and kernels | `src/parallel_block_inverse.f90` |
| Native velocity kernels/programs/RK values | `src/parallel_block_velocity.f90` |
| Native mass restriction/divergence/boundary/RK | `src/parallel_block_mass.f90` with orchestration in `multi_level.f90` |
| New compact scalar records | `src/parallel_block_scalar_storage.f90` |
| Profiling | `src/parallel_block_profile.f90`, instrumentation in producer/transport callers |

Stages 168–172 established temperature production and scalar transport/native
inverse-boundary cuts; read the corresponding stage notes for the exact scope.
Stage 173 removed repeated whole-column thermodynamic reconstruction at edge
lookups. Stage 175 replaced Domain velocity source/gradient/tendency consumers
with native geometry-owner execution and publication. Stage 176 removed full
producer-side scalar records/repacking in production. Stage 177 replaced the
mass restriction, boundary, divergence and RK dependency chain together, without
also replaying completed mass work on final owners. Stage 178 compiled inverse
transport/stencil addresses and removed some packing/synchronization work.

These changes do **not** mean all Domain production work is gone. Shared
PV/pressure/KE/Bernoulli/Exner and external physics/diffusion interfaces, solution
compatibility at remaining consumers, adaptation/remap/interfaces, and ownership
publication still require an explicit current call-chain audit before claiming
their removal. Some retained primitive computation is necessary work executed
once, not a duplicate oracle. Separate this from avoidable representation
transfers. Do not infer costs from a routine name containing “Domain”.

Production profiles must continue to show zero legacy velocity-source/gradient
calls and zero removed mass-flux compatibility work. Stage 176 producer full
oracle-record counts must remain zero. Stage 183 production scalar `%full`
allocations are zero. Oracle mode intentionally has additional reference work.

## Stage 183 storage contract

The existing scalar logical record has 50 slots. Production now stores:

- 17 field-dependent values per scalar/vertical/node sample;
- 33 geometry/metadata values per horizontal node, with two classes: physical
  and inactive. The second class preserves local soil/surface zero geometry
  versus transported inactive geometry, without expansion across all fields.
- Oracle mode retains independent full 50-slot records, not aliases of compact
  production storage. Test both modes; an oracle-only pass cannot certify the
  compact production allocation/indexing path.

Field slots are `1:6, 9:11, 15:20, 34:35`; all other slots are shared. The new
storage API preserves the logical addressing contract through scalar/range/
indexed reads and writes. `scalar_extent` is a logical extent, **not bytes**;
use `scalar_capacity` for actual capacity. Patch layout is patch-major, then
field, then 16 nodes. Boundary/ghost layouts are field-major over their nodes.
Canonical shared writes prevent one field's zero initialization from erasing
another's geometry. Allocation/mode/topology changes rebuild affected storage.

Producer geometry is seeded once per patch, remote geometry once per node;
field buffers are initialized separately. Three/four live-field transfers avoid
reconstructing full records. Existing MPI wire formats, numerical formulas,
operation order, ownership and generation lifetimes are unchanged. The logical
access API still has indexing/dispatch overhead: removing that from hot column
loops is a plausible next computational target, not an already achieved result.

Exactly four production files changed in Stage 183: `Makefile` plus
`parallel_block_scalar_storage.f90`, `parallel_block_mpi.f90`,
`parallel_block_profile.f90`. Transfer all four when the remote baseline is
otherwise current. `multi_level.f90`, `time_integr.f90` and `shared.f90` did not
change in this stage. The one-time migration helper is not a production step.

## Numerical acceptance: keep these distinctions explicit

1. **Stage 179d fixed real bugs**, committed at
   `1b01510ae655eee0e2aab1a8d911e55c147173f2`: missing block-root forward wavelets;
   different forward/inverse level ranges; incorrect physics boundary input
   epoch; production polar remap missing unless oracle enabled; stale module
   dependencies. Preserve these fixes and their negative/coverage tests.
2. Stage 180 extended validation over two new checkpoints/reloads and remaps.
   It exposed a smaller independent legacy difference despite passing oracles.
   Printed time/DOF/mass agreement and a previous matching checkpoint were
   insufficient. “Byte correct” must never imply untested future states or
   independent legacy equality.
3. Stage 182 causally diagnosed the new discrepancy. At step 260, Domain 77,
   raw node 123, layer 26, U differs by only `8.0984e-14` in double precision
   but crosses a float32 midpoint; physics output differs by two float32 ULPs.
   A second W crossing occurs at step 262, Domain 8, node 123, layer 28,
   from a double difference of `2.9896e-13`.
4. Isolated legacy replay changing **only those two packed inputs**, once each,
   brings both checkpoints inside the **unchanged** original screen. This
   causally explains errors above the screen on this fixture, not every possible
   discrepancy or a universal numerical bound. No replay is linked into current
   production. The user explicitly authorized proceeding after this diagnosis.
5. Original unmodified legacy versus block checkpoint 6 still fails the fixed
   screen: max atmospheric wavelets velocity `5.9972e-9`, mass `5.7042e-9`,
   mass-weighted temperature `2.7876e-6`. Screen limits are `1e-10`, `1e-10`,
   `1e-7`; soil/topology/thresholds/headers exact; no nonfinite values allowed.
   With the two diagnostic interventions maxima become `1.1006e-11`,
   `5.5524e-12`, `5.1652e-9`. Do not silently waive/relabel the original failure.
6. Stage 183 preserves the existing **block baseline** exactly in every compared
   checkpoint field at checkpoints 5 and 6. This is not byte equality of entire
   compressed files and not new independent legacy equivalence.

The strict legacy timing campaign still stops at its original numerical gate.
Do not quietly disable it. Any timing under the user's diagnosed-rounding
authorization must keep that original result visibly separate, require exact
candidate-versus-block equivalence, and report timing-quality limitations.

## Completed Stage 183 validation — do not repeat by default

- Clean optimized RK4 and checked RK3/RK4 builds succeeded.
- 48 Python/helper tests, including compiled checked Fortran storage regression,
  pass. Bounds/FPE checks, signalling NaNs and strict 132-column limit retained.
- Full optimized J6 four-rank production interval from checkpoint 4 to `0.3350`
  completed two checkpoint/reload cycles and remaps: exact block-baseline
  checkpoint-field comparison at checkpoints 5 and 6.
- Census executable preserves those production checkpoint fields exactly.
- Checked RK3 and RK4, each oracles off/on: four **first-step** runs passed.
  These do not establish full checked restart/remap coverage for Stage 183.
- At aligned sample 42, `scalar-storage-ready`, four-rank counted scalar
  workspace falls **6.517 -> 2.329 GiB**. Total across the four instrumented
  modules falls **10.487 -> 6.300 GiB**. Full production scalar-record allocation
  is zero across all sampled phases. All ranks have 55 aligned samples.
- These are allocated-capacity counts, not RSS/whole-process memory or a
  simultaneous wall-clock peak. Other modules, automatic temporaries, allocator,
  MPI/OS overhead and unsampled transients are excluded.
- No clean Stage 183 runtime speedup is established. Local memory compression
  persists; four unpinned local ranks do not predict 83-rank cluster scaling.

## Next sequence and eventual completion

### A. Accept and measure this committed checkpoint

Obtain full standard J7 cluster RK3/RK4 checked oracle and checked-production
restart runs, plus separate optimized timing/profile runs. Cluster Stage 183
acceptance is pending; the user authorized this commit as a local checkpoint.
Verify source and binary identities before interpreting results. Keep legacy
unchanged and compare on compatible CPUs, identical inputs and rank placement.
When feasible alternate repeated legacy/block timing pairs in the same
allocation. If only independent `srun` works, report allocation/placement/load
confounding rather than implying a controlled comparison.

Do not start another broad profiling campaign just to rediscover the layout
problem. Existing exclusive timers, work counters, census and samples already
justify the next narrowly scoped audit. Fresh Stage 183 profiles can change its
priority. Legacy/common-boundary work was not the positive block excess in the
pressure-confounded local sample; communication volume alone is not a reliable
optimization objective.

### B. Proposed next substantial computational target (not implemented)

Compile direct scalar column-stencil execution plans so hot loops avoid repeated
coordinate/owner/storage resolution, logical 50-slot decoding, geometry gathers,
and redundant invariant checks. Use shared horizontal geometry, move invariant
work outside scalar/layer loops, and test activity before reading unused fluxes
where existing semantics permit. Preserve ordered numerical operations and
poison/coverage protection with validation at plan construction and appropriate
execution boundaries. Inspect `block_scalar_record_value`,
`resolve_block_scalar_record`, `locate_block_scalar_record` and their consumers.
Do not introduce stale numeric caches or pretend Stage 178 address caching
already removes all these costs. Require explicit reductions in accessor/work
counts, no new per-stage setup collective, exact block-baseline checkpoint gates,
independent oracle checks and measured self-CPU/wall changes.

### C. Finish the computation refactor by measured dependency chains

Inventory each remaining hot Domain producer and compatibility transfer with:
its inputs, owning representation, numerical producer, every consumer, epoch,
mask/topology lifetime, and exclusive CPU/wall cost. Distinguish necessary
single shared primitive work from duplicate calculations and copy-only bridges.
Replace complete high-cost primitive/physics or inverse/forward boundary chains
with native producers and native consumers, carrying direct column data through
the whole chain. Eliminate the obsolete writeback/callback/re-import and storage
only after its final consumer is replaced. Decide geometry-owner versus
final-owner placement from measured transfer and imbalance costs; avoid a blind
ownership rewrite. Audit repeated algorithmic work in either version separately
in isolated candidates, keeping authoritative main unchanged.

Completion criteria: production authoritative numerical state and hot kernels
no longer require dual Domain/block computation; remaining compatibility is
limited to documented low-cost interfaces; topology/mask/RK/remap/restart/polar
coverage is tested; independent numerical acceptance is explicitly reported;
and repeated matched cluster profiles account for residual overhead. A speed
target is not proof that these architectural criteria have been met, nor vice
versa. This is the roadmap, not a claim the remaining work is trivial or done.

## Build/run practice and operational safeguards

- Use fresh named build/bin directories; serial `make -j1 DEBUG=check|false
  PARAM=param_J5 TEST_CASE=climate BUILD_DIR=... BIN_DIR=...`.
- RK3/RK4 selection is `timeint_type` in `src/shared.f90`. Local RK3 tests used
  an isolated source copy; the working source remains RK4. Never run a copied
  parent binary from an overlay whose source now says RK3; check binary identity.
- Optimized flags use `-march=native`: building on another EPYC model caused
  SIGILL before. Build for the execution CPU or use an explicitly compatible
  target. Do not misdiagnose this as a numerical bug.
- Copy the newly built executable only after build success. A successful
  `BUILD_DIR/BIN_DIR` build does not update `bin/climateJ5` automatically.
- Complete current environment settings and source hashes are in the Stage 183
  note. Explicitly set all three validation switches: dynamics, adaptation,
  **remap**. All on for oracle tests, all off for checked/optimized production.
  Clear inherited diagnostic probes. Profiling is a separate run.
- Cluster pattern: `cp ~/wav/bin/climateJ5 .` then
  `srun -n 83 --partition=bb ./climateJ5 simple.in 2>&1 | tee LOG`.
  Use `set -o pipefail`; do not mistake tee's success for solver success.
- Cluster test data may actually be under
  `/mnt/beeond/kevlahan/hydro/climate/`, not `~/hydro/climate/`.
  Resolve the fixture that exists; prior `--fixture ~/hydro/...` failed.
  Outputs must be fresh sibling directories **outside** the immutable seed,
  never children copied recursively back into themselves. Consult
  `experiment.py --help` for independent-srun mode when allocation mode is down.
- Run local MPI simulations sequentially with `mpirun -n 4`, one thread per
  rank/library. Avoid concurrent compilation during timing. Use isolated inputs,
  not the user's standard restart directory. Local MPI sockets may require
  sandbox escalation. Do not escalate permissions by changing the algorithm.
- Shorter `time_end` reduces elapsed work, not necessarily peak live memory.
  J6 retained serious compression; do not call it a pressure-free benchmark.
  Changing tolerance/resume/refinement changes workload and must be labelled.
- Preserve user untracked archives, history, data/topography, generated physics
  build trees, binaries and helper caches. Never `git add .`, clean, hard-reset,
  or overwrite main. Commit only the scoped implementation/tests/evidence.
- An existing optimized `parallel_block.f90` air-temperature descriptor warning
  predates Stage 183. Do not suppress warnings or relax line length to pass builds.

## Evidence and reproducibility locations

Committed JSON summaries preserve identities, numerical gates and measured
results. Large binaries, checkpoints and raw traces below remain **local temp
artifacts**, not Git contents; check existence before reuse. A new chat does
not make them permanent. If absent, reconstruct from committed helpers and
verified original data, and say which evidence is being regenerated.

| Artifact | Path / identity |
| --- | --- |
| Authoritative main | `b82647b36eb825aed0c0aae4e92adcf4ab46d608` |
| Verified legacy archive | `/private/tmp/wavetrisk-stage179b.UypfhC/legacy-main` |
| Corrected block baseline | `/private/tmp/wavetrisk-stage180-block-opt-01`, numerical commit `1b01510ae655eee0e2aab1a8d911e55c147173f2` |
| Evolved J6 immutable seed | `/private/tmp/wavetrisk-stage180-j6-fixture-01/seed` |
| Baseline full comparison outputs | `/private/tmp/wavetrisk-stage180-j6-campaign-01/{legacy,block}` |
| Final Stage 183 optimized archive | `/private/tmp/wavetrisk-stage183-storage-opt-02` (`bin/climate`, `storage-build.log`) |
| Final full production run | `/private/tmp/wavetrisk-stage183-storage-run-02/storage-results.json` |
| Final census build/run | `/private/tmp/wavetrisk-stage183-storage-census-01`, `/private/tmp/wavetrisk-stage183-storage-census-run-01` |
| Checked RK4 binary | `/private/tmp/wavetrisk-stage183-storage-opt-02/bin-checked/climate` |
| Checked RK3 binary | `/private/tmp/wavetrisk-stage183-storage-rk3-check-01/bin-checked/climate` |
| All four checked first steps | `/private/tmp/wavetrisk-stage183-checked-runs-01/checked-results.json` |
| Two-event causal replay | `/private/tmp/wavetrisk-stage182-two-replay-runs-01/legacy` — diagnostic only |

The J6 fixture is genuinely evolved checkpoint **4**, not renamed checkpoint3:
`test_checkpoint_0004.bin.zst`, SHA256
`77218f07fbdec69e53eaa8087b63bff91b59ae0c5541a5550501e3850e6515ae`.
It uses resume4, max_level6, tol0.025, 30 atmospheric and 10 soil levels plus
surface (41 field levels, -10:30). Full validation ends at0.3350; the four
checked first-step cases ended at0.3105. Preserve the original J7 checkpoint3
fixture for cluster validation. No diagnostic physics replay executable is a
production baseline, even if it happens to pass the legacy numerical screen.

Helper suite: `python3 -m unittest discover -s test/parallel_block_profile
-p 'test_*.py'`. Focused inverse/velocity/profile regressions have their own
`run.sh` files. Use existing helpers rather than inventing a new campaign.
Some builders pin older revisions; inspect their CLI/source before assuming
they compile HEAD. Record the actual compiled source and executable hashes.

## Efficient continuation and knowledge boundary

This handoff plus committed code, tests and stage evidence is the durable
engineering state; do not rely on automatic access to the old conversation.
The new chat should first verify `pwd`, branch, HEAD and status, then report
which acceptance evidence it is using. Work in this same saved project/checkout
to use the same repository and local data. If the app chooses a separate
worktree/host, verify paths and copy no untracked data blindly.

Use one explicit optimization hypothesis and bounded validation plan per stage.
Report diagnosis, implementation and acceptance as separate milestones. If work
is expanding beyond the planned test/investigation budget, report what is
complete, why more work is needed and the next bounded step; do not silently
roll several stages into another two-hour response. Maintain concise progress
updates. Do not promise guaranteed full-refactor completion or a runtime ratio;
preserve enough evidence and constraints to pursue them without rediscovery.
