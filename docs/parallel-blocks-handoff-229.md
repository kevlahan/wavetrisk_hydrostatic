# Stage229 handoff

User authorized preparing the recommended Stage203 requalification campaign,
with the goal of identifying recoverable extra work by execution location.

Current integration: the user authorized promoting the Stage229 result after
review. Production build defaults are now `STAGE202=on STAGE203=on`, on the
Stage226 numerical source (`1d7bc050`). Main must remain `b82647b3`. The inverse
workspace and scalar batching implementations were already committed; only
build defaults and documentation change. `STAGE202=on STAGE203=off` retains
the comparison baseline. Stage227 is not included. Stage228's publication-copy
performance target is closed. The campaign preparation/return history below
describes the earlier frozen baseline, not the new default.

Use `parallel-blocks-stage229.md` for the experiment and
`parallel-blocks-cluster-229.md` for complete commands. Ten fresh binaries;
normal successful path 29 runs, maximum 33 only for disputed/near-threshold
confirmation. J4 41/40 capability first, J5 correctness and 16-task screen,
early 83 pair, fresh ABBA. Target 25 ms with 5 ms near-threshold band. Old
Stage203 gains (40–55 ms) justify retesting, not present-day savings claims.

The fixed coarse location lane is A B L AO BO LO LO BO AO L B A at83. If the
candidate stops at16, omit B/BO and still obtain a current A/L comparison.
Numerical/infrastructure failure stops immediately. Calibration failure retains
raw data and labels attribution unqualified. No automatic code promotion,
variant, extra observer pass or deeper instrumentation.

Review confirmed bare A/B first, then signed coarse A−L, B−L and A−B differences.
Eight disjoint locations are built from six outer phases and WT/tendency batches.
Keep both positive and negative rows and the mean/median/observer bridges. All
rank events/phase sums and four WT/tendency calls per ordinary step must pass.
The tendency contexts have implementation-specific boundaries. A positive
remaining location is not a disposable budget; MPI waiting is not independent
savings. No feature ablation identifies inverse versus field-batching gains.

Retain high-precision logs and actual saving/reload/remap events. Exclude steps
1,8,9; retain 2–7 and10–12. Do not treat timesteps as repetitions or this trajectory
as a continuous no-output run. Same-count block and own-observer comparisons are
exact; block/legacy uses unchanged numerical limits and exact topology/soil.

All builds/runs use the same nk host and physical-core map, full-node256 allocation,
no step --exact/-c. Exclude bb04. Legacy is exported read-only with the old common
checkpoint I/O hardening; main is never checked out or edited. The collection
wrapper creates/verifies/reuses the archive and transfers bbserv → murier.

Preparation status and final hashes belong in `parallel-blocks-stage229-kit.json`;
local solver evidence belongs in `parallel-blocks-stage229-local-results.json`.
Do not claim a prepared, validated kit until those records pass. Local MPI timings
are not performance evidence; initial missing-grid and contended oracle attempts
are preserved separately and excluded from validation claims.

## Prepared and verified

All eight local cases and the fixed legacy screen passed. In particular, the
optimized Stage203 full-oracle run at41 tasks/40 Domains completed all12 steps,
checkpoint loads2/3 and both remaps, with an exact final checkpoint. The checked
four-rank run passed, and all three coarse observers match their own references
exactly and pass all-rank event/phase conservation.

The sealed package passes13 helper tests, the real extraction/launcher/collector
workflow using fake Slurm/scp, and the ten-build orchestration using compiler
stubs. Local real solver builds/runs are the numerical evidence; stubbed builds
are only packaging/selector tests. No cluster job has been submitted.

Desktop kit: `/Users/kevlahan/Desktop/wavetrisk-stage229-cluster-kit.tar.gz`

SHA256: `c031ebecbadce7fa67af08d42f18fa4779f2a65736613da4df61c6b1faf43487`

Bytes: 122710857. Source reproduction and Desktop-copy hash verification passed.
Main and production source remain unchanged.

## Return reviewed — job 54737

See `parallel-blocks-stage229-review.md` and `.json`; the return reviewer is
`test/parallel_block_profile/review_cluster_229.py`. All 4,373 file hashes, ten
build identities, 29 runs, 54 numerical reports and all coarse profiles were
rechecked. Fifty comparisons are exact; four block/legacy comparisons pass the
unchanged tolerance policy. J4 41/40 checked/full-oracle capability passes.

Performance **meets_target**: fresh83 ABBA gives 0.806890 → 0.761174 s, a45.715 ms
(5.666%) gain. Paired gains46.652/44.779 ms exceed25 ms; no extra batch needed.
The16-task single pair gains243.254 ms (9.464%). Block-only observer differences
locate approximately34.095 ms in wavelets and13.353 ms in tendency preparation;
these are not separate feature-ablation results.

Fresh bare diagnostic lane: baseline0.813473, candidate0.758268, legacy0.542883 s;
candidate–legacy gap215.386 ms. **Legacy observer fails overhead calibration**
(+11.223 ms mean versus5.429 ms limit); overall attribution is unqualified.
Provisional candidate–legacy signed locations: wavelets+136.817,
adapt/reconstruct+77.664, otherdynamics+19.765, tendency−27.882, other−4.161 ms.
Keep negative offsets and mean/median/observer bridges. These are source-search
priorities, not recoverable budgets, and do not alter bare candidate acceptance.

Recommend production integration of combinedStage203, then a bounded source
audit of remaining complete wavelet representation transitions and their
consumers, checking reuse across adaptation/reconstruction. No automatic
promotion, commit, solver edit or new campaign was made on this attachment turn.
Production remains1d7bc050, Stage203off; main remainsb82647b3. Do not reopen
closed tiny-copy/cache/RK-closure targets or launch another broad observer pass.
Any next prototype needs a consumer/epoch/replacement-cost argument.

## Authorized production integration

Following the user's explicit commit request, promote both build defaults to on,
retain all override choices and the Stage203/Stage202 compatibility check, and
document the existing sparse-outer incompatibility. Numerical implementations,
oracles and all Stage226 above-Domain capability guards are unchanged. Existing
Stage229 real local/cluster checks qualify the exact enabled configuration;
build-selection verification confirms that the new defaults select those flags
in optimized and checked modes. See `parallel-blocks-stage229-integration.json`.
Use fresh BUILD_DIR/BIN_DIR when switching options; old objects are not rebuilt
automatically when flags change. No new optimization or campaign is included.
