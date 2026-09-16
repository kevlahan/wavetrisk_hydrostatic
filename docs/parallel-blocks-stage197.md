# Stage 197: keep scalar wavelet coefficients resident

Current direction: the user has required diagnosis of the main legacy/block cost causes before any further optimization. The [diagnosis objective and evidence requirements](parallel-blocks-performance-gap-diagnosis.md) supersede optimization recommendations below. Existing results are retained as historical evidence.

Stage 196 baseline: `9bc3a19556dab79b3cba3944cf5df53517389da3`.
Stage 197 is an experimental production change plus focused measurements and
a shortened cluster campaign. The default remains the validated previous path;
set `WAVETRISK_STAGE197_RESIDENT_WAVELETS=1` for the candidate. Selection is
checked collectively and printed once. Cluster scripts set it explicitly.

## Production change

At every native WT call, the old path imports the Domain wavelet field after
WT_after_step has set it to zero. The candidate initializes the native wavelet
arrays locally instead. This eliminates communication, packing and remote
installation of structural zeros while preserving the RK recovery checkpoint.

Scalar transform outputs then stay in compact block storage through compression.
The complete vector transform still produces its Domain-owner stage. A narrowed
reverse import supplies only this vector family: no scalar buffer clearing,
packing, scalar collective or scalar installation is performed. Scalar buffer
positions are still advanced and checked against the existing route extents.

Production no longer publishes the uncompressed scalar family before immediately
reading it back. Compressed complete coefficients are published at the existing
consumer boundary. The native inverse continues to consume the resident compact
coefficients and its existing alias/ghost routes. Oracle-enabled runs additionally
materialize the uncompressed Domain output so the established independent
uncompressed-transform checks remain available.

Per-stage transaction accounting distinguishes this validation publication from
production. A readiness flag and writeback-plan generation guard protect the
resident coefficient lifetime. The flag is consumed after compressed publication
and cleared when plans are retired, including restart/grid reconstruction. The
retained arrays are existing block storage; no second coefficient cache or
long-lived payload pointer is introduced. Scalar/vector arithmetic, masking
thresholds, RK order, checkpoint format and comparison tolerances are unchanged.

This implements the first ownership-lifetime change in the proposal. It does not
eliminate the complete inverse transfer chain. Vector coefficients still require
a transfer from their producer owners, and inverse aliases remain required.
Stage 196's 0.205488 s wavelet transfer excess is the containing target, not a
prediction of the savings from this particular change. Cluster timings determine
whether removing these transitions produces a useful whole-workflow gain.

Production edits are confined to `src/parallel_block_mpi.f90`: initialization,
publication, vector-only import/packing, and transaction counts. The selector-off
path retains the previous payload operations. The independent cluster reference
is built from the untouched Stage 196 commit.

## Focused measurement changes

Observers are generated only into disposable build trees. Bare candidate and
reference binaries do not contain the new measurement module. The same observer
binary with the candidate selector off supplies the diagnostic Stage 196 control;
its checkpoint must match the untouched Stage 196 executable. The legacy
observer is generated from the pinned legacy with the common checkpoint-I/O fix.

The Stage 196 clock/conservation mechanism remains. Additional scopes separate
local coefficient initialization, direct scalar flux, restriction, divergence
and temperature-boundary closure. Generic boundary scan/zero and local geometry
counters previously left behind old-profiler guards are now independently
enabled. Parser coverage checks reject missing generic boundary counters when
send/receive work proves those paths were exercised.

New counters include native coefficient doubles initialized, scalar imports
avoided, uncompressed publications avoided, direct scalar samples, restricted
flux edges, divergence samples/evaluations and temperature flux node evaluations.
In each profiled candidate step, avoided-import/publication counts must equal
the number of dynamics WT calls; the selector-off profile must report zero.
This prevents an accidentally disabled candidate from passing as an optimization.

Legacy temperature flux counts include the two regular step1 computation branches.
Native counts cover the actual valid-temperature node evaluations after the
mass shortcut. These distinguish sampled storage from performed temperature
flux work. Divergence counts are collected directly from production stencil columns
and active masks; the old oracle-only comparison statistics are not used.
They are source-level work counts, not instruction counts or an
exhaustive accounting of special boundary kernels. Legacy direct-flux wall scope
still encloses shared step1 producer work, so its ratio to the native direct-flux
scope is not an arithmetic-only slowdown. Shared legacy fusion limits that
comparison; retain the full WT/outside-WT envelopes and the explicit scope names.
The observations provide a narrower next computational target without rewriting
that kernel in this ownership experiment.

`dependencies197-rank-N.txt` records inverse gather/scatter/alias phase,
scalar/vector component, field family, level range, request count, sent/received
values, local post-to-wait interval and completion interval. Empty completions
are retained separately. The parser checks rank/step coverage, finite values,
sequence and agreement with the independently counted inverse completion scopes.
Local intervals do not establish globally synchronized arrival or useful CPU
work: they may contain instrumented preparation or progress. No barrier or
measurement MPI call is inserted, and no clock is read inside the cell loops.
A bounded buffer fails explicitly on overflow and is written after simulation.

## Short campaign

The new cluster runner defaults to 14 executions: two checked runs, three
reference/warmup runs, two alternating rounds across Stage 196/candidate/legacy,
and one focused profile for each implementation. All exact checkpoint, workload,
post-restart remap, checked/optimized and fixed legacy gates remain. Warmups,
checked runs and profiles are excluded from timing summaries.

Six fresh binaries are built once: untouched Stage 196, bare candidate, checked
candidate, candidate observer, legacy, and legacy observer. The candidate observer
is used twice with explicit selector settings for the baseline/candidate profiles.
Build manifests verify sources, patches, binary identities, flags and common
external physics objects. Provenance is passed by filename, avoiding the previous
argument-size problem. Failures leave partial results and logs for collection.

The two timing rounds are a development screen, not a confidence interval. Compare
whole-step wall times and targeted work reductions. The one-profile-per-variant
overhead comparison is also a screen; extend just the affected measurements if
profiling perturbs the result or the apparent gain overlaps run variability.
The legacy profile contains extra counter updates in leaf routines but no leaf
clocks; inspect its overhead before interpreting small timing differences.

See [cluster commands](parallel-blocks-cluster-197.md). Future results retain the
single inventory-bearing tarball collection and one scp transfer to murier.
Local verification and frozen kit identities are recorded below.

## Completed local validation

Seven J4/J6, four-rank RK4 cases passed, with contraction off: optimized selector
on/off, checked candidate, optimized candidate with validation oracles, block
profiles on/off, and legacy profile. Each completed 12 steps, saved/reloaded
checkpoint 3 after starting from checkpoint 2, and exercised two remaps including
one after restart. All six block cases match the independent Stage 196 checkpoint
exactly; the legacy observer matches the independent legacy checkpoint exactly.
Printed workloads match, and the fixed cross-legacy numerical screen passes.
This is not a claim of bitwise equality between blocks and legacy.

The repaired generic boundary counters report nonzero work. Each block profile
contains 8,736 inverse dependency records, matching its completion scopes.
Candidate counters confirm four avoided scalar imports and four avoided
uncompressed publications per rank per ordinary step; the control reports zero.
Production divergence sampling/evaluation counts are now nonzero and measured
at the kernel, independently of oracle enablement. Timing envelopes and restart/
remap flags pass parser checks. These local runs establish correctness and
measurement coverage, not a performance gain.

Ten automated tests cover build provenance, the 14-run/34-gate mocked campaign,
launch/correctness failures, frozen legacy tolerances, oversized provenance,
compiled observer modes at O0/O2, corrupted records and missing production
counts. Six cluster builds are tested with mocks locally; the seven real MPI
cases above use the locally compiled production/checked/observer executables.
The 23 external physics object hashes match across the local source trees.

The frozen patches reconstruct the recorded candidate, block observer and
legacy observer sources exactly. Those inventories also match the local build
sources. Bundle verification, standalone helper imports, shell syntax and
archive/helper hash checks pass.

- Acceptance evidence: [local results](parallel-blocks-stage197-local-results.json).
- Frozen identities: [kit manifest](parallel-blocks-stage197-kit.json).
- Archive: `/private/tmp/wavetrisk-stage197-cluster-kit.tar.gz`.
- Archive SHA-256: `c62cdea9ea76d15dfc1c02d4ef2fff93dbd8e8d56743cf2b96b2e802de235248`.

## Submission guard repair after job 54662

The observer patch changes the Makefile to compile its measurement module. The
submission helper verified the complete observer source inventory, then wrongly
compared the Makefile to the earlier production inventory. Compare it to
`profiled_sources["Makefile"]` instead. No numerical source, build flags or binary
changes are required. The regression fixture now includes the observer Makefile
change: it reproduces the original failure, accepts the corrected handoff, and
still rejects a subsequent unexpected Makefile edit or binary replacement.
All three build/submission tests pass. The kit archive and helper hash have been
updated; the original archive is preserved separately in /private/tmp.
