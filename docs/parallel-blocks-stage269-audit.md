# Stage269: bounded current-cost audit before the next removal

Prepared and locally validated; bbserv measurements pending. The production
solver remains accepted Stage266 `ff740744`. This first step implements the
[Stage268 prioritization protocol](parallel-blocks-stage268-priorities.md),
measuring all major remaining costs together before selecting one intervention.
The current qualified main gap is 206.695 ms per ordinary step, from job 55398;
this audit establishes no new speedup or Domain independence.

## Implementation

The [campaign tools](../test/parallel_block_profile/cost_audit_269/README.md)
generate separate bounded observer overlays for Stage266 and read-only main
`6a245db7`. Exact uninstrumented controls are bundled and remain source-identical
to the Stage266 campaign, including main's established checkpoint I/O hardening.
Observers modify only timing calls, work accounting, final output and module
build dependencies in their diagnostic copies. No arithmetic, transport,
cache lifetime, boundary requirement or ownership selector changes.

A fixed-size Fortran stack records self/inclusive monotonic wall time, process
CPU and call counts for 27 transaction families. Timers occur at routine
boundaries; there are no node/layer clocks, new synchronization or step-time
output. Per-rank records are flushed after all steps. Scope imbalance, step/
stack overflow and invalid settings fail explicitly. Work counters are added
once per actual cache rebuild, plan rebuild, publication, installation and
inverse exchange. They report rebuilt columns/layers/blocks, logical slots/
values, actual installations, remote requests and logical double payloads.

Measured families cover complete wavelet transactions and Domain/native inverse,
setup/publication/installation/wait, provisional and accepted state refresh,
adaptation, catalogue/source extraction/migration/post-grid reconstruction,
hydrostatic cache and native thermodynamic/velocity recomposition, compression,
forward transforms/ghost refresh, surface pressure and Domain boundaries.
Nested self times conserve the observed local step duration. Inclusive regions
overlap. The post-grid scope's self time excludes separately observed catalogue,
source, migration and other child routines.

## Campaign and decision rule

One exclusive nk/bb13 allocation builds five solvers serially, sharing all 23
physics objects. There are seven correctness runs and eight 83-rank calibration
runs, fifteen total. Correctness includes both 16-rank controls, optimized main/
block observers, checked block observer, and independent block oracles at 16 and
167 ranks. The latter covers zero-Domain ranks. Existing checkpoint, VTK,
trajectory, restart/remap, storage and physical placement gates remain. The
oracle proves all 48 inverse returns and 84 retained installations.

Calibration order is `main, main-observer, main-observer, main, Stage266,
block-observer, block-observer, Stage266`, with two independent runs per binary.
The ordinary-step window excludes steps 1, 8 and 9. Each observer has its own
prespecified overhead gate: mean shift within max(3 ms, 1%), individual paired
shifts within max(5 ms, 2%), and control/observer spread within 2%, relative to
the corresponding control mean. These preserve the earlier Stage229 resolution
policy; no thresholds are chosen after observing this campaign's results.

Both calibrations must pass before region differences inform priority. Failed
calibration remains `observer_costs_unqualified`, while science checks and raw
records remain available. The independent review recomputes raw completeness,
self/inclusive consistency, conservation, transaction counts, ordinary timings
and the signed current Stage266-minus-main attribution. There is no automatic
repeat, optimization acceptance or promotion.

After review, select the largest demonstrably redundant operation with a bounded
replacement and clear consumer/epoch proof. Lazy provisional hydrostatic cache
construction is one hypothesis; wavelet or reconstruction work takes precedence
if it offers a larger recoverable cost. Preserve the Stage267 required shared
readers and full publication until their complete replacements qualify. Test
one intervention against uninstrumented Stage266 and main, using the established
short thirteen-run campaign, before accepting it or stacking another removal.

## Validation and limits

All three new local builds passed with `make -j1`. Four real J4 replays ran
serially with four MPI ranks and one library thread: optimized main observer,
optimized block observer, checked block observer and independent block oracle.
Each matches its source-bound reference exactly, including step states,
checkpoint and all 62 VTK members. Every rank has all 12 × 27 records, conserved
self times and the expected provisional/accepted/installation counts. The
oracle retains all 48 exact-return proofs and 84 installations.

Nine workflow tests cover five shared-physics builds, the complete fifteen-run
schedule and observer environment, duplicate guards, purged status, early
failure/raw-file collection, independent archive review, calibration failure,
signed offsets, conservation, one-rank composition, source scope and early
return wrapping. A real Fortran test at O0 and O2 covers nested conservation,
disabled output, bad selectors, mismatched scopes and buffer overflow. Tests
were repeated from the extracted delivery archive; helper imports and shell
syntax passed. Python tar creation avoids macOS extended-attribute headers.

These local checks establish correctness and observer completeness. They do
not predict bb13 overhead, costs or recoverable gain. Signed mean self-region
differences add to a mean local-duration difference, not to the difference of
run-level ordinary medians. The longest-local-step composition selects one
rank per step; it does not sum independent region maxima or prove the global
critical path. CPU includes MPI spinning. Counters cover named products and
inverse payload, not all solver traffic or physical memory bytes. Profiling
locates hypotheses; an uninstrumented single-change comparison establishes
net savings.

The kit requires only the established J5 fixture on bbserv, not old result or
kit directories. Its launcher records the submitted ID and a unique results
path, rejects duplicates until collection, and collects diagnostics even after
an early failure. Review archives omit source/binary/checkpoint/VTK payloads;
cluster-computed hashes and numerical reports are verified, not rerun locally.
See the [sealed preparation record](parallel-blocks-stage269-preparation.json)
for archive identity and local evidence.
