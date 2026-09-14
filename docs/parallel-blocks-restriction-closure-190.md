# Stage 190: column restriction closure and direct legacy comparison

Baseline: Stage 189 `a0c75adb4dc1741f56ecb572ab1809e8a70a32e4`.
All solver comparisons use `FP_CONTRACT=off`. Cluster campaign 54553 passed; see [the cluster review](parallel-blocks-cluster-190-54553.md).

## Production change

Stage 189 batched only interpolation. Stage 190 also computes partial/coarse
closure across physical columns and writes the resulting parent flux column
through `scalar_write_column`. The original scalar/oracle routines are unchanged.
The only production source change is `src/parallel_block_mpi.f90`.

For each active parent edge, the new closure constructs horizontal indices and
computes overlap geometry, normalization and geometric coefficients once.
Compact production storage guarantees that this geometry is shared across
physical layers. The arithmetic expressions and their ordered products, sums
and differences are retained. Full-record oracle geometry stays on the original
scalar path.

An invocation-local cache holds up to 32 field columns, keyed by horizontal
coordinates and field slot. Repeated flux/tendency consumers reuse the same
column instead of resolving patch/boundary/ghost storage again. For 30 layers,
the cache's field array is 7.5 KiB, plus bounded indices and arithmetic scratch.
There is no persistent cache, generation invalidation, new MPI collective or
additional retained solver state. The existing profile reduction carries four
new counters: closure edges, column locations, reused columns, fallback layers.

The transported-tendency sentinel remains significant: only layers with an
unavailable tendency reconstruct it from the fine-level flux stencil, using the
original scalar divergence routine. Other layers do not acquire extra stencil
reads or boundary dependencies. Reconstructed columns are reused only within
the current edge invocation. The finalized fine-level inputs are read before
parent writes, preserving restriction-level and publication ordering. Finite,
poison, coordinate and extent checks remain. Interpolation outputs are checked
before closure, including values not used by an individual parent edge.

Parent values are written directly by column; the old four-element strided-row
temporary at the scalar closure call is removed. No integrator, precision,
tolerance, ownership, checkpoint format or production I/O policy changes.

## Validation

`test_restriction_closure.py` extracts the actual new closure and unchanged
scalar partial/coarse functions. It uses the actual compact storage API and
restriction index tables with a synthetic location resolver covering patch,
boundary and ghost storage. At O0 and O2 with contraction off, it compares bits
for all three edge orientations, multiple coordinates, scalar variables and
column extents. Mixed available/missing tendencies exercise layer-specific
fallback. Negative cases reject nonfinite/poisoned input, invalid extents and
full-record storage on the compact-only column path.

Local solver acceptance uses the evolved J4/J6, four-rank, 12-step RK4 window,
checkpoint 2 -> 3/reload and subsequent remap. The optimized result must match
Stage 189 exactly; checked output must match optimized output and the archived
independent RK4 oracle exactly. The scalar oracle is unchanged, so comparisons
against its retained output do not claim a freshly executed full oracle run.
A separate checked RK3 smoke window covers checkpoint/reload; its short duration
does not cover a post-restart remap. Separate optimized profiling checks work
counters and exact checkpoint output.

The separate legacy comparison archive changes only checkpoint I/O, using the
already validated Stage 187 fix. The authoritative legacy archive remains
unchanged. Its local checkpoint must match the archived contraction-off legacy
output exactly, and the Stage 190/legacy comparison must pass the existing fixed
screen. No tolerance is calibrated to the candidate.

Results and precise coverage are recorded in
`parallel-blocks-restriction-closure-190-results.json`. Full local solver elapsed
times are validation observations, not controlled block/legacy timing ratios:
VM compression and compilation overlap can affect this machine.

The pre-commit helper suite passes all 86 tests (14 September 2026), including the runner repair. The full optimized/checked RK4 comparisons, independent
oracle comparisons, RK3 smoke comparisons, legacy I/O transparency comparison
and fixed legacy screen pass. The final optimized profile also completes all
12 steps and matches Stage 189 and checked checkpoints exactly. Across its two
windows, 312,800 closure edges use 5,004,800 field-column locations and 1,251,200
cached-column reuses: 16 locations and four reuses per edge. The original scalar
fallback reconstructs 73,680 missing layer tendencies.

The full off-profile runs preceded a one-line diagnostic-counter correction
caught by the interpolation component test: the interpolation read count belongs
to counter 8, not the new closure fallback counter 12. Numerical expressions were
unchanged and that assignment is disabled in off-profile runs. All three final
block binaries were rebuilt; the final full-window profile and the passing full
helper suite cover the corrected source. The result record retains both source
identities rather than relabelling earlier runs as runs of the final binaries.

## Component benchmark

`benchmark_restriction_closure.py` times actual scalar/new closure routines with
synthetic routing and compact storage. Eight alternating trials each evaluate
6,000 edges across 30 layers, including missing-tendency fallback. Checksums
must match. The measured local medians are 0.123847 s scalar and 0.0399885 s
column (67.71% reduction). The memory screen passes: 22.77 MiB system-wide
decompression, no compression/swap activity, and 3.70 MiB sampled benchmark RSS.
This is a bounded component result, not a solver or cluster speedup prediction.

## Cluster campaign

See [the cluster protocol](parallel-blocks-cluster-190.md). It builds on `nk` to
make native CPU targeting explicit, freezes binaries/sources/physics objects,
and compares Stage 189, Stage 190 and legacy Domain with normalized checkpoint
I/O in one 83-rank allocation. It rotates their timing order, requires matching
printed step/dt/Jmax/dof trajectories, retains exact block/oracle/cross-mode
gates and applies the existing legacy numerical screen. It records resource and
placement evidence and reports both block/legacy ratios. Campaign 54553 passes the exact block comparisons and fixed legacy numerical
screen, with a 2.14% elapsed and 4.03% ordinary-step reduction versus Stage 189.
