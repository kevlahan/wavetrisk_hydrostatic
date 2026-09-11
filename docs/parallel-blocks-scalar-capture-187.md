# Stage 187: reduce scalar capture storage traffic

Stage 187 retains the Stage 186 oracle ordering correction and removes its
unpromoted zero-import experiment. It targets the scalar restriction capture
timer, which includes patch physics and boundary updates. Numerical formulas,
checkpoint inputs, storage capacities and comparison tolerances are unchanged.

## Changes

Patch physics capture batches the three physics values for all 16 patch nodes
into one validated storage operation. The compact path maps the logical slots
once and assigns the corresponding field array section. Remote producers and
oracle capture retain their previous paths. The storage API also supports
independent full records, with explicit extent, slot and field-boundary checks.

For an established boundary layout, production updates only the three
restricted-flux slots and, for temperature, the three physics slots. Previously
each update read and rewrote a complete 50-slot record. This reduces logical
slot accesses from 100 to three for mass or six for temperature. Layout rebuilds
still install complete geometry; oracle runs retain independent full records.
No persistent cache or additional communication is introduced.

Profiling counts patch batches, replaced node writes, narrow boundary records
and avoided logical slot accesses using an existing reduction. These are
operation counts, not estimates of measured DRAM traffic.

## Component evidence

Both microbenchmarks use the actual storage module compiled with `-O2
-march=native`, approximately 95 MiB of compact field storage and eight
alternating trials. They require matching checksums on every trial. Common
physics calculations and MPI are excluded.

| Storage operation | Original median wall time | Candidate median wall time |
| --- | ---: | ---: |
| Patch physics writes | 0.0143365 s | 0.0125765 s |
| Temperature boundary update | 0.1627465 s | 0.0170005 s |

The isolated boundary storage operation is about 9.6 times faster; patch
batching reduces its isolated median by about 12%. These are component results,
not whole-solver speedups. The observed ranges for original and candidate
boundary operations do not overlap. VM compression was zero; decompression
was 0.125 MiB for the boundary benchmark and 0.234 MiB for patch batching.

- Patch benchmark: `/private/tmp/wavetrisk-stage187-storage-benchmark-01/`
- Boundary benchmark: `/private/tmp/wavetrisk-stage187-boundary-benchmark-01/`

The benchmark drivers are retained as
`test/parallel_block_profile/scalar_patch_batch_benchmark.f90` and
`test/parallel_block_profile/scalar_boundary_storage_benchmark.f90`. To repeat
either benchmark, compile `src/kind.f90`,
`src/parallel_block_scalar_storage.f90` and that driver in a fresh temporary
directory with `gfortran -O2 -march=native`, then run the resulting executable.
Do not run these alongside solver timing trials.

## Validation protocol

All 56 helper tests pass. New tests compare the actual batched storage API
against per-record writes, and the actual boundary routine against a preserved
historical full-record implementation. They cover compact and full layouts,
rebuild and reuse paths, mass and temperature, partial boundary bands,
neighboring records, unchanged geometry and rejected malformed batches, at O0
and O2 with runtime checks.

The immutable Stage 185 four-rank J4/J6 fixture supplies the full 12-step
checkpoint/reload and subsequent-remap gate. Optimized output is compared
exactly with Stage 185 optimized output. Checked production is compared exactly
with the independent full Stage 186 checked oracle result. The short fixture
also compares current checked RK4 oracle-on/off output and current checked RK3
production against the archived independent RK3 oracle result. Short runs have
checkpoint/reload coverage but no subsequent timestep or remap coverage.

The optimized timing protocol excludes a baseline warmup, then runs baseline,
candidate, candidate, baseline sequentially. The baseline contains the oracle
fix without Stage 187 capture changes. Separate baseline and candidate detailed
profiles are excluded from timing comparisons. Every run uses identical input,
seed, grid, rank and thread settings, and must pass exact checkpoint comparison.

## Results

All five validation runs pass. Full optimized RK4 takes 125.2 seconds and
matches Stage 185 optimized output exactly. Full checked RK4 takes 415.2
seconds and matches the independent archived full oracle checkpoint exactly,
including restart and following remap coverage. Current short checked RK4
production/oracle checkpoints match exactly (82.2/212.7 seconds). Current short
checked RK3 production matches its archived independent oracle exactly (88.2
seconds). Validation elapsed times are not used to claim performance gains.

| Measured optimized run | Elapsed seconds | Median ordinary step seconds | VM screen |
| --- | ---: | ---: | --- |
| Baseline 1 | 131.05 | 10.00 | fail |
| Candidate 1 | 117.98 | 8.99 | fail |
| Candidate 2 | 120.96 | 9.06 | fail |
| Baseline 2 | 128.11 | 9.90 | fail |

Both candidate elapsed times are below both baseline times. Median elapsed time
falls from 129.58 to 119.47 seconds, an observed 7.8% reduction. The median of
ordinary-step medians falls from 9.95 to 9.025 seconds (9.3%). All exact numerical
gates pass. The excluded warmup takes 131.32 seconds; separate baseline and
candidate profiles take 126.10 and 120.71 seconds.

All four measured runs exceed the existing 1 GiB system-wide decompression
threshold (baseline 4.24/2.05 GiB, candidate 2.30/3.70 GiB). These observations
therefore do not establish a clean whole-solver speedup, nor predict 83-task
cluster performance. No threshold was changed to accept the timings.

The separate profiles provide supporting component evidence: summed
maximum-rank scalar capture time across the two reporting windows falls from
13.0953 to 5.5973 seconds (57.3%). Both windows improve: 8.5379 to 3.6716 seconds
and 4.5574 to 1.9257 seconds. This timer includes both patch and boundary capture;
it does not isolate either optimization and is not an additive wall-time budget.

Counters across four ranks verify 6,679,920 patch batches replacing 106,878,720
per-node calls, and 118,240,560 narrow boundary updates avoiding 11,291,973,480
logical slot accesses. The latter counts software operations, not physical
memory bytes. The original full-record oracle and rebuild paths remain intact.

The scalar capture change is retained based on repeated isolated gains, both
profile windows, exact numerical checks and consistent ABBA separation. The
oracle correction remains; the Stage 186 zero-import experiment remains
unpromoted. The machine-readable evidence is in
`docs/parallel-blocks-scalar-capture-187-results.json`.

Raw validation: `/private/tmp/wavetrisk-stage187-validation-01/`.
Raw timing: `/private/tmp/wavetrisk-stage187-timing-01/`.
Drivers: `/private/tmp/wavetrisk-stage187-validation.py` and
`/private/tmp/wavetrisk-stage187-timing.py`.

The original legacy-screen failure and earlier checked-versus-optimized
mismatch remain unresolved and separate from the corrected oracle-on/off
mismatch. Checked solver flags do not extend to the external physics package.
The subsequent standard 83-task campaign passed in job 54452, with 4.71% lower
median elapsed time and 7.88% lower median ordinary-step time across three runs
per version. See [the cluster result](parallel-blocks-cluster-187-54452.md).
