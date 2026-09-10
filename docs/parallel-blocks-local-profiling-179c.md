# Local comparison and experimental legacy instrumentation — 10 September 2026

This is profiling work, not a new numerical/performance implementation. No
commit is made. The authoritative `main` revision remains
`b82647b36eb825aed0c0aae4e92adcf4ab46d608`. The user's clarified permission allows
instrumentation in isolated legacy copies, not edits to the authoritative tree.
The existing Stage179b cluster experiment remains available unchanged.

Durable copies of the result JSON, memory observations, both stack samples and
the legacy instrumentation patch are in `docs/profiling-stage179c/`. Temporary
paths below retain full run directories, executables and simulation outputs.
Executable SHA-256 identities used in these experiments:

- Unmodified legacy: `0f99636210320cdd79cff5b49cb3becd0618faa9f042b963de0589e3990a277e`
- Instrumented legacy: `aaeb264e97a149444493ab5d67d952b6b8e926c49be4ebcc62404e21619b0858`
- Block: `8de9b3ac09433902eca952675d74a03fd7f072a18457a33590845d55e2bb72cc`

## Reproduction and identities

The experimental legacy source/build is
`/private/tmp/wavetrisk-local-179c-legacy-04`. It was copied from the previously
compiled, verified main archive at
`/private/tmp/wavetrisk-stage179b.UypfhC/legacy-main`, not checked out over the
parallel-blocks workspace. `experimental-instrumentation.patch` records its
tracked-source changes. The original archive is verified against all tracked
main blobs before and after preparation. No main references are moved.

New helpers under `test/parallel_block_profile`:

- `prepare_legacy_profile.py`: creates a new external experimental copy of a
  verified built archive; adds anchored coarse scopes and transplants only the
  boundary instrumentation diff. Build that copy with its own Makefile, `-j1`.
- `legacy_profile_io.f90`: uses the same exclusive wall/process-CPU timer core;
  enables it only at simulation start and writes rank-local records after the
  simulation. No added timestep MPI, barriers or per-node clocks.
- `local_experiment.py`: fresh input copies, sequential four-rank runs,
  alternating unchanged-legacy/block pairs and separate instrumentation checks.
- `local_memory.py`: system-wide macOS VM counters and sampled process RSS.
- `local_sample.py`: separate rank-0 macOS call-stack sampling executions,
  excluded from timing pairs.
- `test_legacy_profile.py`: exact-anchor and wrapper transformation tests.

Do not use an instrumented legacy executable as the authoritative comparator.
The experimental helper does not issue a `tracked_contents_unchanged` build
identity for its modified source. Numerical formulas and call ordering are
unchanged by this instrumentation. IDs reuse the block profiler vocabulary,
but a similarly named legacy region is not necessarily the same set of work.
In particular, legacy scalar replay includes Bernoulli/Exner restriction, while
the block scope includes direct-flux/divergence work. Do not subtract individual
same-ID rows as if they were perfectly equivalent kernels.

## Completed repeated short experiment

Raw results: `/private/tmp/wavetrisk-local-179c-runs-01/results.json`.
Memory observations: sibling `memory-observations.jsonl`.

Host: iMac21,1, 16 GiB RAM, four performance and four efficiency cores. Four MPI
ranks, one thread per rank, GNU Fortran optimized builds. Compilation and other
simulations were finished before timing. macOS scheduling/core placement is not
pinned or verified. Ordinary desktop background activity remains.

The frozen checkpoint-3 fixture runs to `time_end=0.3040`, giving three ordinary
RK4 steps, with identical input files and printed numerical states in all nine
executions. Every run gets its own directory and executable copy. Startup is
excluded from the application step times below.

| Pair | Order | Legacy mean step, s | Block mean step, s | Ratio |
|---|---|---:|---:|---:|
| 1 | legacy, block | 4.577 | 32.133 | 7.02 |
| 2 | block, legacy | 5.010 | 31.067 | 6.20 |
| 3 | legacy, block | 4.763 | 23.133 | 4.86 |

These are local observations, not a stable benchmark or an estimate of the
83-rank gap. In particular, the large variation is not explained by compiler
overlap in this experiment. The memory monitor recorded about 283 GB of logical
page decompression across one complete block-run observation interval, versus
about 1.4 GB and 0.06 GB in two legacy intervals. These are SYSTEM-WIDE counter
deltas over different-duration intervals, not measured application DRAM traffic
or a causal attribution of every decompressed page to WAVETRISK. They establish
serious memory-pressure confounding. Sampled RSS excludes compressed footprint.

One experimental-legacy detail-off run averaged 5.230 s/step; detail-on averaged
5.177 s/step. This single pair does not establish a tight instrumentation-overhead
bound, especially relative to the completely uninstrumented executable. The
block detail run averaged 23.767 s/step. Its timings are not an independent
speedup result.

## Matched phase evidence

Legacy self-wall averages below are across three steps and four ranks; block
averages use the same workload. Legacy timestep wall is 5.177 s; block is 23.746 s.
Timer sums conserve each local timestep (block printed-vector rounding below
1.4e-5 s; legacy full-precision errors below 1e-6 s).

The common boundary routines have identical subphase instrumentation and are
the cleanest direct comparison:

| Boundary self work, s/step/rank average | Legacy | Block |
|---|---:|---:|
| All measured boundary subphases including waits | 2.974 | 2.344 |
| Excluding the separately timed boundary wait | 2.126 | 1.549 |
| Boundary MPI wait alone | 0.848 | 0.795 |

The block code's remaining Domain-boundary path is not larger in this sample.
Its transmitted doubles were 106,595,355 versus legacy's 209,004,015 across all
ranks/three steps. Reducing this common path may benefit both versions, but it
does not explain the positive block time excess in this experiment.

Leading block-exclusive self costs are scalar/inner inverse kernels (2.910 s),
scalar replay (2.827 s), block ghost packing/installing (2.167 s), shared producer
residual (1.709 s), and native wavelet work (1.226 s). The first three total
7.904 s, about 33% of the block timestep. They include necessary numerical work
and memory-management consequences; they are not entirely removable overhead.

## Structural candidates, not implemented optimizations

1. **Stop expanding shared geometry into every scalar/layer record.** The
   production transport already separates 33 shared and 17 field-dependent
   values, but the computational workspace expands them back to 50 values per
   field/node. This short run's selected expanded arrays peaked at 7.138 GB
   summed across rank peaks; geometry accounts for 4.711 GB, versus a 57.45 MB
   shared-once estimate. Keeping all 17 field slots and sharing only those 33
   geometry slots estimates 2.484 GB instead of 7.138 GB for these arrays: about
   65% smaller, NOT a 65% runtime saving. Peak sums need not be simultaneous and
   do not represent total process memory. Audit all inverse/restriction/boundary
   consumers; do not simply drop mass/soil slots still used by inverse work.

2. **Compile direct column stencils and move invariants outside field loops.**
   The scalar accessor repeats coordinate/storage resolution, extent checks,
   finite/poison checks and record-stride arithmetic per scalar read. The
   divergence kernel reads six fluxes and area before testing its active flag,
   unlike the legacy early mask test. These are specific opportunities for
   cheaper native execution; retain validation at plan construction and oracle
   boundaries rather than weakening checks without replacement.

3. **Investigate common algorithm repetition independently.** Legacy
   `cpt_or_restr_u_source` can invoke the same whole-patch source calculation
   once per absent child before its restriction pass. The climate velocity
   physics callback also repeats geometry-dependent pentagon-distance work
   across layers. Count actual eligible work, audit callback side effects and
   ordered writes, then validate an isolated candidate against unchanged legacy.
   Do not assume these are large costs or blindly regroup the native ordered
   velocity action tape.

4. **Common boundary route scans and buffer clearing remain secondary targets.**
   They are real costs in legacy too. The 100x-scale logical buffer-clear/receive
   ratios are not hardware bandwidth measurements: repeated writes can hit
   cache. No claim that clearing alone explains the block deficit is justified.

The memory/layout and accessor chain is the strongest block-specific next
investigation. The matched 83-rank run and hardware samples are still required
to establish its priority and attainable speedup on the cluster.

## Separate call-stack samples

Directory: `/private/tmp/wavetrisk-local-179c-samples-01`.
Both executable samples start after the simulation marker and target rank 0.
They are separate from the nine timing/profile executions. Sampling interval is
10 ms, requested duration 20 s; the shorter legacy simulation exits earlier.
These snapshots include sleeping threads, inlining and partially sampled RK
stages. They are triage evidence, not comparable useful-CPU percentages.

In the block sample, `block_scalar_record_value` is the hottest application
leaf (230 observations out of 1,831 main-thread observations). Separate leaves
include `resolve_block_scalar_record` (46), `locate_block_scalar_record` (14),
`install_buffer_ghost` (112) and `fill_boundary_node` (81). This corroborates
the profile's accessor/installation hotspots without adding clocks to those
inner routines. The sample's 1,831 `kevent` observations belong to a background
thread, not 100% idle numerical execution; do not mix thread denominators.

Legacy samples instead prominently show local boundary copying, boundary
packing and MPI progress, consistent with its exclusive timers. Allocation
appears beneath the climate velocity physics callback's automatic layer array;
it is a candidate for measurement, not evidence of a dominant allocation cost.

## Extended restart validation and a comparison limitation

Directory: `/private/tmp/wavetrisk-local-179c-restart-01`.
Unmodified and instrumented legacy were run separately to `time_end=0.3120`.
Both complete ten RK4 steps, remapping and checkpoint/restart; all printed
states agree. Their compressed checkpoint-4 files are byte-identical:

`e28987e6a8e7cd8ae60966e4628545907f998682132d3a41c2cbe2143e697538`.

However, comparing unchanged legacy against the earlier block extended run
reveals different post-restart evolution. The input text is identical, as is
the initial checkpoint-3 hash:
`8c7b651c2d6bc6f481838163591565771e1f9a95fd54ec9adf0ce9b74a117e2d`.
Printed states agree through the restart at 0.3104 d, including 64,222 DOFs.
At 0.3116 d, legacy reports 64,019 DOFs and minimum relative mass 0.95642;
block reports 64,069 and 0.95648. At 0.3128 d they also differ in remap activity.
Their checkpoint-4 hashes differ; matching earlier printed values does not
prove full-field equality before the restart. The cause has NOT been diagnosed.
This is not caused by the new legacy instrumentation, whose own transparency
check passes. Nor do existing block-versus-its-Domain-oracle tests establish
exact equivalence to every algorithm/detail in authoritative main.

Keep strict comparison rejection. For a matched performance experiment, use a
separate frozen seed ending **before** checkpoint 4 (for example `time_end=0.3090`,
seven ordinary steps for this exact checkpoint/dt), and retain the extended
restart test as a separate correctness investigation. Do not silently compare
the different post-restart workloads or relax tolerances to make them match.
The original full-length Stage179b comparison may correctly reject such a run;
the cluster protocol needs this distinction before it is resumed.
