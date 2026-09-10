# Stage 179: attribute computational and representation overhead

Base: accepted Stage 178, commit `2eee9646`. Stage 179 is an uncommitted
profiling change, not another performance claim. Numerical kernels, data
layouts, transport order, validation checks and oracle tolerances are unchanged.

## Stage 179a: cluster compilation correction

The report header was 134 columns long. Split it into concatenated continued
strings; the printed header is identical. The local GNU Fortran 15.2 compiler
accepted the original line with default flags, but rejected it when explicitly
given `-ffree-line-length-132`, reproducing the cluster error. The earlier claim
that local builds enforced 132 columns was incorrect.

The Makefile now explicitly selects `-ffree-line-length-132`. The profiling test
runner also checks all Stage 179 source files for overlong non-comment lines.
The corrected MPI module compiles with explicit 132-column limits and the
checked warning/error flags; profiling/parser tests pass. No numerical changes
or new runtime test requirements are introduced by this formatting correction.
Transfer `src/parallel_block_mpi.f90` and `Makefile` again, then resume the
planned cluster tests. The archive and checksum manifest have been updated.

## Questions this stage measures

1. How much reconstruction time is source serialization proof, received-block
   serialization proof, local-copy byte comparison, real deep copying, actual
   migration packing/unpacking, or MPI?
2. How much inverse time is plan construction versus repeated scaffold staging,
   packing, local aliases, posting/self copies, completion wait, installation
   and scalar/inner/outer kernels?
3. How much scalar time is setup/storage, geometry replication, initial
   transport/install, replay, or ghost communication?
4. How much shared-producer time remains in basic operators, primitive setup,
   Domain boundary packing/local copies/install, native mass/velocity, or waits?
5. How large are the actual expanded final-owner scalar working records, rather
   than only the producer/transport buffers reported by previous stages?

These are hypotheses to measure, not proof that every listed cost is removable.
No validation ablation is enabled in this stage. In particular, all previously
unconditional pack/unpack/byte-comparison checks still execute. Once their cost
is known, a separately validated, single-change ablation can test causality.

## Instrumentation contract

New switch: `WAVETRISK_PROFILE_BLOCK_DETAIL=1`, with
`WAVETRISK_PROFILE_PARALLEL_BLOCKS=1`. Both default off. The new switch accepts
only 0/1; rank disagreement is rejected once at initialization. No added
barriers or collectives occur inside a numerical/communication phase. Four
gathers run only at the existing end-of-window profiling report.

The new report has a stack of nested regions. **Self wall and process-CPU times
exclude nested instrumented regions.** Inclusive times overlap and must not be
summed. Recursive scopes are supported; their inclusive time also overlaps.
The report gives self wall average/minimum/maximum, average self CPU, CPU on the
maximum-wall rank, inclusive wall average, that rank's ID, and global call count.
It also prints all ranks' self-wall vectors, indexed by region ID, so averages
do not hide concentration of packing/kernel work on particular owners.
Calls count region entries across ranks, not MPI messages or completed physical
timesteps. Empty local regions and recursive entries can contribute calls.

Residual rows contain work not split into a more specific instrumented region.
For example, the dynamics residual still includes some RK/wavelet/hydrostatic
work and independent oracle work. Consult the existing inclusive phase report
and CPU samples for those residuals. Domain writeback's detail row primarily
measures its packing/exchange entry point; later Domain staging/commit remains
in enclosing residuals. The existing inclusive writeback row retains full scope.

There are no per-node clocks. Disabled hooks perform no clocks, allocations,
environment reads or reductions, but their small function-call overhead must
still be checked against the Stage 178 binary. Enabled clock/counter overhead
and rank perturbation must likewise be measured, not assumed negligible.

Process CPU includes MPI spinning and library work. Wall-minus-CPU is not a
network-cost estimate. A large MPI wait can mean late arrival of a rank doing
more copying/computation. The maximum rank can differ between phases; do not
sum per-phase maxima or treat a window's maximum as a reconstructed critical path.

Startup is excluded from the first detail timestep window. Work between
timestep scopes can still be attributed to its own instrumented region. The
analysis helper reports any difference between summed self times and the
timestep inclusive total rather than hiding it by normalization.

### Column geometry and storage

Three high-water samples cover the allocated 50-double scalar records in
patch, boundary and ghost storage:

- actual expanded record bytes;
- the 33 geometry/metadata slots within those records;
- a shared-once estimate dividing those geometry slots by scalar count and
  field-level count, retaining horizontal halo copies.

This is neither total RSS nor a guarantee all 33 slots can be removed without
new metadata. Summed rank peaks need not be simultaneous. The estimate includes
soil/padded field slots in the present allocation and must not be confused with
30 atmospheric layers alone. Shared horizontal geometry means physical vertical
layers share a stencil; different horizontal refinement levels do not share
the same grid. Local and received geometry-copy counters are separate and count
actual assignment extents, not hardware memory traffic.

## Files to transfer

From the repository root:

```sh
scp Makefile bbserv:~/wav/
scp src/parallel_block_profile.f90 src/parallel_block.f90 \
  src/parallel_block_build.f90 src/parallel_block_inverse.f90 \
  src/parallel_block_mpi.f90 src/comm_mpi.f90 src/multi_level.f90 bbserv:~/wav/src/
```

Keep an independently named, optimized Stage 178 executable before building.
The new Makefile must accompany the new module. A source manifest is provided
in `docs/parallel-blocks-stage179.sha256`. No changes to `time_integr.f90`,
`ops.f90` or physics sources are needed.

To transfer the optional analysis/sampling helpers (the cluster's `~/wav/test`
directory must exist):

```sh
scp -r test/parallel_block_profile bbserv:~/wav/test/
```

## Builds and numerical checks

Use fresh build directories and `make -j1`. Compile on a CPU compatible with
the execution partition: optimized builds use `-march=native`. Select RK3
in `src/shared.f90` for its build, then restore RK4 for both RK4 builds.

```sh
# timeint_type = "RK3"
make -j1 DEBUG=check PARAM=param_J5 TEST_CASE=climate \
  BUILD_DIR=build/stage179-rk3-check BIN_DIR=bin/stage179-rk3-check
# timeint_type = "RK4"
make -j1 DEBUG=check PARAM=param_J5 TEST_CASE=climate \
  BUILD_DIR=build/stage179-rk4-check BIN_DIR=bin/stage179-rk4-check
make -j1 DEBUG=false PARAM=param_J5 TEST_CASE=climate \
  BUILD_DIR=build/stage179-rk4-opt BIN_DIR=bin/stage179-rk4-opt
```

All five environment settings should be explicit for each run:

| Run | DYNAMICS | ADAPTATION | DETAILED_DIAGNOSTICS | PROFILE_PARALLEL_BLOCKS | PROFILE_BLOCK_DETAIL |
|---|---:|---:|---:|---:|---:|
| RK3 checked oracle | 1 | 1 | 0 | 1 | 1 |
| RK4 checked oracle | 1 | 1 | 0 | 1 | 1 |
| RK4 optimized timing | 0 | 0 | 0 | 0 | 0 |
| RK4 existing profile | 0 | 0 | 0 | 1 | 0 |
| RK4 detailed profile | 0 | 0 | 0 | 1 | 1 |

The full variable names for the oracle rows are:

```sh
export WAVETRISK_VALIDATE_BLOCK_DYNAMICS=1
export WAVETRISK_VALIDATE_BLOCK_ADAPTATION=1
export WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS=0
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=1
export WAVETRISK_PROFILE_BLOCK_DETAIL=1
set -o pipefail
# From your established restart test directory; copy only after build success.
cp ~/wav/bin/stage179-rk3-check/climate ./climateJ5
srun -n 83 ./climateJ5 simple.in 2>&1 | tee RK3_debug_oracle.log
cp ~/wav/bin/stage179-rk4-check/climate ./climateJ5
srun -n 83 ./climateJ5 simple.in 2>&1 | tee RK4_debug_oracle.log
```

For the detailed optimized run:

```sh
export WAVETRISK_VALIDATE_BLOCK_DYNAMICS=0
export WAVETRISK_VALIDATE_BLOCK_ADAPTATION=0
export WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS=0
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=1
export WAVETRISK_PROFILE_BLOCK_DETAIL=1
cp ~/wav/bin/stage179-rk4-opt/climate ./climateJ5
set -o pipefail
srun -n 83 ./climateJ5 simple.in 2>&1 | tee RK4_detail.log
```

Then repeat with both profile switches zero for `RK4_production.log`, and with
only `WAVETRISK_PROFILE_PARALLEL_BLOCKS=1` for `RK4_profile.log`. Debug timings
are not production performance measurements. Require completed runs, zero
cache/final oracle mismatches, identical printed states and checkpoint payloads.

## Matched performance experiment

Preserve named legacy, Stage 178 and Stage 179 binaries plus their hashes,
source revisions, compiler versions/flags, input/checkpoint hashes, Slurm job
ID, allocated hostnames and rank/CPU binding. Use the same 83-task placement,
restart and numerical settings. Do not compare binaries compiled for different
EPYC instruction sets, or mix measurements from different busy allocations.

1. Within one allocation, run legacy / 178 / 179-detail-off, then reverse
   the order. Obtain at least three observations of each. All oracles off.
2. Pair 179-detail-off / 179-existing-profile / 179-detail-profile and reverse
   the order. This quantifies instrumentation overhead and variability.
3. Compare ordinary matching timestep records separately from the checkpoint/
   restart spike. Also inspect dynamics, adaptation and output separately.
   Never subtract inclusive overlapping phase rows to invent an exclusive sum.
4. Once the detail profile identifies hot ranks/regions, sample the same ranks
   in legacy and block binaries as described below. A selected rank is not a
   complete job profile; include rank 0, a hot rank, and a contrasting owner.

`test/parallel_block_profile/analyze.py` reads logs and writes JSON to stdout:

```sh
python3 test/parallel_block_profile/analyze.py RK4_production.log RK4_profile.log RK4_detail.log
```

It compares printed state/mass records without CPU fields, separates records
preceded by checkpoint writes, lists every timing sample and ranks exclusive
regions. Medians are descriptive; no significance/speedup is inferred from a
single short run. Checkpoint hash equality is an additional independent check.

## CPU sampling and hardware counters (cluster, optional)

The provided `sample_rank.sh` wraps selected **application ranks**, not `srun`.
It requires the cluster's `perf` tool and permission to use it. It does not
change kernel settings, ask for sudo, or silently downgrade failed sampling.
Keep all application profiling/oracle switches off for this experiment.

After copying `test/parallel_block_profile/` to the cluster, choose ranks based
on the detailed report. Example (0,42,82 are examples, not assumed hot ranks):

```sh
export WAVETRISK_VALIDATE_BLOCK_DYNAMICS=0
export WAVETRISK_VALIDATE_BLOCK_ADAPTATION=0
export WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS=0
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=0
export WAVETRISK_PROFILE_BLOCK_DETAIL=0
export WAVETRISK_SAMPLE_RANKS=0,42,82
export WAVETRISK_SAMPLE_MODE=record
export WAVETRISK_SAMPLE_DIR="$PWD/stage179-samples-01" # new shared directory
set -o pipefail
srun -n 83 sh ~/wav/test/parallel_block_profile/sample_rank.sh \
  ~/wav/bin/stage179-rk4-opt/climate simple.in 2>&1 | tee RK4_sample.log
perf report --stdio --no-children --sort comm,dso,symbol \
  -i "$WAVETRISK_SAMPLE_DIR/rank-42.record"
```

The default is low-frequency user-space CPU-clock sampling with DWARF stacks;
the existing optimized build already includes `-g`. Samples include startup,
I/O and child processes: keep executable/MPI symbols distinct from compression
children and do not equate whole-run sample percentages with ordinary-step
percentages. Inspect call paths into allocation, serialization, lookup helpers,
array-copy routines, native arithmetic and MPI progress. `perf` distinguishes
Self versus Children attribution; see its [official tutorial](https://perfwiki.github.io/main/tutorial/)
and [kernel tool tips](https://github.com/torvalds/linux/blob/master/tools/perf/Documentation/tips.txt).

For a separate counter run, set `WAVETRISK_SAMPLE_MODE=stat` and a new output
directory. Hardware event availability and access vary across chips/partitions;
retain any unsupported/multiplexed-event warnings. Compare instructions, cycles
and cache events on matched ranks, but do not interpret generic cache misses
as measured DRAM bytes. If profiling access is denied, return that error and
the application detail logs; do not change cluster security policy.

## Decision before another performance change

Return a ranked attribution: measured self time, imbalance, sampled hot call
paths, bytes/work counts, evidence of block-specific or duplicate work, and a
realistically removable share. Prefer a complete computational dependency or
layout replacement with a numerical test, not merely smaller MPI messages.

Candidate outcomes include gating serialization proofs behind an oracle;
sharing horizontal geometry in working storage and using column-oriented
operators; preserving unaffected blocks across adaptation; or removing an
owner/representation handoff. These remain candidates until measurements
select one. This stage alone cannot promise proximity to legacy speed.

## Local validation

Fixture root: `/private/tmp/wavetrisk-stage179.yxQh6r`. The user's standard test
directory and `bin/climateJ5` are untouched.

- Fresh serial GNU Fortran RK3/RK4 checked and RK4 optimized builds passed with
  the local compiler's default line width (see the Stage 179a correction above).
  Checked bounds/FPE/SNaN flags and warnings remain.
  Only the pre-existing optimized `air_temperature` descriptor warning occurs.
- Profiler tests pass at checked O0/O2: disabled mode, nested/recursive self
  wall/CPU conservation, counter accumulation, peak samples and resets.
- Existing inverse-route pack/install tests pass at checked O0/O2.
- Log parser tests pass, including checkpoint separation, exclusive totals,
  CPU-independent state extraction and failure detection. The Slurm sampler
  passes shell syntax and unselected-rank execution checks; Linux `perf` itself
  cannot be exercised on this macOS host.
- Final-source RK3/RK4 debug oracles pass all three/four substages plus final
  inverse; all printed state/mass records match Stage 178. Cache/final mismatch
  counts remain zero. Both runs exercise the new detailed profiling scopes.
- Final optimized, detailed-profile RK4 passes ten steps, remapping, checkpoint
  4 write/reload and post-restart steps. All printed state/mass records match
  Stage 178. Checkpoint 4 is byte-identical, SHA-256:
  `6086e1fab875369095b1dea9800af7a39261ab751d6f10916718ed3c081179f3`.
- Both optimized report windows contain all 42 regions. Summed self wall times
  agree with timestep inclusive totals within printed rounding (under 0.4 ms
  in the first window). Removed Domain velocity-source/gradient counters and
  scalar producer oracle-record counts remain zero in production.
- A separate optimized one-step RK4 run with both oracles and both profiling
  switches off passes and matches the corresponding Stage 178 state/mass
  records. It produces no detailed profiling report. This is a disabled-path
  correctness check, not a statistically controlled overhead measurement.

Two overlapping local jobs were deliberately stopped to relieve memory pressure
and rerun serially; their interrupted logs are not passing tests or timing data.
Local validation timings are not a cluster performance comparison. The first
serial optimized window samples 7.576 GB of expanded scalar working records,
of which 5.000 GB are the replicated geometry/metadata subset; the shared-once
estimate is 60.980 MB. This storage observation does not assign execution cost
or prove an equivalent allocation can simply be deleted.
