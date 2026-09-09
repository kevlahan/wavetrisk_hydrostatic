# Stage 173: reusable thermodynamic columns

Base: accepted Stage 172, commit `abab57a2`.
Accepted for checkpoint commit after 83-rank cluster validation.

## Cluster acceptance

RK3 and RK4 checked runs complete through checkpoint/restart with both
dynamics and adaptation oracles enabled. All 42 RK3 and 36 RK4 cache/final
comparison records report zero mismatches. RK4's printed mass, timestep and
DOF sequence matches the optimized production and profiling runs.

Production time is 22.414 s versus Stage 172's 25.123 s (10.8% lower).
Profiled time is 22.590 s versus 25.610 s (11.8% lower). These are encouraging
single-run comparisons, not controlled repeated paired measurements.
First-eight-step work counters equal the local counters below; inverse
gather/alias/scatter calls, messages and bytes remain unchanged. The retained
Domain tendency compatibility pass still executes 32 times and costs 3.8031 s
rank-average in that window.

An initial optimized executable failed with SIGILL at startup after being
built on a different EPYC CPU using `-march=native`. Recompiling for the
execution platform resolved it. Fresh build directories and compatible CPU
targets are required; that failure did not require a numerical source change.

## Scope and computation removed

This is the approved narrower computational stage, not removal of the entire
Domain mass/velocity dependency chain. Only `src/parallel_block_mpi.f90`
changes relative to Stage 172.

Previously each edge endpoint and physical level independently reconstructed
the entire column pressure, subtracted preceding layers, and evaluated Exner.
Stage 173 computes pressure, dynamic Exner and potential temperature once per
referenced column, then reuses the values for adjacent edges and levels:

- Native thermodynamic-gradient evaluation: one lazy column build per storage
  node per block kernel invocation, across all its patches and edges.
- Domain compatibility-residual packing: one lazy column build per geometry
  node per residual transaction, shared by local packing and outgoing payloads.

For Z layers, each column now needs Z pressure additions and Z-1 subtractions,
plus Z Exner exponentiations. The previous cost repeats Z+k-1 pressure terms
and one exponentiation for every endpoint lookup at layer k. This eliminates
quadratic-in-Z pressure work and repeated endpoint exponentiations in these
two paths. It does not eliminate the required Domain tendency pass itself.

Buffers retain capacity; validity is reset before every producer transaction
or native block invocation. Unreferenced columns are not evaluated. No values
are reused across RK stages, blocks, adaptation, or restart. Native workspace
reuse relies on the existing serial per-rank block-kernel loop; future threaded
kernel execution would need separate workspaces. There is no new communication
during an RK stage. Profiling adds one counter reduction at report time only.

## Numerical and oracle contract

The original accumulation and subtraction order is retained, including
`lower - 0.5*g*mass`. The existing hydrostatic Exner cache is deliberately not
substituted: its midpoint expression has different rounding. The original
slow Domain/reference reconstruction remains independent and unchanged.

With the dynamics oracle enabled, every built column's Exner and temperature
values are additionally compared bit-for-bit with the original calculation.
Existing full tendency, boundary, and adaptation checks retain their existing
tolerances. The production cache executes in both oracle and non-oracle runs.

## Profiling

Two new lines report global sums for native and residual column processing:

```text
thermodynamic native: columns lookups old/new pressure terms = ...
thermodynamic residual: columns lookups old/new pressure terms = ...
```

`old` is the equivalent pressure-addition/subtraction count for the original
algorithm at the observed lookups; `new` counts that work in built columns.
These are algorithmic work counters, not timings or complete floating-point
operation counts. They exclude independent oracle work. Old exponentiation
count is `lookups`; new count is `columns * Z`. Existing inclusive phase timers
must not be added together as disjoint costs. Transport volumes and Domain
compatibility writeback counts should remain unchanged from Stage 172.

## Transfer and fresh builds

Transfer the changed source and manifest (the cluster is assumed to have the
accepted Stage 172 dependency set):

```sh
scp /Users/kevlahan/wavetrisk_hydrostatic/src/parallel_block_mpi.f90 bbserv:~/wav/src/
scp /Users/kevlahan/wavetrisk_hydrostatic/docs/parallel-blocks-stage173.sha256 bbserv:~/wav/docs/
```

On the cluster:

```sh
cd ~/wav
sha256sum -c docs/parallel-blocks-stage173.sha256
```

Select RK3 using `timeint_type` in `src/shared.f90`, then build into fresh
directories:

```sh
make -j1 DEBUG=check PARAM=param_J5 TEST_CASE=climate \
  BUILD_DIR=build/stage173-rk3-check BIN_DIR=bin/stage173-rk3-check
```

Switch `timeint_type` back to RK4 before these builds:

```sh
make -j1 DEBUG=check PARAM=param_J5 TEST_CASE=climate \
  BUILD_DIR=build/stage173-rk4-check BIN_DIR=bin/stage173-rk4-check
make -j1 DEBUG=false PARAM=param_J5 TEST_CASE=climate \
  BUILD_DIR=build/stage173-rk4-opt BIN_DIR=bin/stage173-rk4-opt
```

These commands do not update `~/wav/bin/climateJ5`. Copy directly from the
intended build directory, and record the executable hash. Use the same input,
starting checkpoint, task count and allocation for comparisons.

## Complete settings for each cluster test

Run in the test directory. No new feature flag is needed.

RK3 debug oracle:

```sh
cp ~/wav/bin/stage173-rk3-check/climate ./climateJ5
sha256sum ./climateJ5
export WAVETRISK_VALIDATE_BLOCK_DYNAMICS=1
export WAVETRISK_VALIDATE_BLOCK_ADAPTATION=1
export WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS=0
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=0
set -o pipefail
srun -n 83 ./climateJ5 simple.in 2>&1 | tee RK3_debug_oracle.log
```

RK4 debug oracle:

```sh
cp ~/wav/bin/stage173-rk4-check/climate ./climateJ5
sha256sum ./climateJ5
export WAVETRISK_VALIDATE_BLOCK_DYNAMICS=1
export WAVETRISK_VALIDATE_BLOCK_ADAPTATION=1
export WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS=0
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=0
set -o pipefail
srun -n 83 ./climateJ5 simple.in 2>&1 | tee RK4_debug_oracle.log
```

RK4 production timing:

```sh
cp ~/wav/bin/stage173-rk4-opt/climate ./climateJ5
sha256sum ./climateJ5
export WAVETRISK_VALIDATE_BLOCK_DYNAMICS=0
export WAVETRISK_VALIDATE_BLOCK_ADAPTATION=0
export WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS=0
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=0
set -o pipefail
srun -n 83 ./climateJ5 simple.in 2>&1 | tee RK4_production.log
```

RK4 optimized profiling:

```sh
cp ~/wav/bin/stage173-rk4-opt/climate ./climateJ5
sha256sum ./climateJ5
export WAVETRISK_VALIDATE_BLOCK_DYNAMICS=0
export WAVETRISK_VALIDATE_BLOCK_ADAPTATION=0
export WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS=0
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=1
set -o pipefail
srun -n 83 ./climateJ5 simple.in 2>&1 | tee RK4_profile.log
```

Compare first-eight-step profiles separately from the post-restart window.
Use repeated paired Stage 172/173 unprofiled runs before claiming an 83-rank
speedup. This stage targets actual computation, but remaining inverse,
restriction, Domain compatibility, and topology costs still limit total speed.

## Local validation

Fresh RK3/RK4 `DEBUG=check` and RK4 `DEBUG=false` builds succeed with the
normal 132-column source limit. Checked builds retain bounds checks, signaling
NaN initialization, floating-point traps, and warnings as errors. The optimized
build retains the existing `air_temperature` descriptor warning in unchanged
`parallel_block.f90`; no warning suppression was added.

Both four-rank checked oracle runs complete one restart timestep, including
all three RK3/four RK4 stages and the final transform. The new bitwise column
checks and existing tendency/adaptation checks pass. Cache/final-shadow records
have zero mismatches. These are short local oracle tests, not full 83-rank
cluster acceptance runs.

Isolated local test directory: `/private/tmp/wavetrisk-stage173.fgCDIJ`.
Oracle logs are `rk3-oracle/RK3_debug_oracle.log` and
`rk4-oracle/RK4_debug_oracle.log`. The user's standard test directory and
`bin/climateJ5` are unchanged.

Checkpoint 4 is byte-for-byte identical to Stage 172 (`cmp` succeeds), with
SHA-256 `6086e1fab875369095b1dea9800af7a39261ab751d6f10916718ed3c081179f3`.
The optimized non-oracle run completes all ten local timesteps, including
adaptation, vertical remapping, checkpoint write/reload and two post-restart
steps. Its complete printed mass/time/dt/refinement-level/DOF sequence matches
Stage 172. Log: `restart/RK4_optimized_restart_profile.log`.

First-eight-timestep optimized four-rank work counters (30 physical layers):

| Path | Columns | Exner lookups before | Pressure terms before | Pressure terms now |
|---|---:|---:|---:|---:|
| Native | 1903592 | 258001920 | 11481085440 | 112311928 |
| Residual | 1681088 | 172001280 | 7654056960 | 99184192 |

The algorithmic pressure-term reductions are 99.02% and 98.70%; Exner
exponentiation reductions are 77.87% and 70.68%, respectively. These counts
describe removed source-algorithm work, not measured hardware instructions.
Gather/scatter/alias calls, messages and bytes exactly match Stage 172 in this
window, as do the 179 Domain compatibility writebacks.

The same window's complete-timestep rank-average time is 314.13 s versus the
archived Stage 172 run's 314.89 s. Block tendency/hydrostatic time is 98.603 s
versus 97.099 s. This unpaired comparison shows **no demonstrated overall
speedup**. The work reduction must not be presented as a corresponding runtime
reduction; the broader production dependency-chain work remains necessary.
The full ten-step profiled total is 374.04 s versus the earlier Stage 172
run's 366.49 s; these were not paired runs under controlled machine load.
