# Stage 172: resident, phase-specific inverse transport

Base: accepted Stage 171, commit `bea072fb`.
Stage 172 accepted for checkpoint commit on 2026-09-09 after cluster validation.

## Cluster acceptance (83 ranks)

User-supplied RK3 and RK4 debug-oracle logs complete successfully with 42 and
36 cache/final-shadow records respectively, all reporting zero mismatches.
Both optimized runs complete through checkpoint/restart, and their printed
mass/time/dt/refinement-level/DOF sequences match the RK4 debug oracle.

The matching first-eight-timestep profile confirms the Stage 172 transport:

| Metric | Stage 171 | Stage 172 |
|---|---:|---:|
| Gather bytes | 11456052480 | 4835316480 |
| Scatter bytes | 29114847360 | 19475907360 |
| Combined gather/scatter bytes | 40570899840 | 24311223840 |
| Alias messages | 412084 | 273640 |
| Alias bytes | 7694283360 | 7694283360 |
| Overall compatibility writeback calls | 179 | 179 |

Combined traffic falls 40.1% and alias messages 33.6%. The 50% combined
traffic target is not met. The four inverse Domain bridge operations remain
absent. Unprofiled time is 25.123 s versus Stage 171's 25.284 s; the profiled
total is 25.610 s. Inverse rank-average time is 4.1847 s versus 4.2090 s.
These single-run differences do not establish a meaningful speedup under
variable cluster conditions. Acceptance is for correctness and reduced
transport work, not demonstrated execution-speed improvement.

The next substantial stage should target removal of the retained production
Domain tendency compatibility computation and its re-import/recomposition
dependencies. That pass still runs once per RK stage without the oracle;
its 3.9504 s rank-average inclusive time is not all redundant or automatically
removable. Required mass/velocity/physics inputs must be produced natively
before disabling their Domain producers. Success must be judged by removed
production kernel/traversal work and paired runtime, not bytes alone.

## Production change

The inverse boundary workspace now retains unchanged intermediate state.
Ordinary scalar and inner-vector synchronization gathers only the level that
the preceding block kernel can modify, plus interior alias destinations that
must be restored from their authoritative block values. Before outer-vector
reconstruction and the provisional final closure, no block kernel has changed
the corresponding interiors since their last synchronization: these gathers
only restore alias-overwritten interiors. Fixed-coarse inputs are still reset
at the same points as Stage 171.

Scatter uses separate output plans. Ordinary phases install the active level,
all-level local alias destinations and fixed-coarse dependencies. The outer
phase installs both coarse and fine levels, including off-level plus-side and
pentagon targets. This avoids resending unrelated remote boundary levels.
Initial seeding and final compatibility handoffs remain complete.

The subsets are conservative: they retain entire scalar/vector node records
on selected levels, not just individual numerically changed entries. Local
aliases are deliberately not reduced to the remote phase's level range.
This is a safe level-based invalidation contract, not floating-point value
comparison, a generic halo substitution, or speculative stale-cache reuse.

Both ends compile matching persistent subplans from metadata appended to the
existing dependency-key exchange (six integers instead of four). There is no
additional count collective or per-phase route discovery. Plans and numeric
buffers are rebuilt only with the existing topology generation. All transfer
subplans share one maximum-sized numeric buffer pair; subplans retain only
integer routing metadata. Alias subplans likewise share one numeric buffer
pair, including combined multi-level exchanges. The increased metadata has a setup/memory cost to
measure without multiplying numeric buffer storage by the number of phases.

Multi-level alias operations now send one combined payload per peer. All
remote values are snapshotted before ordered local copies, and a precompiled
installation permutation preserves the previous level-major receive order.
Single-level operations and dependencies between reconstruction phases are
not fused. The outer-vector operation tape and numerical kernels are unchanged.

## Correctness contract

The same native transport executes with or without the oracle. Existing
independent Domain phase comparisons still check the entire workspace,
including values retained rather than transferred. Existing tolerances are
unchanged. Mass/velocity compatibility, scalar restriction, forward transforms
and the single final inverse compatibility writeback are outside this change.
The four inverse Domain bridge rows removed by Stage 171 must remain absent
in non-oracle profiles.

## Transfer and build

Only these source files differ from Stage 171:

```sh
scp \
  /Users/kevlahan/wavetrisk_hydrostatic/src/parallel_block_inverse.f90 \
  /Users/kevlahan/wavetrisk_hydrostatic/src/parallel_block_mpi.f90 \
  bbserv:~/wav/src/
scp \
  /Users/kevlahan/wavetrisk_hydrostatic/docs/parallel-blocks-stage172.sha256 \
  bbserv:~/wav/docs/
```

The delivery archive also includes Stage 171's unchanged Makefile,
parallel_block.f90, adapt.f90 and time_integr.f90 to permit verification of
the complete dependency set. Do not substitute files from an earlier stage.

On the cluster:

```sh
cd ~/wav
sha256sum -c docs/parallel-blocks-stage172.sha256
```

Select RK3 or RK4 using the existing `timeint_type` in `src/shared.f90` before
each build. Use distinct fresh build directories. For the RK3 debug oracle,
after selecting RK3:

```sh
make -j1 DEBUG=check PARAM=param_J5 TEST_CASE=climate \
  BUILD_DIR=build/stage172-rk3-check BIN_DIR=bin/stage172-rk3-check
```

After switching back to RK4:

```sh
make -j1 DEBUG=check PARAM=param_J5 TEST_CASE=climate \
  BUILD_DIR=build/stage172-rk4-check BIN_DIR=bin/stage172-rk4-check
make -j1 DEBUG=false PARAM=param_J5 TEST_CASE=climate \
  BUILD_DIR=build/stage172-rk4-opt BIN_DIR=bin/stage172-rk4-opt
```

Copy directly from the intended build directory into the test directory;
these commands do not update `~/wav/bin/climateJ5` automatically. Record the
copied binary's hash. Use the same input, starting checkpoint and 83 ranks
for all four runs, resetting the starting checkpoint/input between runs.

## Complete oracle settings and commands

Run these in the cluster test directory. No new feature flag is required.

RK3 debug oracle:

```sh
cp ~/wav/bin/stage172-rk3-check/climate ./climateJ5
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
cp ~/wav/bin/stage172-rk4-check/climate ./climateJ5
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
cp ~/wav/bin/stage172-rk4-opt/climate ./climateJ5
sha256sum ./climateJ5
export WAVETRISK_VALIDATE_BLOCK_DYNAMICS=0
export WAVETRISK_VALIDATE_BLOCK_ADAPTATION=0
export WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS=0
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=0
set -o pipefail
srun -n 83 ./climateJ5 simple.in 2>&1 | tee RK4_production.log
```

RK4 production profiling (same optimized executable):

```sh
cp ~/wav/bin/stage172-rk4-opt/climate ./climateJ5
sha256sum ./climateJ5
export WAVETRISK_VALIDATE_BLOCK_DYNAMICS=0
export WAVETRISK_VALIDATE_BLOCK_ADAPTATION=0
export WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS=0
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=1
set -o pipefail
srun -n 83 ./climateJ5 simple.in 2>&1 | tee RK4_profile.log
```

## Performance acceptance

Compare the first-eight-timestep window separately from the post-restart
window. Stage 171's 83-rank baseline is:

| Native phase | Calls | Messages | Bytes |
|---|---:|---:|---:|
| Gather | 496 | 67022 | 11456052480 |
| Aliases | 568 | 412084 | 7694283360 |
| Scatter | 568 | 76751 | 29114847360 |

The target is at least 50% lower combined gather/scatter bytes (baseline
40570899840), plus fewer alias messages. Alias bytes should be unchanged:
batching changes message granularity, not the requested alias values.
The phase call counts need not decrease. Monitor inverse setup/import and
total inverse cost as well as the overall dynamics and timestep costs;
timers are inclusive and must not be added as disjoint contributions.

This target is not a claimed result. Four-rank local and 83-rank cluster
traffic differ because local alias destinations remain mandatory and route
ownership differs. Use repeated paired Stage 171/172 timing runs within the
same allocation before claiming a speedup under variable cluster load.

## Local validation

Fresh delivery builds succeeded for RK3/RK4 `DEBUG=check` and RK4
`DEBUG=false`, using the normal strict 132-column source limit, bounds checks,
signaling-NaN initialization and floating-point traps in checked builds.
No warning suppression or relaxed oracle tolerance was added. The existing
optimized `air_temperature` descriptor warning is unchanged. An intermediate
incremental build rejected a stale Fortran module after a derived-type layout
change; all delivery binaries were subsequently built in fresh directories.

Both final delivery debug oracles completed one restart timestep, including
all three RK3/four RK4 stages and the final transform. Cache/final-shadow
records report zero mismatches; full-workspace inverse phase comparisons and
the existing tendency/adaptation assertions passed. Printed mass, timestep,
refinement level and DOF match the corresponding Stage 171 local oracle runs.
Logs: `rk3-oracle/RK3_final_oracle.log` and `rk4-oracle/RK4_final_oracle.log`.

The final optimized non-oracle run completed all ten local test timesteps,
including vertical remapping, topology changes, checkpoint 4 write/reload
and two subsequent steps. The full printed mass/time/dt/refinement-level/DOF
sequence matches Stage 171. Log: `restart/RK4_optimized_restart_profile.log`.
`cmp` confirms checkpoint 4 is byte-for-byte identical to Stage 171, with
SHA-256 `6086e1fab875369095b1dea9800af7a39261ab751d6f10916718ed3c081179f3`.

Matching first-eight-timestep local profile:

| Metric (four ranks) | Stage 171 | Stage 172 |
|---|---:|---:|
| Gather bytes | 1000546560 | 427119360 |
| Scatter bytes | 2403561600 | 1972149120 |
| Combined gather/scatter bytes | 3404108160 | 2399268480 |
| Alias messages | 8844 | 5568 |
| Alias bytes | 1692104640 | 1692104640 |
| Overall compatibility writeback calls | 179 | 179 |

Gather traffic is 57.3% lower; scatter traffic is 17.9% lower; combined
traffic is 29.5% lower. Alias messages are 37.0% lower with exactly unchanged
alias bytes. All four removed inverse Domain bridge rows remain absent.
The 50% combined-traffic target has **not** been achieved on four ranks;
the 83-rank result remains to be measured. All-level local alias destinations
and complete initial/final handoffs intentionally limit the current reduction.

Inverse setup/import rank-average time was 9.1882 s versus Stage 171's
9.3686 s. The full profiled run took 366.49 s versus the archived Stage 171
303.66 s. These unpaired historical timings do not demonstrate a speedup;
the higher total must not be presented as a performance win or assumed to
be entirely due to load. Cluster acceptance requires the four test runs and
paired timing comparison above. No 83-rank execution is claimed locally.

Test root: `/private/tmp/wavetrisk-stage172.Q0eHeE`.
The user's standard test directory and `bin/climateJ5` are not modified.

Delivery binary SHA-256 (local macOS executables, not cluster binaries):

```text
db87e9650906c46d6b5804c994f553f1c16e5570f78ae2227bdc98a0ba1bc1e0  stage172-delivery-rk3-check/climate
739cc0f1a7ec7b62e9292c099397ffc45dc089857b68229ccf9cfb15d27c33e8  stage172-delivery-rk4-check/climate
86bfe9272b7ef49adb710bdb1d4e329745254a40d34975e44e9f6ffa337199d9  stage172-delivery-rk4-opt/climate
```
