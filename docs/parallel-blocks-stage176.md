# Stage 176: direct shared-producer input for native scalar restriction

Base: validated Stage 175, commit `bfc02792`.
Accepted for checkpoint commit after local and 83-rank cluster validation.

## Cluster acceptance

RK3 and RK4 debug oracles complete through restart, with zero cache/final
mismatches in 42 and 36 reports. The three RK4 runs agree in all 20 printed
numerical state/mass records. Non-oracle profiles show zero interior oracle
records and zero Domain velocity source/gradient calls before and after restart.

Aggregate targeted scratch is 4157990400 -> 202557312 bytes (95.1% smaller);
post-restart it is 2514316800 -> 124822656 bytes (95.0% smaller). This is not
whole-application memory. Production time is 21.839 s versus Stage 175's
23.593 s; profiled time is 22.054 s versus 22.498 s. Both are lower, but these
single runs do not establish a repeatable speedup. Retained Domain tendency
compatibility is essentially unchanged at 3.6514 s rank-average over eight
steps, versus 3.6635 s. The next cut must replace the mass computational and
boundary/RK dependency chain rather than another producer-buffer refinement.

## Delivery boundary

This targets the repeated representation transfers in the shared primitive /
mass dependency chain. It is **not** a claim that the entire shared Domain
primitive or mass-boundary computation has been removed.

Previously, on every production plan rebuild, the shared producer populated
50-value interior records for every scalar and field level. Finalization read
those records back and packed 33 shared geometry values per node plus three
physics values per field. The completed Domain flux/divergence snapshots in
the intermediate records were not production inputs to the native kernels.

Now production uses:

`generation-scoped geometry + once-evaluated physics -> native transport stream -> native scalar restriction`

- Geometry is compiled once per writeback-plan generation, using the existing
  patch/ownership manifest. It contains no solution, mass flux or tendency.
- The shared producer writes live physics values directly into the actual
  transport layout, for both rebuilds and reused RK stages.
- Non-oracle producer-side full-record and repacking buffers have zero extent.
  The receive buffer is sized for the actual production stream, not 50 slots.
- Production no longer traverses interior Domain subtrees to capture direct
  flux, restricted flux, divergence and repeated geometry per physical field.
- Local records receive cached geometry directly. Existing final-owner kernel
  storage and scalar arithmetic remain unchanged.
- No additional setup collective is introduced. Message counts, bytes and
  the independent full-oracle wire contract are unchanged.

This removes allocation, initialization, traversal and copying work rather
than merely reducing network volume. The magnitude of any runtime improvement
must be measured; this does not eliminate the full retained tendency pass.

## Correctness and lifecycle

The producer keeps the original sender block order, including its existing
identity manifest. Field-major physics offsets account for the geometry record
inserted only at physical level 1 of the first scalar. Scaffold physics is zero
seeded each transaction. Every producing patch must refresh every physical
level before finalization; geometry coverage cannot substitute for this check.

Oracle mode retains the original complete capture. It independently compares
cached geometry at each actual level capture and compares the new outgoing
geometry/physics stream bit-for-bit with original records before MPI. No
tolerance is relaxed. The production path cannot call the oracle repacker.
Checked non-oracle execution additionally catches any accidental access to
the now zero-sized full-record buffers.

The existing generation invalidates geometry and address plans on ownership /
topology changes. Restart and adaptation must be exercised by acceptance tests.
Domain mass restriction, physical boundary closure, shared PV/pressure/KE,
mass RK compatibility and final-owner 50-slot kernel storage remain. They
are explicitly outside the producer-side intermediate removed here.

## Profiling acceptance

New rows:

```text
scalar producer: geometry nodes physics values oracle records = G P 0
scalar producer scratch bytes: old-equivalent current = OLD NEW
```

The first two counters count actual compiler/producer work; oracle records
must be zero in non-oracle runs. Scratch bytes are aggregate resident sizes at
the reporting point, not cumulative allocation volume or whole-process RSS.
`OLD` is the equivalent Stage 175 interior transport buffers for this plan;
`NEW` includes the new geometry cache as well as current transport buffers.
Oracle mode deliberately retains reference buffers and can use more memory.
Boundary buffers and final-owner kernel storage are excluded from both totals.

Stage 175's Domain velocity source/gradient counters must remain zero in
production. Scalar and inverse message counts/bytes should remain unchanged.
Use paired runs on the same CPUs and allocation, not debug/profile times as
performance comparisons. Inclusive phase timers must not be added as disjoint
budgets.

## Local validation

Artifacts: `/private/tmp/wavetrisk-stage176.0iX7hV`. Tests use four MPI ranks
and isolated copies of the established checkpoint-3 inputs; the user's
standard run directory and `bin/climateJ5` are untouched.

- Fresh serial checked RK3/RK4 and optimized RK4 builds: PASS. The final
  rebuild logs contain no new warnings. Checked flags retain the 132-column
  limit, warnings-as-errors, bounds checks and floating-point traps.
- Existing numeric velocity tests at checked `-O0` and `-O2`: PASS.
- Final checked RK4 with both oracles disabled: PASS, all four substages.
  Native velocity Domain-source/gradient calls are 0/0. Scalar producer
  counters are `43744 31495680 0`. Targeted scratch bytes are
  `324326400 27613440`, a 91.5% reduction including cached geometry.
- Final RK3 oracle: PASS, one timestep / three substages, including the new
  bitwise producer checks. Printed numerical state agrees with Stage 175.
- Final RK4 oracle: PASS, one timestep / four substages, with the same exact
  producer checks and matching Stage 175 printed numerical state.
- Full optimized RK4: PASS, ten timesteps including vertical remap,
  checkpoint/restart and two post-restart timesteps. All 22 printed numerical
  state/mass records match Stage 175. Checkpoint 4 is byte-identical, SHA-256
  `6086e1fab875369095b1dea9800af7a39261ab751d6f10916718ed3c081179f3`.
- First-eight-step producer counters: `358336 258001920 0`; post-restart:
  `71616 51563520 0`. Domain velocity source/gradient calls remain zero.
- Targeted scratch bytes in the two optimized profile windows:
  `360537600 -> 30282624` and `231436800 -> 21350784`, including geometry cache.
- Profiled call/message/byte records match the corresponding Stage 175 run;
  this stage changes producer memory work, not the communication contract.

These memory figures describe only the targeted buffers, not total RSS.
Checked execution time is not a performance benchmark.

### Matched local timing

The saved Stage 175 binary (SHA-256
`462213e349f8e58ad494e7f13b131ed5ff058af29daa81dc18d2f6308179c874`)
ran immediately after Stage 176 with the same input, checkpoint, four ranks
and profiling flags. All 22 printed numerical records and checkpoint 4 match.

Reported totals: **Stage 175 382.13 s; Stage 176 380.39 s** (0.46% lower).
This single profiled pair establishes **no significant overall speedup**.
First-eight-step scalar initial transport falls from 5.0404 s to 4.3921 s
rank-average (12.9%); other inclusive timings fluctuate. This stage is a
verified producer-memory and traversal reduction, not completion of the
shared primitive/mass computational replacement. Repeated unprofiled cluster
measurements remain necessary to assess execution-speed benefit.

Logs beneath the artifact directory:

```text
checked-production/RK4_checked_production.log
rk3-oracle/RK3_final_debug_oracle.log
rk4-oracle/RK4_final_debug_oracle.log
production/RK4_optimized_restart_profile.log
baseline/RK4_stage175_baseline_profile.log
```

Final binary hashes:

```text
b80e800567224299b4acde1870909ddda489eb95d61b5bb2cd048d9fc2053187  bin/stage176-rk3-check/climate
e54647aa0af1d48234d11ff05093e10e4eafd3f6af85d71017d6b4048fb399f0  bin/stage176-rk4-check/climate
467bb168aea390c7e26244ba92d927923bf93d04b43944221d866808905d476d  bin/stage176-rk4-opt/climate
```

## Transfer and builds

Only `src/parallel_block_mpi.f90` changes from Stage 175. Transfer it and the
manifest together; the remaining Stage 175 dependencies must already match.

```sh
scp /Users/kevlahan/wavetrisk_hydrostatic/src/parallel_block_mpi.f90 bbserv:~/wav/src/
scp /Users/kevlahan/wavetrisk_hydrostatic/docs/parallel-blocks-stage176.sha256 bbserv:~/wav/docs/
```

On the cluster, from `~/wav`:

```sh
sha256sum -c docs/parallel-blocks-stage176.sha256
```

Use **fresh serial build directories**, including all transitive Fortran module
dependencies. Build on CPUs compatible with every run node (`-march=native`).
Select RK3 in `src/shared.f90` before the first build, then restore RK4 for the
other two builds:

```sh
make -j1 DEBUG=check PARAM=param_J5 TEST_CASE=climate \
  BUILD_DIR=build/stage176-rk3-check BIN_DIR=bin/stage176-rk3-check
# Restore timeint_type to RK4 before continuing.
make -j1 DEBUG=check PARAM=param_J5 TEST_CASE=climate \
  BUILD_DIR=build/stage176-rk4-check BIN_DIR=bin/stage176-rk4-check
make -j1 DEBUG=false PARAM=param_J5 TEST_CASE=climate \
  BUILD_DIR=build/stage176-rk4-opt BIN_DIR=bin/stage176-rk4-opt
```

### RK3 debug oracle

```sh
cp ~/wav/bin/stage176-rk3-check/climate ./climateJ5
export WAVETRISK_VALIDATE_BLOCK_DYNAMICS=1
export WAVETRISK_VALIDATE_BLOCK_ADAPTATION=1
export WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS=0
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=0
set -o pipefail
srun -n 83 ./climateJ5 simple.in 2>&1 | tee RK3_debug_oracle.log
```

### RK4 debug oracle

```sh
cp ~/wav/bin/stage176-rk4-check/climate ./climateJ5
export WAVETRISK_VALIDATE_BLOCK_DYNAMICS=1
export WAVETRISK_VALIDATE_BLOCK_ADAPTATION=1
export WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS=0
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=0
set -o pipefail
srun -n 83 ./climateJ5 simple.in 2>&1 | tee RK4_debug_oracle.log
```

### RK4 optimized timing

```sh
cp ~/wav/bin/stage176-rk4-opt/climate ./climateJ5
export WAVETRISK_VALIDATE_BLOCK_DYNAMICS=0
export WAVETRISK_VALIDATE_BLOCK_ADAPTATION=0
export WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS=0
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=0
set -o pipefail
srun -n 83 ./climateJ5 simple.in 2>&1 | tee RK4_production.log
```

### RK4 optimized profile

```sh
cp ~/wav/bin/stage176-rk4-opt/climate ./climateJ5
export WAVETRISK_VALIDATE_BLOCK_DYNAMICS=0
export WAVETRISK_VALIDATE_BLOCK_ADAPTATION=0
export WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS=0
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=1
set -o pipefail
srun -n 83 ./climateJ5 simple.in 2>&1 | tee RK4_profile.log
```

These builds do not update `bin/climateJ5`; copy the explicit intended binary.
No new environment flag is required.

## Following substantial dependency replacement

The next complete computational cut must address the shared producer and mass
consumers together: native pressure/geopotential/PV/kinetic-energy and direct
mass flux, ordered coarse/fine mass restriction and divergence, boundary
completion, and native mass RK closure. Shared fields must be produced once
and consumed directly by scalar and velocity kernels. Moving only one callback
or reusing completed Domain mass tendencies would leave the dependency intact.
Geometry-owner execution is a supported starting point; final-owner placement
can follow once that full computational chain is independent.
