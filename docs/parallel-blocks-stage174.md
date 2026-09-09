# Stage 174: native velocity-residual construction

Base: accepted Stage 173, commit `6694a7d9`.
Accepted for checkpoint commit after 83-rank cluster validation.
Only `src/parallel_block_mpi.f90` changes.

## Cluster acceptance

Both checked oracle runs complete through checkpoint/restart with dynamics
and adaptation validation enabled: 42 RK3 and 36 RK4 cache/final records,
all with zero mismatches. RK4's printed mass/time/dt/refinement-level/DOF
sequence matches both optimized runs.

The non-oracle profile confirms zero residual-column work before and after
restart. Raw values and native residuals both total 129000960 in the first
eight steps and 12827520 in the post-restart step. Native thermodynamic work
and inverse communication counts/volumes are unchanged from Stage 173.

Production total is 22.418 s versus Stage 173's 22.414 s; profiled total is
23.061 s versus 22.590 s. There is **no demonstrated cluster speedup**.
Acceptance is for correctness and elimination of duplicate production
residual computation, not for the unpaired local runtime reduction below.
The retained Domain tendency compatibility pass still runs 32 times in the
first-eight-step window, taking 3.8775 s rank-average.

## Production dependency removed

Stage 173 still built a second set of thermodynamic columns on geometry owners
to pack velocity residuals. Stage 174 removes that entire production residual
reconstruction and transports the already available raw compatibility velocity
tendency instead, using the same persistent route and payload size.

Previously:

`Domain trend + Domain-reconstructed Exner term -> transport -> subtract native Exner term`

Now:

`Domain trend -> transport -> add native Exner term -> subtract native Exner term`

The native Exner term is already required and computed by the block kernel;
forming the residual adds no new column integration or exponentiation there.
The add-then-subtract order is deliberately retained, not simplified to raw
trend, because the intermediate rounding is part of the accepted result.

This is removal of duplicate thermodynamic residual computation, **not**
removal of the full Domain velocity gradient or mass/velocity tendency pass.
In particular, restricted Domain Exner remains authoritative in that retained
gradient. No assumption that restricted Exner equals reconstructed Exner is
made. The validated mass, scalar restriction, physics and inverse paths remain
unchanged. No new production communication or feature flag is introduced.

## Independent oracle and validity

Both modes snapshot the same raw compatibility tendency before the complete
Domain oracle can overwrite it. Oracle mode additionally constructs and
transports the old Stage 173 residual into separate reference storage. Before
subtracting the native Exner term, the native residual must match this reference
**bit-for-bit**. This checks the exact intermediate whose equality guarantees
unchanged Stage 173 recomposition; final tendency tolerances are not loosened.
Original per-column bitwise checks and complete legacy tendency comparisons
remain enabled. A mismatch aborts with block/patch/node/edge/level and values.

Raw and reference payloads have separate readiness state; reference storage is
allocated only in oracle runs and cannot overwrite raw inputs. Both refresh
every transaction, including after adaptation/restart. The old non-compatibility
diagnostic path retains its original residual contract. There is no stale-cache
fallback. Oracle mode does incur an extra reference exchange by design.

## Work acceptance

Non-oracle profiles must contain:

```text
thermodynamic residual: columns lookups old/new pressure terms = 0 0 0 0
native velocity residual: raw values formed oracle-checked = N N 0
```

For the matching first-eight-timestep test window, Stage 173's residual path
built 1681088 columns, made 172001280 Exner lookups, and evaluated 50432640
Exner exponentiations plus 99184192 pressure-integration terms. All this
production residual-path work should disappear. Native thermodynamic work
should remain unchanged. Raw payload values and formed residual counts should
match. With profiling also enabled during oracle tests, the third counter
must equal formed residuals; residual columns are then expected for references.

Inverse gather/alias/scatter counts and bytes, and the 179 compatibility
writebacks in the first eight steps, should remain unchanged. Counters are
algorithmic work, not hardware instructions or a promised speedup. Compare
paired Stage 173/174 unprofiled runs on identical CPUs and allocations; inspect
the first-eight-step window separately from the post-restart window. Existing
timers are inclusive and cannot be added as disjoint costs.

## Transfer and build

```sh
scp /Users/kevlahan/wavetrisk_hydrostatic/src/parallel_block_mpi.f90 bbserv:~/wav/src/
scp /Users/kevlahan/wavetrisk_hydrostatic/docs/parallel-blocks-stage174.sha256 bbserv:~/wav/docs/
```

On the cluster, verify the accepted dependency set:

```sh
cd ~/wav
sha256sum -c docs/parallel-blocks-stage174.sha256
```

Build on a CPU compatible with **all** execution nodes: optimized builds use
`-march=native`. A binary built on a different EPYC generation may SIGILL.
Use fresh build directories rather than reusing architecture-specific objects.
Select RK3 using `timeint_type` in `src/shared.f90`, then:

```sh
make -j1 DEBUG=check PARAM=param_J5 TEST_CASE=climate \
  BUILD_DIR=build/stage174-rk3-check BIN_DIR=bin/stage174-rk3-check
```

Switch back to RK4 before:

```sh
make -j1 DEBUG=check PARAM=param_J5 TEST_CASE=climate \
  BUILD_DIR=build/stage174-rk4-check BIN_DIR=bin/stage174-rk4-check
make -j1 DEBUG=false PARAM=param_J5 TEST_CASE=climate \
  BUILD_DIR=build/stage174-rk4-opt BIN_DIR=bin/stage174-rk4-opt
```

These commands do not update `~/wav/bin/climateJ5`. Copy the intended binary
directly for each test. Keep the same starting checkpoint and input.

## Complete cluster test commands

Run in the test directory. RK3 debug oracle:

```sh
cp ~/wav/bin/stage174-rk3-check/climate ./climateJ5
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
cp ~/wav/bin/stage174-rk4-check/climate ./climateJ5
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
cp ~/wav/bin/stage174-rk4-opt/climate ./climateJ5
sha256sum ./climateJ5
export WAVETRISK_VALIDATE_BLOCK_DYNAMICS=0
export WAVETRISK_VALIDATE_BLOCK_ADAPTATION=0
export WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS=0
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=0
set -o pipefail
srun -n 83 ./climateJ5 simple.in 2>&1 | tee RK4_production.log
```

RK4 production profiling:

```sh
cp ~/wav/bin/stage174-rk4-opt/climate ./climateJ5
sha256sum ./climateJ5
export WAVETRISK_VALIDATE_BLOCK_DYNAMICS=0
export WAVETRISK_VALIDATE_BLOCK_ADAPTATION=0
export WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS=0
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=1
set -o pipefail
srun -n 83 ./climateJ5 simple.in 2>&1 | tee RK4_profile.log
```

## Local build validation

Fresh RK3/RK4 checked and RK4 optimized builds succeeded. A new optimized
warning in an oracle guard was corrected by evaluating the oracle-enabled
function in a separate conditional; the fallback diagnostic entry point was
also made to initialize its independent reference through the same wrapper.
All three delivery executables were rebuilt after these changes. Final rebuild
logs contain no warnings/errors; the initial full optimized build also shows
the previously documented `air_temperature` descriptor warning in unchanged
`parallel_block.f90`. No warning suppression or relaxed tolerance was added.

The normal 132-column source limit, bounds checking, signaling NaN initialization
and floating-point traps are retained in checked builds. `git diff --check`
and the source manifest pass. Local tests use isolated directories under
`/private/tmp/wavetrisk-stage174.j2yt7o`; the standard test directory and
`bin/climateJ5` are untouched.

Both final four-rank debug-oracle executables complete one restart timestep
with dynamics and adaptation validation enabled. Profiling is also enabled
locally to count the new exact checks:

| Test | Raw values | Native residuals | Bitwise matches |
|---|---:|---:|---:|
| RK3 | 11810880 | 11810880 | 11810880 |
| RK4 | 15747840 | 15747840 | 15747840 |

All cache/final-shadow records have zero mismatches; original tendency and
adaptation assertions also pass. Logs: `rk3-oracle/RK3_debug_oracle.log` and
`rk4-oracle/RK4_final_debug_oracle.log`. These short local oracle tests do not
replace the full 83-rank cluster oracle runs.

The optimized run writes and reloads checkpoint 4. `cmp` confirms that file
is byte-for-byte identical to Stage 173, with SHA-256
`6086e1fab875369095b1dea9800af7a39261ab751d6f10916718ed3c081179f3`.

First-eight-timestep non-oracle profile:

```text
thermodynamic native: columns lookups old/new pressure terms = 1903592 258001920 11481085440 112311928
thermodynamic residual: columns lookups old/new pressure terms = 0 0 0 0
native velocity residual: raw values formed oracle-checked = 129000960 129000960 0
```

This confirms complete removal of the production residual-column work, with
unchanged native-column work. Inverse gather/alias/scatter calls, messages
and bytes exactly match Stage 173: 496/806/427119360, 568/5568/1692104640,
and 568/1605/1972149120 respectively. The Domain compatibility writeback
phase still has 179 calls. The retained Domain tendency pass still has 32.

The first-eight-step rank-average time is 240.80 s, versus the earlier Stage
173 run's 314.13 s. This is an **unpaired local comparison**, not an established
83-rank speedup or a precise attribution of the timing difference. Cluster
validation and matched timing remain required.

The final optimized executable completes all ten local timesteps, including
vertical remapping, topology changes, checkpoint/restart and two post-restart
steps. Its entire printed mass/time/dt/refinement-level/DOF sequence matches
Stage 173. The post-restart window also has zero residual-column work and
25781760 raw values / 25781760 native residuals, with unchanged inverse
transport and compatibility-writeback counts. Full profiled total: 287.40 s
versus the earlier Stage 173 run's 374.04 s, subject to the same unpaired-run
caveat. Log: `restart/RK4_optimized_restart_profile.log`.
