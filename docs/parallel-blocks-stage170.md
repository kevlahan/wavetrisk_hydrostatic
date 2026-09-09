# Stage 170: production scalar-restriction transport contract

Base checkpoint: `e01809db` (validated Stage 168/169 temperature chain).
Stage 170 accepted after 83-rank cluster validation (2026-09-09).

RK3/RK4 oracle runs completed with zero cache/final-boundary mismatches;
production and profiled RK4 completed with matching printed numerical results.
Unprofiled runtime was 24.993 s versus 26.878 s for Stage 169 (7.0% lower in
these single runs). The matching eight-step profile reports exactly
10,999,192,704 restriction bytes versus 38,222,384,832, a 71.2% reduction.
Temperature-edge traffic fell 36.3%; rank-average restriction time fell 24.2%.

## Scope

This stage removes oracle-shaped production **wire payloads**, deduplicates
native temperature edge traffic, and narrows temperature publication. It does
not change the numerical restriction algorithm, RK coefficients, boundary
producer ordering, final-owner mapping or oracle tolerances.

Only `src/parallel_block_mpi.f90` changes from the checkpoint.

### Production boundary contract

On each topology generation, the boundary stream contains 33 shared geometry /
mask / provenance slots **once per node**, then three live values per scalar
field: mass compatibility flux or temperature closure. Previously it carried
all 50 slots for every field. Subsequent RK stages carry just the three live
values. Capture writes the compact outgoing stream directly, eliminating the
old in-place full-record repack. The persistent send/receive allocations use
the reduced wire sizes. The independent oracle retains its original records.

### Production ghost contract

The generation's full ghost stream contains shared geometry once per node,
plus three native positive-edge fluxes and one native divergence per field.
Reference flux, reference divergence and other non-consumed field slots are
not transported. Missing slots are poisoned on installation so an accidental
production dependency fails instead of silently receiving zero. Local ghost
copies use the same codec as remote copies. Existing per-level dynamic flux /
divergence refreshes, including final-boundary cache service, are unchanged.

For F scalar/vertical fields, full per-node payload sizes are:

| Route | Stage 169 | Stage 170 |
|---|---:|---:|
| Boundary | 50F | 33 + 3F |
| Ghost | 50F | 33 + 4F |

For the current two-scalar, 30-atmospheric-level, 10-soil-level layout (F=82),
these full payloads shrink by 93.2% and 91.2%, respectively. This is a byte-count
reduction, not a promised timestep speedup. Existing interior physics transport
and dynamic ghost refreshes still contribute to the aggregate restriction row.

### Unique native temperature edge requests

Compile one request per destination rank and `(final block, patch, edge,
evaluation level)` key. A receiver-side fanout map preserves every boundary
destination and sign. Hash collisions are resolved by exact key comparison.
Level remains part of identity: different restriction stages must not share
temporally different producer values. Existing batched request discovery is
retained; no extra collective is introduced, and per-level numeric schedules
remain persistent within the topology generation. This routing is checked by
both the oracle and production paths.

### Temperature-only publication

The native RK publication route transmits only atmospheric temperature, not
mass, soil/surface slots or velocity. It uses the established patch manifest
and existing buffers; received temperature is expanded backwards into the
validated Domain staging addresses. The resident Domain staging interface and
mass/velocity compatibility computations are otherwise unchanged.

## What remains

The 50-slot local kernel workspace remains as an address-compatible storage
layout; this stage removes its redundant wire representation, not all local
storage. Shared geometry is expanded into that workspace on a generation
rebuild. The once-evaluated physics interface, mass/velocity computation, native
inverse-transform boundary round trip, adaptation and IO remain. No new
production fallback to Domain temperature arithmetic has been introduced.

## Validation

Local validation is complete in `/private/tmp/wavetrisk-stage170.fo7zyw`.
Runs use separate output directories and four MPI ranks; the user's standard
restart directory and `bin/climateJ5` are untouched. Checked oracle and production
tests cover one complete timestep for each integrator; the optimized RK4 test
covers the full restart regression. Local timing is not an 83-rank performance test.

The combined checked RK4 production smoke test completes all four stages.
Its native temperature edge traffic is 19,922,880 bytes, compared with
46,086,720 bytes in the same local run before edge deduplication (56.8% less).
The publication change reduces the inclusive writeback byte count from
211,893,248 to 206,408,704; most of that inclusive row is other writeback work.
Both comparisons use identical one-step topology. The first development run
already included the new boundary/ghost codec, so it is not a Stage 169 traffic
baseline for those two components.

Final optimized RK4 completes ten timesteps, two remaps, checkpoint 4, restart
and continuation. Every printed mass, timestep, level and degree-of-freedom
signature matches the Stage 169 local reference. Checkpoint 4 is byte-identical:
`6086e1fab875369095b1dea9800af7a39261ab751d6f10916718ed3c081179f3`.
Log: `restart/RK4_optimized_restart.log` in the development directory above.
Final optimized executable SHA-256:
`7203e634fe17d13852499d86e3d7e9ff8642b45d6530dac5d517f4effdec299a`.

This four-rank profiled run reports 378.09 s, versus 327.53 s for the earlier
unprofiled Stage 169 local reference. These are not controlled timing runs
(profiling settings differ and machine load was not controlled); **no local
speedup is claimed**.
The 83-rank communication-heavy workload must establish performance acceptance.

Final RK4 strict debug oracle: **PASS**, all four stages, 2,268,960 final-boundary
samples matched and zero cache/final mismatches. Log: `rk4-oracle/RK4_debug_oracle.log`.
Final RK3 strict production: **PASS**, all three stages; printed mass, timestep,
level and degree-of-freedom signatures match the Stage 169 RK3 oracle reference.
Log: `rk3-production/RK3_checked_production.log`.
Final RK3 strict debug oracle: **PASS**, all three stages, 1,701,720 final-boundary
samples matched and zero cache/final mismatches. Log: `rk3-oracle/RK3_debug_oracle.log`.
The RK3 production/oracle printed signatures agree exactly. Existing observational
owner/dscalar-shadow mismatch counters are not acceptance assertions; those can
remain nonzero just as at the validated checkpoint.

Final strict builds use `-Werror`, 132-column source limits, bounds checking,
signaling initialization and invalid/zero/overflow traps. The changed source
also compiles without optimized-build warnings. The fresh optimized build
reports the pre-existing `air_temperature` descriptor warning in unchanged
`src/parallel_block.f90`; no warning suppression was added.

## Transfer and clean builds

On a cluster already containing the Stage 169 checkpoint, transfer:

Source SHA-256:
`4fced299ef2561a00169aa3e52240c43096d722f0c5ad4c1759cf31c8faf127a`.

```sh
scp /Users/kevlahan/wavetrisk_hydrostatic/src/parallel_block_mpi.f90 bbserv:~/wav/src/
```

Select `timeint_type="RK3"` or `"RK4"` in `src/shared.f90`. Use a fresh build
directory for each integrator/mode: the changed Fortran plan type can make an
incremental build encounter stale transitive module metadata. Make does not
track compiler-flag changes either.

```sh
# RK3 selected in src/shared.f90:
make -j1 DEBUG=check PARAM=param_J5 TEST_CASE=climate \
  BUILD_DIR=build/stage170-rk3-check && mv bin/climate bin/climateJ5

# RK4 selected in src/shared.f90:
make -j1 DEBUG=check PARAM=param_J5 TEST_CASE=climate \
  BUILD_DIR=build/stage170-rk4-check && mv bin/climate bin/climateJ5

# RK4 selected, optimized:
make -j1 DEBUG=false PARAM=param_J5 TEST_CASE=climate \
  BUILD_DIR=build/stage170-rk4-opt && mv bin/climate bin/climateJ5
```

After every build copy the revised executable into the execution directory.
Use `&&` so a failed build cannot rename or run an earlier binary.

## Complete test settings

| Run | Dynamics oracle | Adaptation oracle | Detailed diagnostics | Profiling |
|---|---:|---:|---:|---:|
| RK3 checked oracle | 1 | 1 | 0 | 0 |
| RK4 checked oracle | 1 | 1 | 0 | 0 |
| RK3 checked production | 0 | 0 | 0 | 0 |
| RK4 checked production | 0 | 0 | 0 | 0 |
| RK4 optimized timing | 0 | 0 | 0 | 0 |
| RK4 optimized profile | 0 | 0 | 0 | 1 |

For each checked oracle run:

```sh
export WAVETRISK_VALIDATE_BLOCK_DYNAMICS=1
export WAVETRISK_VALIDATE_BLOCK_ADAPTATION=1
export WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS=0
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=0
cp ~/wav/bin/climateJ5 .
set -o pipefail
srun -n 83 ./climateJ5 simple.in 2>&1 | tee RK4_debug_oracle.log
```

Use the RK3 executable and `RK3_debug_oracle.log` for RK3. For checked production
and optimized timing, explicitly change both oracle variables to zero:

```sh
export WAVETRISK_VALIDATE_BLOCK_DYNAMICS=0
export WAVETRISK_VALIDATE_BLOCK_ADAPTATION=0
export WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS=0
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=0
cp ~/wav/bin/climateJ5 .
set -o pipefail
srun -n 83 ./climateJ5 simple.in 2>&1 | tee RK4_production.log
```

Use separate logs for each integrator/build. For a separate optimized profile,
change only `WAVETRISK_PROFILE_PARALLEL_BLOCKS=1` and use `RK4_profile.log`.
Do not enable the optional Stage 169 diagnostic compiler macros for timing.

Both checked production tests are essential: oracle runs retain the independent
full-record restriction transport. Include adaptation, remapping, checkpoint /
restart and continuation, rather than testing only the first stage.

## Performance acceptance

Compare equivalent profile windows, not the eight-step window against the
one-step post-restart window. Expect reduced bytes in `restriction initial full`
and `restriction ghost full`, reduced `temperature edge exchange` bytes, and
reduced temperature-publication contribution to `Domain compatibility writeback`.
Dynamic restriction refresh payloads and peer-message counts need not shrink.
Inclusive parent/child timers overlap and must not be added.

Applying only the wire-width changes to the preceding 83-rank eight-step
profile predicts aggregate restriction traffic of about 11.00 GB instead of
38.22 GB for identical topology/field counts (about 71% lower). This excludes
the separate deduplicated temperature edge and publication savings. Validate
the actual counts on the cluster, then compare repeated **unprofiled** timings
with the 26.878 s Stage 169 run under identical placement/input/compiler flags.
