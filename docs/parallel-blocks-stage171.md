# Stage 171: native inverse-transform boundary dependency chain

Base: accepted Stage 170, commit `c5360456`.
Stage 171 accepted for checkpoint commit on 2026-09-09 after cluster validation.

## Cluster acceptance (83 ranks)

The RK3 and RK4 debug-oracle logs completed with 42 and 36 cache/final-shadow
records respectively, all with zero reported mismatches. Corrected optimized
logs supplied after a filename mix-up confirm the Stage 171 production path.
Both optimized runs completed through checkpoint/restart; printed mass,
timestep, refinement-level and DOF sequences match the RK4 debug oracle.

In the matching first-eight-timestep profile, all four old inverse Domain
bridge rows are absent. Overall compatibility writeback calls fall from 499
to 179. Native gather/aliases/scatter are present with respectively 496/568/568
calls, 67022/412084/76751 messages, and
11456052480/7694283360/29114847360 reported bytes across ranks.

Unprofiled runtime is 25.284 s versus the historical Stage 170 result of
24.993 s; the corrected profiled run is 25.381 s. The user reports slower
cluster conditions, so these unpaired measurements do not establish either
a speedup or a meaningful regression. This is an accepted correctness and
architectural checkpoint, not a demonstrated cluster performance win.
The new setup and repeated native transfers remain the next optimization
target. Timing rows are inclusive and must not be added as disjoint costs.

## Production contract

The inverse transform no longer publishes complete Domain fields between
phases, invokes a Domain boundary callback, or imports the result back into
blocks. This includes compressed-wavelet boundary setup, scalar lifting and
reconstruction boundaries, outer-vector/pentagon reconstruction, inner-vector
boundaries, the provisional RK handoff, and the final transform handoff.

The new `parallel_block_inverse_mod` compiles a sparse dependency workspace,
signed alias schedules and per-level outer-edge operation lists once per
writeback-plan generation. Scalar and inner-vector arithmetic remains on
final-owner blocks. Required interior values come directly from those blocks;
native boundary values and outer-edge results return through sparse routes.
Numeric transfers use persistent peer schedules and reusable buffers. There
is no new per-phase collective. Metadata discovery is confined to plan setup.

The outer-edge operation lists preserve plus-side writes, source-expression
grouping, in-place producer order and pentagon corrections. Local aliases
retain the legacy **all-level** copy order. Remote aliases use the requested
level range and are packed before local copies. These distinctions must not
be replaced by a generic nearest-owner halo copy.

The fixed coarse input is also explicit: the old writeback copied
uncatalogued interior patches from **global `sol`**, even when its destination
was a provisional RK scaling array. Stage 171 captures that fixed input once
per transform and restores it at the equivalent gather points. Existing
boundary/scaffold values are seeded once; Domain topology and geometry are
read during compilation. There is no duplicated Domain inverse computation
on the production path.

One final full interior publication remains for compatibility consumers.
Native auxiliary/boundary values are also published at the consumer handoff.
The inverse's production writeback count is checked to be **exactly one**,
instead of `4 * (jmax-jmin) + 2`. The provisional all-level boundary operation
is now native. Its caller does not invoke `update_bdry` or re-import boundaries.
The final handoff retains native solution/wavelet interiors and boundaries;
it refreshes ghosts and derived hydrostatic state without full Domain imports.
Redundant ghost-verification exchanges run only for the oracle on this path.

## Oracle independence

With the dynamics oracle enabled, the original writeback, geometric boundary
callbacks and Domain-array outer-edge routine execute as references. Native
wavelet setup and each scalar/vector boundary phase are compared against
them before the native results are consumed. Outer/pentagon coverage counts
must agree. Existing tendency, boundary, hydrostatic and adaptation checks
remain enabled. Oracle tolerances have not been relaxed to admit a mismatch.

The first local development test exposed the fixed-coarse input distinction
above. Its mismatch was corrected by matching the input contract, not by
changing the comparison tolerance.

## Remaining limitations

This is not removal of all Domain usage in WAVETRISK. Forward-transform,
mass/velocity compatibility, physics, adaptation and IO still have their
existing interfaces. The sparse outer/boundary workspace executes on the
original geometric route owners; full redistribution of those operations to
final block owners is not part of this stage. Native dependency gathering,
alias communication and scattering still cost time. Their three new profile
rows expose that replacement cost rather than treating the deleted Domain
bridge time as a predicted speedup.

## Files to transfer

Transfer `Makefile` to the repository root, and these five source files to
`src/` (the new source module is mandatory):

```sh
scp /Users/kevlahan/wavetrisk_hydrostatic/Makefile bbserv:~/wav/
scp \
  /Users/kevlahan/wavetrisk_hydrostatic/src/parallel_block.f90 \
  /Users/kevlahan/wavetrisk_hydrostatic/src/parallel_block_inverse.f90 \
  /Users/kevlahan/wavetrisk_hydrostatic/src/parallel_block_mpi.f90 \
  /Users/kevlahan/wavetrisk_hydrostatic/src/adapt.f90 \
  /Users/kevlahan/wavetrisk_hydrostatic/src/time_integr.f90 \
  bbserv:~/wav/src/
```

Do not reuse a build directory from an older stage. RK3/RK4 selection remains
the existing `timeint_type` setting in `src/shared.f90`; it is not an oracle
environment variable. Select the intended scheme before building, and use
separate fresh build directories for each scheme/build mode. For example,
after selecting RK4:

```sh
make -j1 DEBUG=check PARAM=param_J5 TEST_CASE=climate \
  BUILD_DIR=build/stage171-rk4-check BIN_DIR=bin/stage171-rk4-check
cp bin/stage171-rk4-check/climate bin/climateJ5
```

For RK3 use the RK3 source setting and `stage171-rk3-check` directories. For
optimized RK4, select RK4 and use `DEBUG=false` with `stage171-rk4-opt`
directories. Copy the newly selected binary to the run directory every time.
Record source and executable hashes with the logs.
The delivery source manifest is `docs/parallel-blocks-stage171.sha256`;
from the repository root, check it with
`sha256sum -c docs/parallel-blocks-stage171.sha256` on the cluster or
`shasum -a 256 -c docs/parallel-blocks-stage171.sha256` on macOS.

## Complete test settings

Use the same restart input, 83 tasks and checkpoint as Stage 170. These four
variables are the complete settings for this stage; no new feature flag is
required.

| Run | Build | DYNAMICS | ADAPTATION | DETAILED_DIAGNOSTICS | PROFILE_PARALLEL_BLOCKS |
|---|---|---:|---:|---:|---:|
| RK3 oracle | `DEBUG=check`, RK3 | 1 | 1 | 0 | 0 |
| RK4 oracle | `DEBUG=check`, RK4 | 1 | 1 | 0 | 0 |
| RK4 production timing | `DEBUG=false`, RK4 | 0 | 0 | 0 | 0 |
| RK4 production profile | same optimized binary | 0 | 0 | 0 | 1 |

The full variable names and debug-oracle invocation are:

```sh
export WAVETRISK_VALIDATE_BLOCK_DYNAMICS=1
export WAVETRISK_VALIDATE_BLOCK_ADAPTATION=1
export WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS=0
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=0
cp ~/wav/bin/climateJ5 .
sha256sum ./climateJ5
set -o pipefail
srun -n 83 ./climateJ5 simple.in 2>&1 | tee RK4_debug_oracle.log
```

Use the RK3 binary and `RK3_debug_oracle.log` for the RK3 oracle. For production:

```sh
export WAVETRISK_VALIDATE_BLOCK_DYNAMICS=0
export WAVETRISK_VALIDATE_BLOCK_ADAPTATION=0
export WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS=0
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=0
cp ~/wav/bin/climateJ5 .
sha256sum ./climateJ5
set -o pipefail
srun -n 83 ./climateJ5 simple.in 2>&1 | tee RK4_production.log
```

Then reset the same test to the same starting checkpoint, set only
`WAVETRISK_PROFILE_PARALLEL_BLOCKS=1`, and run to `RK4_profile.log`. Do not use
the profiled run as the timing comparison. Debug oracles must pass before
acceptance; production must also complete through checkpoint/restart.

In the production profile, `inverse boundary writeback`, `inverse Domain
boundary callback`, `inverse full interior import` and `inverse boundary
import` should all be absent (zero calls). `inverse native gather`, `inverse native aliases`
and `inverse native scatter` report the new cost and peer payload volume.
Setup includes dependency-plan compilation and initial state seeding. Timers
remain inclusive and must not be added as disjoint costs.

Stage 170's 83-rank unprofiled reference was 24.993 s. Compare paired repeated
runs, ordinary-step timings and the matching eight-step profile window;
checkpoint/output costs and rank imbalance remain in the total.

## Local validation

Validation uses isolated inputs/outputs in
`/private/tmp/wavetrisk-stage171.ZHjzmv`. The user's standard restart directory
and `bin/climateJ5` are untouched. No 83-rank cluster execution is claimed.

Fresh GNU Fortran builds succeeded for checked RK3/RK4 and optimized RK4.
The checked builds use the normal strict flags, including 132-column source
limits, `-Werror`, bounds checks, signaling-NaN initialization and floating
point traps. No warning suppression or increased line limit was added.
The pre-existing optimized `air_temperature` descriptor warning is unchanged.

The checked RK3 and RK4 delivery oracle runs each completed one restart
timestep, including all RK sub-stages and the final fixed-coarse transform.
The checked RK4 non-oracle run also completed, with printed numerical output
matching its oracle run. Logs include
`rk3-oracle/RK3_delivery_oracle.log`,
`rk4-oracle/RK4_delivery_oracle.log` and
`rk4-production/RK4_checked_production.log`.

The optimized non-oracle run completed all ten local test timesteps, including
vertical remapping, topology changes, checkpoint 4 write/reload and subsequent
steps. Log: `restart/RK4_optimized_restart_profile.log`.
`cmp` confirms checkpoint 4 is **byte-for-byte identical** to Stage 170:

```text
6086e1fab875369095b1dea9800af7a39261ab751d6f10916718ed3c081179f3
```

The matching first-eight-step profile confirms 499 → 179 full compatibility
writeback calls, with all four old inverse bridge rows absent. The inverse
setup/seeding row increased from 3.4415 to 9.3686 rank-average seconds; the
new workspace and its setup are not free. The inverse row now includes the
provisional all-level closure that previously ran outside it, so it is not
an unchanged timing boundary.

The local profiled total was 303.66 s versus the archived Stage 170 run's
378.09 s (19.7% lower). These are historical single four-rank measurements,
not a controlled paired benchmark or a prediction of 83-rank speedup. The
cluster oracle, timing and profile runs above remain the acceptance tests.

Delivery binary SHA-256 (local macOS executables, not cluster binaries):

```text
22e79a71fde1bcbbb5528c93b8f4fd96da34d5462aa300beb3ebf06414b0bcef  stage171-delivery-rk4-check/climate
884a466831a0475a04afb70665f1ec331f4ac965d3cbc90e979058096bdade93  stage171-delivery-rk3-check/climate
6552d34d71f2fbc76f9615fa33f84d3638de2a1d0961389dcd2457980fe61e4d  stage171-delivery-rk4-opt/climate
```
