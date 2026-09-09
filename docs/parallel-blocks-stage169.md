# Stage 169: temperature-only production chain

Status: local and 83-rank cluster validation complete; checkpoint accepted.
Base: the cluster-tested, uncommitted Stage 168 sources on `5e1ffecc`.
This document accompanies the Stage 168/169 checkpoint commit.

## Cluster acceptance (2026-09-09)

RK3 and RK4 debug oracle logs complete through adaptation and restart, with zero
asserted cache/final-boundary mismatches. Non-oracle RK4 completes with matching
printed numerical signatures: 26.878 s without profiling and 27.257 s with
profiling. These single-run timings are not a controlled speedup measurement.
The matching eight-step profile window reports rank-average dynamics of
16.998 s versus 19.676 s in the preceding profile. Existing restriction transport
still reports 38,222,384,832 aggregate bytes, unchanged, with the new temperature
edge exchange reporting another 1,613,906,880 bytes. Production transport is the
next structural optimization target.

## What this changes

Non-oracle RK3/RK4 now execute this temperature chain:

1. Capture the once-evaluated physics residual and sparse boundary closures.
2. Compute interior advective flux on final-owner blocks.
3. Exchange signed native temperature edge fluxes through persistent per-level plans.
4. Restrict those fluxes, deliver native boundary divergence, and assemble temperature tendency.
5. Advance the native temperature RK state, including the extra plus-side boundary cells.
6. Publish temperature into the existing wavelet/halo input fields.

Production does not evaluate the duplicate Domain temperature advective-flux,
restriction, divergence or RK-tendency arithmetic. It also omits the unused full
Domain tendency-reference exchange. Mass and velocity compatibility remain.

The boundary part of step 5 is essential: the legacy RK kernel updates some
plus-side cells outside the catalogued patch interiors. Ordinary halo exchange
does not replace every one of these local producer values. Merely publishing
native interiors left stale temperature at the next wavelet stage. This caused
production-only drift even when the interior tendency oracle passed.

The new boundary plan records the exact producer patch, coordinates and level.
Native divergence is evaluated when that level's flux is final, not reconstructed
later from overwritten flux storage. Shared aliases use the last producer in
legacy traversal order. The compact results are returned to resident storage
and combined with the timestep-start state before an in-place final RK update.
An oracle-only assertion compares the independent legacy and native boundary
candidates after halo completion.

Topology plans are rebuilt only when the writeback generation changes and reused
through RK stages. Edge-request discovery is batched across levels, with no
per-level discovery collective. The RK-boundary adapter has one batched request
plan and one compact numerical exchange per stage.

## Remaining Domain use

The temperature physics interface still computes diffusion inputs once,
including `cal_Laplacian_scalars` for biharmonic diffusion, and supplies sparse
boundary geometry/state. This stage does not claim all physics or wavelet
computation is Domain-free. Domain-shaped prognostic staging, adaptation, and
mass/velocity compatibility remain.

Temperature publication uses existing scalar transport buffers, which carry
both scalar slots, but commits only physical temperature. It does not overwrite
mass, exchange velocity, or overwrite soil/surface temperature slots.

Profile rows `temperature edge plan`, `temperature edge exchange`, and
`temperature RK boundary` expose setup and communication costs. These are
inclusive child timers; do not add them to their enclosing phase totals again.
No 83-rank performance improvement is claimed before cluster measurement.

## Local validation

All runs use isolated directories, four MPI ranks, checkpoint 3, J5–J7 and
30 atmospheric levels. Short tests end at day 0.3010; the full restart test ends
at day 0.3120. MPI runs are serialized.

- Clean optimized non-oracle RK4: **PASS**, ten timesteps, adaptation, vertical
  remapping, checkpoint 4, restart and continuation. All printed mass, timestep,
  level and degree-of-freedom signatures match Stage 168.
- Checkpoint 4: **byte-identical to Stage 168**, SHA-256:
  `6086e1fab875369095b1dea9800af7a39261ab751d6f10916718ed3c081179f3`.
- First-timestep temperature physics, native tendency and boundary-flux
  fingerprints match the independent oracle in all four RK stages.
- Strict RK3/RK4 builds: pass bounds checks, signaling initialization,
  floating-point traps and warnings-as-errors compilation.
- Final RK4 checked oracle: **PASS**, all four stages and boundary-state comparisons.
- Final RK3 checked oracle: **PASS**, all three stages and boundary-state comparisons.
- Final RK4 checked non-oracle poisoned-workspace test: **PASS**, all four stages.
- Final RK3 checked non-oracle poisoned-workspace test: **PASS**, all three stages.

Existing tendency and integrated-state tolerances are not weakened. The older
observational boundary-owner/reconstruction counters can still report differences
for storage aliases; the asserted native comparisons and final-cache checks pass.

Full production log:
`/private/tmp/wavetrisk-stage169.1oG2Mz/final-restart/RK4_final_optimized_restart.log`.
This run used the exact clean delivery build, with both oracles, detailed
diagnostics and profiling disabled, and no diagnostic preprocessor macros.
Executable: `bin/stage169-final-release/climate`, SHA-256
`bfa38865dd7fe178b716e800f094f7f599c0009b8f91963c89ad74bdb8959037`.
An earlier fingerprint-enabled full run passed independently. Local timings
are not controlled performance measurements: only four ranks were used and
some earlier runs overlapped with serial compilation.

The fresh optimized build reports an existing `air_temperature` descriptor
warning in unchanged `src/parallel_block.f90`. The modified sources compile
without new warnings. The strict builds use `-Werror` and pass.

Final oracle logs:
`oracle/RK4_boundary_debug_oracle.log` and
`rk3-oracle/RK3_boundary_debug_oracle.log` under the development directory.
Final poisoned-workspace production logs:
`production/RK4_boundary_poisoned_production.log` and
`rk3-oracle/RK3_boundary_poisoned_production.log` in the same directory.

Development logs and the original Stage 168 source backups are under
`/private/tmp/wavetrisk-stage169.1oG2Mz`. Earlier failed development logs are
retained there; they must not be mistaken for delivery-test results. The standard
restart directory and `bin/climateJ5` have not been replaced.

Optional compile-time diagnostics, **not runtime oracle settings**:

- `WAVETRISK_TEST_TEMPERATURE_CUT`: poison Domain temperature tendency through
  native tendency/RK assembly; poison transport scratch after each diffusion
  prepass and verify that transport leaves it untouched. Restore physics scratch
  before the following prepass.
- `WAVETRISK_TEST_TEMPERATURE_FINGERPRINT`: print global fingerprints of physical
  temperature physics residuals, native divergence and boundary fluxes.

Do not enable either macro for performance measurements. Use separate build
directories when changing preprocessor flags; Make does not track flag changes.

## Transfer and build

Source SHA-256 values:

```text
cd250611bc357273ae8cd0b3680d0fe710a7a7ad56206819e5a98718925273f3  Makefile
c961f557de7f656a534e2fedf996f2c91bfce683bc329969fd432880c4e972db  src/ops.f90
e87e6e219ddf651eba52156b093d291a22637eef5a789d1b0e55f0312b2f2c68  src/multi_level.f90
9c38b841a192c7ec4630ab7845a913c2eabae4668381e29c8db58eaf37b0d960  src/parallel_block_mpi.f90
e847694a51cc70cfe10a1eb4bbea111dd4822efe71f6ad731eaefe60fc8cd755  src/time_integr.f90
```

Transfer all four source files and the Makefile dependency update:

```sh
scp \
  /Users/kevlahan/wavetrisk_hydrostatic/src/ops.f90 \
  /Users/kevlahan/wavetrisk_hydrostatic/src/multi_level.f90 \
  /Users/kevlahan/wavetrisk_hydrostatic/src/parallel_block_mpi.f90 \
  /Users/kevlahan/wavetrisk_hydrostatic/src/time_integr.f90 \
  bbserv:~/wav/src/
scp /Users/kevlahan/wavetrisk_hydrostatic/Makefile bbserv:~/wav/Makefile
```

Select `timeint_type="RK3"` or `"RK4"` in `src/shared.f90`; the input file
does not select the integrator. Use a fresh build directory per integrator/mode:

```sh
# With RK3 selected:
make -j1 DEBUG=check PARAM=param_J5 TEST_CASE=climate \
  BUILD_DIR=build/stage169-rk3-check && mv bin/climate bin/climateJ5

# With RK4 selected:
make -j1 DEBUG=check PARAM=param_J5 TEST_CASE=climate \
  BUILD_DIR=build/stage169-rk4-check && mv bin/climate bin/climateJ5

# With RK4 selected, optimized:
make -j1 DEBUG=false PARAM=param_J5 TEST_CASE=climate \
  BUILD_DIR=build/stage169-rk4-opt && mv bin/climate bin/climateJ5
```

After **every** build, copy the new executable into the run directory. Use
`&&`, not `;`, so a failed build cannot rename or execute a stale binary.

## Complete cluster test settings

| Run | Dynamics oracle | Adaptation oracle | Detailed diagnostics | Profiling |
|---|---:|---:|---:|---:|
| RK3 debug oracle | 1 | 1 | 0 | 0 |
| RK4 debug oracle | 1 | 1 | 0 | 0 |
| RK3 checked production | 0 | 0 | 0 | 0 |
| RK4 checked production | 0 | 0 | 0 | 0 |
| RK4 optimized timing | 0 | 0 | 0 | 0 |
| RK4 optimized profile | 0 | 0 | 0 | 1 |

For each RK3/RK4 debug-oracle run:

```sh
export WAVETRISK_VALIDATE_BLOCK_DYNAMICS=1
export WAVETRISK_VALIDATE_BLOCK_ADAPTATION=1
export WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS=0
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=0
cp ~/wav/bin/climateJ5 .
set -o pipefail
srun -n 83 ./climateJ5 simple.in 2>&1 | tee RK4_debug_oracle.log
```

Use the RK3 executable and `RK3_debug_oracle.log` for RK3.

For each checked-production run and the optimized timing run:

```sh
export WAVETRISK_VALIDATE_BLOCK_DYNAMICS=0
export WAVETRISK_VALIDATE_BLOCK_ADAPTATION=0
export WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS=0
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=0
cp ~/wav/bin/climateJ5 .
set -o pipefail
srun -n 83 ./climateJ5 simple.in 2>&1 | tee RK4_production.log
```

Use distinct log names for RK3, RK4 checked, and RK4 optimized runs. For a separate
optimized profile, change only `WAVETRISK_PROFILE_PARALLEL_BLOCKS=1`.

Both oracle and non-oracle tests are required: oracle success alone does not
prove that a removed production dependency is complete. Keep the same restart,
input, 83-rank placement and compiler options for timing comparisons. Include
remapping and checkpoint/restart, not just the first timestep.
