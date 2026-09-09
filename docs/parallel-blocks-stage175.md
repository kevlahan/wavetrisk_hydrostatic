# Stage 175: native geometry-owner velocity computation

Base: Stage 174, `d75e4634`.
Accepted for checkpoint commit after local and 83-rank cluster validation.

## Cluster acceptance

The supplied RK3 and RK4 checked oracle logs complete through restart with
dynamics and adaptation oracles enabled. All 42 RK3 and 36 RK4 cache/final
comparison reports have zero mismatches. The RK4 checked, production and
profiled runs agree in all 20 printed numerical state/mass records.

Non-oracle profiles show zero Domain velocity source/gradient calls before
and after restart. Native publication counts are 129000960 and 12827520;
native RK values are 323910720 and 31723200. Native plan/source/gradient cost
is 0.0813394 s rank-average for the first eight steps.

Production total is 23.593 s versus recorded Stage 174's 22.418 s; profiled
total is 22.498 s versus 23.061 s. These opposite timing movements demonstrate
no established overall cluster speedup. Acceptance is for correctness and
the verified production dependency replacement. Further work should target
the shared primitive/mass chain and repeated representation transfers.

## Agreed delivery boundary

The user approved **native geometry-owner execution first**, deferring migration
of these kernels to final owners. This supersedes the placement requirement in
the original Stage 175 proposal, not the requirement to remove completed Domain
velocity tendencies as production inputs.

The production chain is now:

`shared primitives -> native Qperp/source -> native child-edge restriction -> native gradient -> native tendency publication/RK`

The existing shared mass-compatible pass still produces mass flux, PV,
pressure/geopotential, kinetic energy, Bernoulli and Exner once. Its arrays are
read through a zero-copy primitive handoff. This stage does **not** claim removal
of all shared Domain computation, all dual storage, or final-owner execution of
velocity kernels. The pressure-layer helper in the numeric module remains a
tested building block, not a second production pressure producer.

## Production changes

- `parallel_block_velocity.f90` contains Domain-free numerical kernels, ordered
  per-level operation programs and separately owned numerical/RK workspaces.
- `multi_level.f90` compiles addresses, static geometry weights, masks and
  direct/restriction schedules. Plans persist across RK substages. Ownership or
  topology generations invalidate them; exact mask-predicate comparisons also
  catch mask-only changes. No new MPI collective is needed for plan rebuilding.
- Qperp and the direct source execute on flat numerical arrays. The retained
  read-only physics interface is evaluated once per needed node/physical/grid
  level. Absent-child duplicate evaluations are collapsed; mixed-mask direct
  writes and child sums retain their original order. Untouched source scaffold
  values start at zero, as in the reference.
- The complete native gradient uses **final restricted** Bernoulli/Exner and
  the original edge signs, theta averages, metrics and arithmetic order. It
  does not substitute the separately reconstructed dynamic Exner cache.
- The existing owner-publication route carries native-workspace tendency,
  never `trend(S_VELO)`. Stage 174's accepted add/subtract residual rounding is
  preserved on final owners, with the independent old residual oracle intact.
- Velocity values needed by RK compatibility, including boundary/scaffold
  slots, are integrated from native tendency in persistent RK buffers. They
  are computed **before** the independent reference RK update because the last
  substage aliases input and output. Final-owner integrated publication and
  native inverse-transform processing remain unchanged.
- Non-oracle execution no longer calls `velocity_trend_source` or
  `velocity_trend_grad`. Guards reject such calls during the production
  transaction. The original routines remain available for independent oracles.

Mass compatibility, the existing temperature chain, external physics,
adaptation, remapping and checkpoint interfaces retain their accepted behavior.

## Independent checks and proof of the cut

Oracle runs compare compiled native source after every grid level and native
gradient/RK values bit-for-bit. The previous independently replayed source
oracle remains in addition to the original reference. No tolerance was relaxed.
Native tendencies are never initialized from captured Domain source/gradient
answers, and full oracle evaluation cannot overwrite the native workspace.

`WAVETRISK_TEST_VELOCITY_CUT` is a compile-time test, not a production switch. In
a non-oracle run it fills the Domain velocity tendency with NaNs and verifies
that it remains poisoned. Successful finite native integration then exercises
publication and RK without a usable completed Domain tendency.

New profiling records:

```text
native velocity plan
native velocity source
native velocity gradient
native velocity: plan builds direct restrictions gradients physics = ...
native velocity: published RK-values Domain-source Domain-gradient = N M 0 0
```

`Domain velocity compatibility` must have zero calls/ be absent in non-oracle
profiles. The last two counters must be zero; they are nonzero for deliberate
oracle evaluation. `native velocity residual: native values formed
oracle-checked` replaces the old `raw values` label. Publication message sizes
are unchanged. Timers are inclusive and must not be added as disjoint budgets.

## Local validation artifacts

All runs use four MPI ranks, isolated inputs, and the original checkpoint 3:
`/private/tmp/wavetrisk-stage175-production.Xxvrzf`.

Final builds use separate `build/stage175-final-*` and `bin/stage175-final-*`
directories. `src/shared.f90` remains RK4; RK3 uses a separate source overlay.
The user's standard run directory and `bin/climateJ5` are untouched.

The standalone tests run with `sh test/parallel_block_velocity/run.sh` and cover
all stencil inputs, operation order, mixed masks, inactive/poisoned operands,
gradient orientation and physical-layer pressure arithmetic at checked `-O0`
and `-O2`.

### Accepted local correctness results

- Fresh strict checked RK3 and RK4 builds: PASS, normal 132-column limit,
  warnings as errors, bounds checks and invalid/zero/overflow traps.
- Standalone numerical/program tests at `-O0` and `-O2`: PASS.
- Final RK3 oracle: PASS, one timestep / three substages, including exact
  compiled source, gradient, velocity RK and retained residual comparisons.
- Final RK4 oracle: PASS, one timestep / four substages, with the same checks.
- Checked non-oracle velocity-poison test: PASS, four substages, finite final
  state while the Domain velocity tendency remains NaN. Profile:
  `published RK-values Domain-source Domain-gradient = 15747840 39378240 0 0`.
- Full optimized RK4: PASS, ten timesteps, vertical remap, checkpoint 4,
  restart and two post-restart timesteps. Checkpoint 4 is byte-identical to
  the accepted Stage 174 checkpoint, SHA-256:
  `6086e1fab875369095b1dea9800af7a39261ab751d6f10916718ed3c081179f3`.
- First-eight-step native publication count: 129000960; native RK values:
  323910720; production Domain source/gradient calls: **0 / 0**. Post-restart:
  25781760 published, 63892800 RK values, **0 / 0** Domain calls. The Domain
  velocity compatibility phase is absent in both production profile windows.

Logs under the isolated directory:

```text
rk3-oracle/RK3_final_debug_oracle.log
rk4-oracle/RK4_final_debug_oracle.log
poison/RK4_poison.log
production/RK4_optimized_restart_profile.log
baseline/RK4_stage174_baseline_profile.log
```

### Matching local performance comparison

The accepted Stage 174 optimized binary and Stage 175 ran sequentially with
the same four-rank input, checkpoint and profiling flags. Reported total time
was **355.43 s -> 344.01 s**, approximately **3.2% lower**. All 22 printed
mass/state records agree after removing CPU timings, and the two checkpoint 4
files are byte-identical.

For the first eight steps, rank-average velocity kernel time fell from
**5.1806 s** (Domain velocity compatibility) to **3.5915 s** (native plan,
source and gradient combined), approximately **30.7% lower**, including native
plan setup. Publication and inverse-transform communication counts are
unchanged in this matched comparison.

This is one noisy, profiled local pair, not a repeated unprofiled benchmark or
an 83-rank speedup claim. The remaining shared primitive/mass pass and other
compatibility work still limit total savings. Cluster correctness tests and
matched production timing/profile runs are the next acceptance gate.

Tested final binaries (SHA-256):

```text
49e64e145bb0689c72ad67e026f251d4028a8d3125db6c60127406e41b1570ec  bin/stage175-final-rk3-check/climate
2c2dc059b0705adb9511c4e9770f7ec655bda3f868fde75f6dcc838078493975  bin/stage175-final-rk4-check/climate
462213e349f8e58ad494e7f13b131ed5ff058af29daa81dc18d2f6308179c874  bin/stage175-final-rk4-opt/climate
39db6be8a5350cbe1eef2ef424422a5d4cc18937f7c60f40a4f918951c78d264  bin/stage175-final-poison/climate
```

The poison build adds `-DWAVETRISK_TEST_VELOCITY_CUT` to the normal checked
preprocessor flags and runs with both oracles disabled. It is not the binary
to use for performance measurements. The optimized build retains an existing
GNU warning about the `air_temperature` allocation descriptor in unchanged
`parallel_block.f90`; the new code builds without warnings.

## Transfer and cluster builds (after local acceptance)

Five implementation files must be transferred together:

```sh
scp /Users/kevlahan/wavetrisk_hydrostatic/Makefile bbserv:~/wav/
scp \
  /Users/kevlahan/wavetrisk_hydrostatic/src/parallel_block_velocity.f90 \
  /Users/kevlahan/wavetrisk_hydrostatic/src/multi_level.f90 \
  /Users/kevlahan/wavetrisk_hydrostatic/src/parallel_block_mpi.f90 \
  /Users/kevlahan/wavetrisk_hydrostatic/src/time_integr.f90 \
  bbserv:~/wav/src/
scp /Users/kevlahan/wavetrisk_hydrostatic/docs/parallel-blocks-stage175.sha256 bbserv:~/wav/docs/
```

Then run `sha256sum -c docs/parallel-blocks-stage175.sha256` from `~/wav` on
the cluster. The manifest also verifies the unchanged accepted block, inverse
and adaptation dependencies. `shared.f90` is intentionally not in the manifest
because RK3/RK4 selection differs between builds.

Build on CPUs compatible with every execution node: optimized builds use
`-march=native`. Select RK3 in `src/shared.f90`, then:

```sh
make -j1 DEBUG=check PARAM=param_J5 TEST_CASE=climate \
  BUILD_DIR=build/stage175-rk3-check BIN_DIR=bin/stage175-rk3-check
```

Restore RK4 before:

```sh
make -j1 DEBUG=check PARAM=param_J5 TEST_CASE=climate \
  BUILD_DIR=build/stage175-rk4-check BIN_DIR=bin/stage175-rk4-check
make -j1 DEBUG=false PARAM=param_J5 TEST_CASE=climate \
  BUILD_DIR=build/stage175-rk4-opt BIN_DIR=bin/stage175-rk4-opt
```

These builds do not update `bin/climateJ5`. Copy the specified binary for each
run. Keep the same input/checkpoint and 83-rank allocation.

### RK3 debug oracle

```sh
cp ~/wav/bin/stage175-rk3-check/climate ./climateJ5
export WAVETRISK_VALIDATE_BLOCK_DYNAMICS=1
export WAVETRISK_VALIDATE_BLOCK_ADAPTATION=1
export WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS=0
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=0
set -o pipefail
srun -n 83 ./climateJ5 simple.in 2>&1 | tee RK3_debug_oracle.log
```

### RK4 debug oracle

```sh
cp ~/wav/bin/stage175-rk4-check/climate ./climateJ5
export WAVETRISK_VALIDATE_BLOCK_DYNAMICS=1
export WAVETRISK_VALIDATE_BLOCK_ADAPTATION=1
export WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS=0
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=0
set -o pipefail
srun -n 83 ./climateJ5 simple.in 2>&1 | tee RK4_debug_oracle.log
```

### RK4 production timing

```sh
cp ~/wav/bin/stage175-rk4-opt/climate ./climateJ5
export WAVETRISK_VALIDATE_BLOCK_DYNAMICS=0
export WAVETRISK_VALIDATE_BLOCK_ADAPTATION=0
export WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS=0
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=0
set -o pipefail
srun -n 83 ./climateJ5 simple.in 2>&1 | tee RK4_production.log
```

### RK4 production profiling

```sh
cp ~/wav/bin/stage175-rk4-opt/climate ./climateJ5
export WAVETRISK_VALIDATE_BLOCK_DYNAMICS=0
export WAVETRISK_VALIDATE_BLOCK_ADAPTATION=0
export WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS=0
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=1
set -o pipefail
srun -n 83 ./climateJ5 simple.in 2>&1 | tee RK4_profile.log
```

Assess speed using paired unprofiled Stage 174/175 runs on the same CPUs and
allocation. Shared primitive work, owner publication and native residual
recomposition still cost time; this stage makes no predetermined speedup claim.
