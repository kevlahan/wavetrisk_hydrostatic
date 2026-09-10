# Stage 177: native mass restriction, boundary and RK chain

Base: accepted Stage 176, commit `ec522cab`. Stage 177 is accepted for commit.

## Accepted cluster results (2026-09-10)

The supplied 83-task RK3 and RK4 debug oracle runs both completed, including
restart, with zero cache/final mismatches. Optimized production and profiling
runs also completed and agree with RK4 oracle printed numerical state.
Production time was 21.284 s (Stage 176: 21.839 s); profiling time was
21.381 s (Stage 176: 22.054 s). These single runs suggest only a modest
2.5–3.1% improvement, not a controlled repeatability claim.

The first-eight-step native inverse-transform time remains 3.8343 s
rank-average, including 2.0146 s boundary synchronization. Repeated native
gather/alias/scatter and inverse kernel work are the next target. Inclusive
timers overlap and must not be added as independent costs.

## Production cut

The once-evaluated shared primitive/physics producer writes mass flux directly
to native numerical storage on geometry owners. Ordered, precompiled mass
restriction and divergence stencils operate on that storage. Per-level MPI
plans complete signed edge and scalar-node boundaries without Domain field
writeback, boundary callbacks, or re-import. Velocity consumes the native flux
at its original fine-to-coarse execution phase.

The native mass RHS supplies mass RK closure and final-owner publication.
Production **does not also replay mass direct flux, restriction or divergence
on final owners**. It reuses the existing scalar wire channel: mass carries its
completed native RHS in the first of three dynamic slots (the others are zero),
while temperature still carries three physics inputs. No additional runtime
collective or routing-setup collective is introduced, apart from an optional
profiling-counter reduction. The wire geometry and
temperature contract remain unchanged.

This removes the legacy production mass restriction/divergence and mass RK
tendency reads together. It is not a completed removal of every Domain field:
shared PV/pressure/KE, physics/diffusion prepasses, solution compatibility at
primitive consumers, adaptation and external interfaces remain. Native mass
and velocity execution still resides on geometry owners, as previously agreed.
The existing scalar message layout retains unused mass slots; this stage targets
computation rather than another traffic-volume reduction.

## Numerical and lifecycle contracts

- Restriction preserves child/traversal order and partial/coarse/small-flux
  arithmetic. Immutable interpolation/overlap coefficients and addresses are
  compiled once per generation, with exact active/restriction mask checks.
- Native boundary plans reproduce remote level filtering, edge signs and
  ordered same-rank copies. In particular, legacy same-rank copies touch all
  levels during a single-level exchange; that behavior is not silently narrowed.
- Flux and RHS are consumed at their original level phases, not reconstructed
  from an end-of-pass snapshot after storage aliases have changed.
- Native mass RK values are formed before the last RK substage aliases its
  input/output arrays. Mass RHS readiness is checked before publication and RK.
- Generation changes rebuild workspaces and boundary plans after adaptation,
  ownership changes and restart. Mask-only changes rebuild affected stencils.

## Oracle and poison checks

Oracle mode retains the separate Domain restriction/divergence plus the full
independent `trend_ml` oracle. The new native phase checks compare flux inputs,
divergence and mass RK values **bit-for-bit**; existing final-owner oracle
tolerances are unchanged. The shared direct primitive/physics computation is
evaluated once and supplies both phase paths; the complete legacy tendency
oracle remains independent.

`WAVETRISK_TEST_MASS_CUT` is a compile-time local test, not a production runtime
switch. It poisons Domain mass tendency and the advective mass-flux scratch
after the retained diffusion prepass, checks that neither is written by the
native chain, and restores only the diffusion scratch for the next layer.
The checked non-oracle run must complete without using those poisoned values.

New profiling rows:

```text
native mass plan
native mass restriction/boundary
native mass: plans restricted-edges divergence boundary published RK-values = ...
```

Production `Domain mass-flux compatibility` must have zero calls (the row may
be omitted); existing Domain velocity source/gradient counters remain zero.
Plan and kernel counters must be nonzero where there is corresponding work.
Inclusive profile times must not be summed as disjoint costs.

## Build and cluster tests

Transfer **Makefile plus all five source files**, not only the old three-file set:

```sh
scp Makefile bbserv:~/wav/
scp src/parallel_block_mass.f90 src/multi_level.f90 src/ops.f90 \
  src/parallel_block_mpi.f90 src/time_integr.f90 bbserv:~/wav/src/
```

Use fresh build directories and serial builds. Compile on a CPU compatible
with the execution partition: the optimized build uses `-march=native`.
Select RK3 in `src/shared.f90` for the first build, then restore RK4 for the
remaining builds. Do not overwrite an executable until the build succeeds.

```sh
# With timeint_type = "RK3":
make -j1 DEBUG=check PARAM=param_J5 TEST_CASE=climate \
  BUILD_DIR=build/stage177-rk3-check BIN_DIR=bin/stage177-rk3-check
# Restore timeint_type = "RK4":
make -j1 DEBUG=check PARAM=param_J5 TEST_CASE=climate \
  BUILD_DIR=build/stage177-rk4-check BIN_DIR=bin/stage177-rk4-check
make -j1 DEBUG=false PARAM=param_J5 TEST_CASE=climate \
  BUILD_DIR=build/stage177-rk4-opt BIN_DIR=bin/stage177-rk4-opt
```

From the established cluster test directory, retaining the checkpoint-3
restart input for every run:

### RK3 debug oracle

```sh
cp ~/wav/bin/stage177-rk3-check/climate ./climateJ5
export WAVETRISK_VALIDATE_BLOCK_DYNAMICS=1
export WAVETRISK_VALIDATE_BLOCK_ADAPTATION=1
export WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS=0
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=0
set -o pipefail
srun -n 83 ./climateJ5 simple.in 2>&1 | tee RK3_debug_oracle.log
```

### RK4 debug oracle

```sh
cp ~/wav/bin/stage177-rk4-check/climate ./climateJ5
export WAVETRISK_VALIDATE_BLOCK_DYNAMICS=1
export WAVETRISK_VALIDATE_BLOCK_ADAPTATION=1
export WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS=0
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=0
set -o pipefail
srun -n 83 ./climateJ5 simple.in 2>&1 | tee RK4_debug_oracle.log
```

### RK4 optimized production timing

```sh
cp ~/wav/bin/stage177-rk4-opt/climate ./climateJ5
export WAVETRISK_VALIDATE_BLOCK_DYNAMICS=0
export WAVETRISK_VALIDATE_BLOCK_ADAPTATION=0
export WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS=0
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=0
set -o pipefail
srun -n 83 ./climateJ5 simple.in 2>&1 | tee RK4_production.log
```

### RK4 optimized profiling

```sh
cp ~/wav/bin/stage177-rk4-opt/climate ./climateJ5
export WAVETRISK_VALIDATE_BLOCK_DYNAMICS=0
export WAVETRISK_VALIDATE_BLOCK_ADAPTATION=0
export WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS=0
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=1
set -o pipefail
srun -n 83 ./climateJ5 simple.in 2>&1 | tee RK4_profile.log
```

Both oracle tests, checked non-oracle tests, and the optimized restart must
pass before accepting this stage. Use paired same-allocation runs against
Stage 176 for speed claims, not absolute times from different cluster loads.

## Local validation

Artifacts: `/private/tmp/wavetrisk-stage177.xWaLB8`. The standard user test
directory and `bin/climateJ5` are untouched. No cluster run is claimed.

- Standalone mass kernel tests at checked `-O0` and `-O2`: PASS. These cover
  all nonempty edge masks, asymmetric impulses, ordered overlapping writes,
  and inactive divergence with poisoned geometry. Existing velocity tests:
  PASS at both optimization levels.
- RK4 debug oracle: PASS, one timestep/all four substages, including bitwise
  native mass phase/RK checks and the existing complete tendency oracle.
- Checked non-oracle RK4 with `WAVETRISK_TEST_MASS_CUT`: PASS, all four substages.
  No Domain advective mass flux/tendency use or write was detected. State output
  matches the accepted baseline. Published native mass RHS count is 5249280;
  native mass RK-value count is 13126080.
- Full optimized RK4 restart: PASS, ten timesteps including vertical remapping,
  checkpoint 4, restart and two post-restart steps. All 22 printed numerical
  state/mass records match Stage 176. Checkpoint 4 is byte-identical:
  `6086e1fab875369095b1dea9800af7a39261ab751d6f10916718ed3c081179f3`.
- Production has no Domain mass-flux compatibility calls and zero Domain
  velocity source/gradient calls before and after restart. Native mass RHS
  publication counts are 43000320 and 8593920 in the two profile windows.
- Optimized RK4 oracle cross-check: PASS, all four substages, including the
  exact native mass comparisons under `-O2 -march=native`.
- Final RK3 debug oracle: PASS, all three substages, exact mass phase/RK
  checks, zero cache/final mismatches, and baseline-identical printed state.

The optimized development run reported 241.77 s, versus the earlier Stage 176
run's 380.39 s. This is encouraging, **not a controlled speedup measurement**:
the runs were not paired, and builds overlapped part of the Stage 177 run.
First-eight-step native mass plan compilation is 0.071752 s rank-average;
native mass restriction/boundaries/divergence total 5.1703 s. The old Domain
mass-flux compatibility row disappears. These are inclusive measurements.
Repeat same-allocation cluster timings to establish the actual improvement.

Logs:

```text
rk4-oracle/RK4_second_oracle.log
checked-production/RK4_mass_cut.log
production/RK4_optimized_restart_profile.log
rk3-oracle/RK3_debug_oracle.log
optimized-oracle/RK4_optimized_oracle.log
```

The initial local build/run caught and corrected a missing import and an
index-name collision in the plan compiler. An optimized impure-function
elimination warning introduced by this stage was removed. Final rebuilds have
no new warnings. The fresh optimized build also reported the pre-existing
`air_temperature` maybe-uninitialized warning in `parallel_block.f90`.

Executed binaries (SHA-256), preserved alongside their logs:

```text
40257f6549c4f7d4e92dbb1a14d5fc90d5fb026be9cd8fe6b7551feb7bda795d  rk4-oracle/climate
40257f6549c4f7d4e92dbb1a14d5fc90d5fb026be9cd8fe6b7551feb7bda795d  checked-production/climate
cf5e5d55ba1e46b59dacdd451e4996930726a45db3ba2fa844daa5104db67de6  production/climate
fe711a954f8dac53324188a6cf241333129db6b82b01eb0322dc1085cb941abe  rk3-oracle/climate
1f7146df97c7f05aea3c9980ba00cd79886f8fe3f3a35143cbc66bd8e1c0df66  optimized-oracle/climate
```

The first three runs precede final profiling-accounting/comment cleanup and
an MPI profiling error-check placement correction; the RK3 and optimized
oracle tests use the final source. No numerical kernel changed after the
successful checked RK4/poison runs. Oracle runs performed concurrently are
correctness tests, not timing comparisons.
