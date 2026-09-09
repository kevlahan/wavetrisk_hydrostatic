# Stage 168: single-evaluation scalar physics capture

Base: `5e1ffecc` (Stage 167, native scalar restriction authoritative).

## Change and scope

Stage 167 moved native scalar restriction ahead of physical recomposition.
Production nevertheless rebuilt the interior Domain reference records on every
RK stage. This evaluated `physics_scalar_flux` a second time and repeatedly
constructed geometry and six-edge stencils that the compact transport discarded.

Stage 168 removes this duplicated capture work:

- `step1` exports its already computed positive-edge physics residuals through
  one patch-sized scratch array. Production no longer calls the physics function
  a second time to obtain those residuals.
- A map built with the writeback plan gives each producing Domain patch its
  retained-block or outgoing-stream address. It requires no extra collective.
- After geometry has been initialized for a plan generation, interior capture
  touches only the three physics residual values per scalar/vertical-level/node.
  The interior direct/restricted-flux reference walks are skipped.
- Remote compact residuals are written directly into an existing transport
  workspace. The separate gather from 50-value records is removed.
- Initial production interior transport also drops unused oracle field values:
  for `F` stored scalar/vertical-level fields per node, its payload falls from
  `33 + 17*F` doubles to `33 + 3*F`. The 33 shared geometry values and explicit
  received-block manifest are retained. Oracle payloads retain their old layout.
- Compact boundary capture refreshes the three live flux values without
  reconstructing boundary geometry.

Native arithmetic and the restriction sweep retain Stage 167's ordering. The
initial geometry capture is still required after adaptation or repartitioning.
Domain scalar restriction/trend and velocity compatibility still execute:
`Qperp` consumes restricted mass flux, and replacing the boundary flux producer
requires signed **edge** ownership, which differs from the validated scalar-node
ownership used to route dscalar. Reusing that node route for flux would be wrong.

This stage targets capture CPU time and memory traffic. It does not establish
that the branch is within 20% of legacy performance. Assess end-to-end speed
using repeated, unprofiled runs on the same 83-task allocation.

## Validation

With the dynamics oracle enabled, the original independent capture traversal
recomputes direct flux, edge metrics, and physics residuals and compares them
with the fused capture before the existing native tendency checks. The new map
checks duplicate assignments during construction. Every run also checks that
every mapped physical patch/vertical level was refreshed exactly once per RK
stage, including when its geometry was cached.

Local runs use four MPI tasks and the existing climate checkpoint 3, in
isolated temporary output directories. The repository's RK4 selection is
unchanged; a temporary copy of `shared.f90` selects RK3 for its separate test
build.

- GNU Fortran `DEBUG=check` build passed with `-Werror`, bounds checks, and
  invalid/zero/overflow traps. The optimized build also passed.
- Short RK4 debug oracle restart passed all four stages (one timestep,
  `time_end=0.3010`). Log:
  `/private/tmp/wavetrisk-stage168-test.35Cia2/RK4_debug_oracle.log`.
- Short RK3 debug oracle restart passed all three stages (one timestep,
  `time_end=0.3010`). Log:
  `/private/tmp/wavetrisk-stage168-rk3.qeqtYg/RK3_debug_oracle.log`.
- Final RK4 checked production restart passed all four stages with both oracles
  off, bounds checks and floating-point traps on (`time_end=0.3010`). This covers
  the reduced initial payload and subsequent compact refreshes. Log:
  `/private/tmp/wavetrisk-stage168-checked-production.oxWRd3/RK4_checked_production.log`.
- The optimized reference built from `5e1ffecc` completed ten timesteps through
  checkpoint 4 and restart (`time_end=0.3120`). Final Stage 168 completed the same
  ten timesteps successfully. All printed mass minima, simulated times,
  timesteps, Jmax values and degree-of-freedom counts match the reference.
- `cmp` confirms that the optimized Stage 168 and reference checkpoint 4 files
  are byte-for-byte identical, including the saved fields and adaptive topology.
  Both compressed checkpoint files have SHA-256:
  `6086e1fab875369095b1dea9800af7a39261ab751d6f10916718ed3c081179f3`.

The oracle test executables preceded the final production-only transport
reduction. Those changes are inactive with the oracle enabled, whose payload
format is unchanged. The production tests exercise the reduced initial stream
and compact direct-send path.

The first eight timesteps (the profiling interval ending at checkpoint 4) give:

| Measurement | Stage 167 | Stage 168 |
|---|---:|---:|
| Capture seconds, rank average | 25.023 | 18.597 |
| Capture seconds, maximum rank | 27.978 | 19.709 |
| Initial full-transfer bytes, global | 2,300,317,568 | 2,050,071,936 |
| Compact-transfer bytes, global | 518,937,984 | 518,937,984 |

These are local profiling observations, not a controlled speed benchmark.
Independent background processes were consuming substantial CPU during the
Stage 168 run, and compilation overlapped part of the reference run. No overall
speedup can be established from this pair. The byte-count reduction is exact;
the unchanged compact byte count is expected because direct packing removes
CPU/memory work without changing that message format. Boundary compatibility
records still account for much of the full-transfer traffic.
The reported timestep totals were 387.22 s (reference) and 418.84 s (Stage 168),
so these runs do not demonstrate an end-to-end improvement. Use the cluster
timing procedure below to measure that under comparable conditions.

Optimized logs:

- Reference: `/private/tmp/wavetrisk-stage168-baseline.InYaql/RK4_baseline_profile.log`
- Stage 168: `/private/tmp/wavetrisk-stage168-production.tc4RzD/RK4_stage168_profile.log`

Source SHA-256 values for transfer verification:

```text
430d12066405482d1df3e4d5b6139bc423d3fed5c7c234fb9dc8c1d767a1d8b4  src/ops.f90
d1ebfa94553fe9e1ddb198ff24fefb54d2c9ec41a8817f26944c15c5372896e7  src/multi_level.f90
ef7bebdbfde0adceea0cb4a0b7dbd68098101a3c5d11488b371965782f50d6c3  src/parallel_block_mpi.f90
```

## Files to transfer

```sh
scp \
  /Users/kevlahan/wavetrisk_hydrostatic/src/ops.f90 \
  /Users/kevlahan/wavetrisk_hydrostatic/src/multi_level.f90 \
  /Users/kevlahan/wavetrisk_hydrostatic/src/parallel_block_mpi.f90 \
  bbserv:~/wav/src/
```

`step1` has a new optional argument, so rebuild its callers as well. Use a fresh
build directory for each integrator/build mode. Select `timeint_type = "RK3"`
or `"RK4"` in `src/shared.f90` before the corresponding build, as in the existing
workflow; the climate input file does not select the integrator.

For example, with RK4 selected:

```sh
make -j1 DEBUG=check PARAM=param_J5 TEST_CASE=climate \
  BUILD_DIR=build/stage168-rk4-debug && mv bin/climate bin/climateJ5
```

Use `build/stage168-rk3-debug` with RK3 selected. For optimized RK4:

```sh
make -j1 DEBUG=false PARAM=param_J5 TEST_CASE=climate \
  BUILD_DIR=build/stage168-rk4-opt && mv bin/climate bin/climateJ5
```

The `&&` prevents renaming an old executable after a failed build. Copy the
new executable into the test directory after every build and check the printed
`timeint_type` before interpreting the results.

## Run settings

All four relevant environment variables are specified explicitly below. No new
stage-specific switch is needed.

RK3 and RK4 debug oracle runs:

```sh
export WAVETRISK_VALIDATE_BLOCK_DYNAMICS=1
export WAVETRISK_VALIDATE_BLOCK_ADAPTATION=1
export WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS=0
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=0
cp ~/wav/bin/climateJ5 .
set -o pipefail
srun -n 83 ./climateJ5 simple.in 2>&1 | tee RK4_debug_oracle.log
```

Use `RK3_debug_oracle.log` for the RK3 executable. Also run the RK4 debug binary
with dynamics/adaptation validation both set to `0` to exercise the compact
production path under bounds checks and floating-point traps.

Optimized RK4 timing:

```sh
export WAVETRISK_VALIDATE_BLOCK_DYNAMICS=0
export WAVETRISK_VALIDATE_BLOCK_ADAPTATION=0
export WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS=0
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=0
cp ~/wav/bin/climateJ5 .
set -o pipefail
srun -n 83 ./climateJ5 simple.in 2>&1 | tee RK4_optimized_timing.log
```

Optimized RK4 profiling uses the same settings except:

```sh
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=1
srun -n 83 ./climateJ5 simple.in 2>&1 | tee RK4_optimized_profile.log
```

Compare `restriction capture` and full/compact initial transport, then dynamics and
the complete timestep. Fused capture now occurs inside the basic-operator pass,
so the inclusive operator-compatibility timer includes it. Inclusive parent and
child timings must not be added together.

Accept this checkpoint after both debug oracles, checked production, and the
optimized restart test pass. The numerical trajectory should agree with Stage
167. Commit only after reviewing those results.
