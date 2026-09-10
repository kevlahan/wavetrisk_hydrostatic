# Stage 180: correctness-gated local profiling

This stage changes profiling/test helpers, not the numerical algorithm. It builds
the corrected block checkpoint `1b01510ae655eee0e2aab1a8d911e55c147173f2` and compares
against unchanged authoritative main `b82647b36eb825aed0c0aae4e92adcf4ab46d608`.
Experimental legacy instrumentation is confined to a separate archive.

## What changes

- Pin the common-boundary instrumentation to committed revisions. Preparing an
  instrumented legacy archive no longer depends on an uncommitted Git diff.
- Generate a genuinely evolved J6 checkpoint using unchanged legacy's supported
  refinement cap. Keep 30 atmospheric levels, 10 soil levels, tolerance 0.025,
  and the original J7 checkpoint-3 fixture unchanged. Do not rename checkpoint 3
  or edit its binary header to manufacture checkpoint 4.
- Hash seed input parameters, checkpoint contents, and all grid files, following
  nested directory links with cycle detection. Bind validation to binaries,
  rank count and single-thread settings; record clean build/source identities.
- Require complete semantic checkpoint comparison over at least two newly
  written checkpoints, actual reloads and a remap before launching timing pairs.
  Compare unchanged legacy with block production, block with all three oracles,
  experimental legacy with instrumentation disabled/enabled, and block detail.
- Retire the old printed-state-only local timing entry point. Printed agreement
  remains a secondary timing-run check, never the full-field acceptance gate.
- Separate ordinary, remap and checkpoint steps. Profile attribution uses each
  timestep's actual slowest rank, not the sum of unrelated region maxima.
- Sample simulation-child RSS and system VM counters. Screen out pressure-affected
  timing pairs; these counters are not application memory bandwidth measurements.

## Numerical acceptance and limitations

The fixed local screen compares every supported checkpoint coarse/wavelet field,
polar mass/temperature, required-child flags, headers and thresholds. Atmospheric
absolute limits are velocity/mass `1e-10`, temperature `1e-7`; soil, thresholds and
topology must agree exactly. Nonfinite values fail. Limits are declared before
the experiments, not fitted to candidate errors. Experimental legacy requires
zero numerical field difference from unchanged legacy.

This is **not byte equality**, a universal error bound, or a comparison of every
transient work array at every RK substage. Two complete checkpoint/reload cycles
exercise persistent-state evolution; the existing independent RK/adaptation/remap
oracles remain complementary. The local campaign uses the optimized RK4 build
with runtime oracles on/off; it does not replace checked RK3/RK4 cluster tests.

Four unpinned local ranks cannot predict 83-rank cluster scaling. Reducing maximum
refinement removes level-7 work, so retain the original J7 case for later cluster
confirmation. No numerical optimization or shared-geometry layout rewrite is
included in this stage.

## Commands

Use new, external output directories. The existing verified legacy build may be
reused; its tracked source is verified against main before and after execution.

```sh
python3 test/parallel_block_profile/build_block_profile.py \
  --repo "$PWD" --out /private/tmp/wavetrisk-stage180-block-opt-01
python3 test/parallel_block_profile/prepare_legacy_profile.py \
  --repo "$PWD" --baseline /path/to/verified-legacy-main \
  --out /private/tmp/wavetrisk-stage180-legacy-profile-01
make -C /private/tmp/wavetrisk-stage180-legacy-profile-01 \
  -j1 DEBUG=false PARAM=param_J5 TEST_CASE=climate
python3 test/parallel_block_profile/prepare_light_fixture.py \
  --repo "$PWD" --legacy-source /path/to/verified-legacy-main \
  --fixture /path/to/original-checkpoint3-fixture \
  --out /private/tmp/wavetrisk-stage180-j6-fixture-01
python3 test/parallel_block_profile/local_campaign.py \
  --repo "$PWD" --legacy-source /path/to/verified-legacy-main \
  --instrumented-source /private/tmp/wavetrisk-stage180-legacy-profile-01 \
  --block-source /private/tmp/wavetrisk-stage180-block-opt-01 \
  --fixture /private/tmp/wavetrisk-stage180-j6-fixture-01/seed \
  --out /private/tmp/wavetrisk-stage180-j6-campaign-01
python3 test/parallel_block_profile/summarize_local_campaign.py \
  /private/tmp/wavetrisk-stage180-j6-campaign-01/results.json
```

The macOS runner launches simulations sequentially via `mpirun -n 4`, from fresh
run directories. It explicitly clears inherited `WAVETRISK_*` switches, disables
all three oracles for production/timing, enables all three for oracle validation,
and enables both profiling switches only in detail runs. No authoritative main
file is modified. A failed gate stops before timing; do not widen the limits to
make an unexplained mismatch pass.

## Local test outcome (2026-09-10)

Both clean builds succeeded. The helper suite passes 37 tests, and the checked
Fortran profiler regression passes (nesting, exclusive accounting, recursion,
boundary attribution, reset and overflow). No numerical Fortran source changed.
The rebuilt experimental legacy archive also passes two full-interval
transparency checks: instrumentation disabled and detail instrumentation enabled
both preserve all compared checkpoint fields exactly against unchanged legacy.
The enabled run's four rank-local profile files pass exclusive wall/CPU
conservation checks. These separate diagnostics do not override the failed
block-versus-legacy numerical gate.

The derived checkpoint contains 640 level-5 patches and 1,270 level-6 patches,
with nonzero evolved wavelets and no level-7 patches. The original fixture and
authoritative main remain unchanged.

The first campaign deliberately **did not launch timing pairs**. Unchanged
legacy and production block completed the longer validation interval, but the
checkpoint-6 comparison exceeded the fixed numerical screen:

| Maximum atmospheric wavelet difference | Checkpoint 5 | Checkpoint 6 |
| --- | ---: | ---: |
| Velocity | 9.305e-12 | 5.997e-9 |
| Mass | 2.468e-12 | 5.704e-9 |
| Mass-weighted temperature (not temperature in kelvin) | 2.409e-9 | 2.788e-6 |

Headers, thresholds, required-child topology and soil values agree exactly;
there are no nonfinite samples. The largest checkpoint-6 differences in these
fields occur in Domain 77 at vertical level 26. Checkpoint 5 passes the screen.
This is not yet a diagnosed algorithm defect: the existing single-precision
physics conversion is a possible amplifier, but needs localization rather than
an automatic tolerance increase. The new fixture extends the independent
comparison beyond the earlier original-case checkpoint-4 validation.

There is also a separate timing-quality failure. Sampled simultaneous child RSS
peaks are approximately 3.18 GiB (legacy) and 8.53 GiB (block), excluding compressed
pages. During block execution, system-wide decompression totals about 468 GiB,
versus about 22 MiB during legacy execution. This is repeated VM activity, **not
468 GiB of allocated application memory or a measurement of application DRAM
traffic**. No swapout was observed; compression alone can still confound timing.
The J6 cap therefore does not establish a pressure-free local benchmark.

Evidence directories:

- `/private/tmp/wavetrisk-stage180-j6-fixture-01`: generation log, immutable seed
  and identity/topology manifest.
- `/private/tmp/wavetrisk-stage180-j6-campaign-01`: raw legacy/production logs,
  checkpoints 5/6, per-run identities, memory samples and failed gate report.
- `parallel-blocks-local-protocol-180-results.json` in this documentation
  directory: compact, reproducible summary including oracle and legacy
  instrumentation comparisons. Raw runs remain outside the repository.

The separate optimized RK4 oracle diagnostic completed with all three runtime
oracles enabled and no assertion failure. Its checkpoints pass the numerical
screen against oracle-disabled block execution. The largest on/off field
difference is 1.756e-9 (mass-weighted temperature); the maximum differences from
legacy listed above are unchanged. Thus there is no evidence here of an oracle
repairing the larger discrepancy. This does not identify its cause.

The next numerical diagnostic is localization of the first enlarged difference
between checkpoints 5 and 6, particularly the physics precision conversion and
remap/restart transitions around the affected column.
Do not infer a new speedup, waive the screen, or start a numerical optimization
from these validation-run timings. The paired timing campaign remains
unaccepted until both numerical and memory quality gates pass.
