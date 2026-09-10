# Restart correctness corrections — 179d

## Scope

Correctness gate before further profiling, smaller fixtures, or optimizations.
The authoritative `main` revision remains
`b82647b36eb825aed0c0aae4e92adcf4ab46d608`. All legacy instrumentation and builds
are isolated copies. The original checkpoint-3 J7/30-atmospheric-layer fixture
and adaptation tolerance are unchanged. Diagnostic simulations use four local
MPI ranks sequentially and are **not timing benchmarks**.

The local correctness gates below passed. This checkpoint is ready to commit
before any smaller climate configuration or further optimization is attempted.
It does not claim validation of the corrected code on the 83-rank cluster yet.

## Corrections

1. **Scalar wavelet coverage at block roots.** Iterate eligible target patches,
   rather than only children of retained parents. Roots whose parent is outside
   the compact block must participate in the final forward transform. The
   provisional transform still excludes them. The independent oracle traverses
   the Domain target set, checks expected coverage, and evaluates the legacy
   scalar stencil independently of native producer routes.
2. **Forward and inverse transform ranges are distinct.** The final forward
   transform includes the coarse parent at `level_start-1`; the legacy inverse
   starts at `level_start`. Using the forward range for both reconstructed an
   additional level in native execution.
3. **Physics boundary input epoch.** Native inverse kernels require refreshed
   aliases after reconstruction phases. Legacy physics instead observes scalar
   aliases from the post-prolongation phase and finest-level velocity aliases
   from the post-outer/pre-inner phase. Retain those values in compact alias-only
   storage and publish them after the fresh-state mirror checks, immediately
   before the Domain physics consumer. Never re-import this older consumer image
   as a native kernel dependency. This preserves legacy behavior; changing that
   behavior would require a separate numerical change, not a performance fix.
4. **Production polar remap coverage.** Native remapping visits compact patch
   interiors only. Legacy `apply_no_bdry2` also visits polar columns. Previously
   only enabling `WAVETRISK_VALIDATE_BLOCK_REMAP` supplied that missing work.
   Production now remaps those polar scalar columns explicitly, before native
   interior writeback, using a column-local mass snapshot. It neither repeats
   the interior Domain remap nor copies the full Domain mass field. Velocities
   at the pole remain untouched, as in legacy.
5. **Incremental build dependencies.** Declare the affected Fortran module
   dependencies so a private derived-type change cannot leave stale transitive
   module consumers in a warm build.

## Evidence collected during localization

- Initial checkpoint-3 interiors and masks were exact between legacy and block.
- Missing scalar roots accounted for 230,400 nonzero coefficients per scalar.
- After corrections 1–3, sampled first-step dynamics, physics, adaptation, and
  the step ending with vertical remap agree with legacy to floating-point scale;
  sampled soil fields are exact.
- The first remaining large discrepancy appears in the **first dynamics step
  after remap** (cumulative step 257), not at checkpoint serialization. The
  earlier interior-only snapshots did not include polar storage.
- The focused polar callback regression extracts both actual routines and
  compares every resulting field exactly at `-O0` and `-O2`, with bounds checks.
  It covers both polar addresses, active/inactive masks, and unchanged velocity.
- The original-resolution, first-step RK4 `DEBUG=check` dynamics/adaptation
  oracle passed all four sub-stages after corrections 1–3.
- The helper suite has 22 passing tests, including root-omission
  negative control, snapshot/checkpoint parser checks, and polar remap coverage.

## Production checkpoint validation

The corrected optimized RK4 run completed checkpoint 4, its immediate reload,
and the following time step with exit status 0. Its original input remains
J7, `tol=2.5e-2`, 30 atmospheric layers, and 10 soil layers plus surface.
All three oracles were explicitly disabled.

Checkpoint headers, thresholds, required-patch topology, and all soil values
are exact. No nonfinite values or missing patch keys occurred.

| Atmospheric quantity | Maximum absolute difference |
| --- | ---: |
| Checkpoint coarse velocity | 1.027e-12 |
| Checkpoint coarse mass | 1.592e-12 |
| Checkpoint coarse mass-weighted temperature | 4.366e-10 |
| Checkpoint velocity wavelets | 6.687e-12 |
| Checkpoint mass wavelets | 3.402e-12 |
| Checkpoint temperature wavelets | 1.784e-9 |
| Checkpoint polar mass | 1.137e-13 |
| Checkpoint polar mass-weighted temperature | 2.911e-11 |
| Reloaded active velocity | 8.296e-12 |
| Reloaded active mass | 3.752e-12 |
| Reloaded active mass-weighted temperature | 1.459e-9 |

This is numerical agreement, **not bitwise equality**. Immediately before
checkpointing, 38 active velocity solution samples differ by more than 1e-9,
with maximum absolute difference 9.542e-7. These differences are approximately
one single-precision ULP (maximum 1.00051 ULP units), consistent with the
single-precision physics conversion acting on slightly different double
precision inputs. The velocity wavelets and serialized/reloaded state agree
at double-precision floating-point scale. Soil remains exact. No oracle
tolerance or numerical threshold was changed to accept these results.

The first post-remap step also agrees after the polar correction: active
velocity/mass/temperature maxima are respectively 8.271e-12, 3.297e-12,
and 1.019e-9, with exact soil fields. The previous O(1e-3–1e2) errors are gone.

Local evidence: `/private/tmp/wavetrisk-restart-probe-runs-09` compared against
the numerical legacy reference in `/private/tmp/wavetrisk-restart-probe-runs-06/legacy`.
The former contains compact JSON reports as well as raw diagnostic snapshots.

## Checked validation

- RK4, original-resolution first step, clean `DEBUG=check` build: all four
  dynamics/adaptation oracle sub-stages passed; exit 0. This run preceded the
  polar correction and did not reach a remap event.
  Binary SHA256: `081bee8e08e92f4740f157bfd05e86e4370199321bed3f6b4f548cc83d55d60a`.
- RK3, original-resolution first step, clean `DEBUG=check` build including the
  polar correction: all three dynamics/adaptation oracle sub-stages passed;
  exit 0. The remap switch was explicitly enabled, but the first step only
  exercised its minimum-mass check, not remapping itself.
  Binary SHA256: `1a7324c511e30840ea92d27934bb31f99ee0996cc8f1df20f727f285c0fc18c9`.
- Profiling instrumentation and inverse-route Fortran regressions passed at
  `-O0` and `-O2`; the 132-column source check passed. The focused polar callback
  test also enables bounds checks, signalling-NaN initialization and FPE traps.
- RK4 `DEBUG=check`, five original-resolution steps through the actual remap
  event: minimum-mass checks and the native/legacy remap comparison passed;
  exit 0. Dynamics/adaptation oracles were off and remap oracle was on, so this
  independently exercises production dynamics leading into the remap check.
  Binary SHA256: `4624b5ca60460517b10855eb8ceb34106743d292cbfb8fc8de5d70408dc40a17`.
  Run: `/private/tmp/wavetrisk-restart-probe-rk4-remap-check-01/block/run.log`.

## Checkpoint and next gate

Production restart comparison, checked RK3/RK4 sub-stages, checked remap,
regression tests, and source-width checks have passed. Binary identities and
numerical differences are recorded above and in
`restart-mismatch-179d/corrected-rk4-summary.json`. Source transfer checksums are
in `parallel-blocks-stage179d.sha256`.

Commit this corrected/profiling checkpoint before proceeding. The next work may
use reduced-memory fixtures for profiling, but should retain this original-case
correctness gate and repeat the cluster RK3/RK4 oracle and production restart
tests when the 83-rank allocation is available. No new speedup is claimed here.

## Oracle settings for isolated validation

`run_restart_probe.py --oracle` now explicitly sets all three independent
oracles to 1:

```sh
export WAVETRISK_VALIDATE_BLOCK_DYNAMICS=1
export WAVETRISK_VALIDATE_BLOCK_ADAPTATION=1
export WAVETRISK_VALIDATE_BLOCK_REMAP=1
export WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS=0
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=0
export WAVETRISK_PROFILE_BLOCK_DETAIL=0
```

The production comparison explicitly sets all three oracles and all profiling
switches to 0. `--remap-oracle` enables only the remap oracle, independently of
dynamics and adaptation. Before this correction the diagnostic helper's
`--oracle` switch enabled only dynamics/adaptation; earlier one-step results
must not be described as remap-oracle validation.
