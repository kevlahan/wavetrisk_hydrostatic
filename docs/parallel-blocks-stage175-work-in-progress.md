# Stage 175 implementation status — internal gate, NOT a performance release

Historical internal-gate record. The user subsequently approved native
geometry-owner execution first. Production replacement work and current
acceptance status are now documented in `parallel-blocks-stage175.md`.

Base: Stage 174, `d75e4634`. Changes below are uncommitted.

The complete delivery contract remains in `parallel-blocks-stage175-plan.md`.
This working tree does **not** meet that contract yet. Do not transfer these
changes as a completed Stage 175 or interpret passing component checks as a
production speedup. Stage 174 remains the accepted cluster checkpoint.

## Implemented

- `src/parallel_block_velocity.f90`: a numerical module depending only on
  `kind_mod`, with no Domain, MPI, callback, or mutable-workspace dependencies.
  It contains Qperp geometry weights, the ten-term Qperp stencil, direct
  physics/source combination, per-edge child-source restriction, complete
  Bernoulli/Exner gradient and compressible pressure/geopotential layer update.
  Operation ordering follows the original routines, including the negative
  diagonal gradient orientation and pressure midpoint calculation.
- `src/multi_level.f90`: an oracle-only adapter and independently populated
  source workspace. Fine sources are computed from primitive mass flux/PV,
  geometry and measured physics. Coarse sources use native child values,
  including absent-child and mixed-mask cases. No completed Domain source is
  copied into that workspace. Every scheduled level is compared bit-for-bit
  against the original source; the final gradient uses this independent source
  and the final restricted reference B/Exner fields. Nonfinite values fail even
  if their bit patterns happen to agree.
- `Makefile`: the new module and explicit build dependencies.
- `test/parallel_block_velocity/`: standalone strict numerical tests at `-O0`
  and `-O2`, exercising all 27 Qperp input positions, dense ordered summation,
  asymmetric geometry, all eight restriction-mask combinations, unused poisoned
  operands, edge signs/metrics, inactive cells, and 30 ordered physical layers.

The oracle adapter still uses Domain topology traversal and shared Domain
primitive inputs. That is intentional for this internal arithmetic check;
it is **not** a final-owner production plan. Physics is taken from the existing
measurement, not evaluated again by the new kernel. Pressure integration has
standalone tests only; it is not yet a native authoritative producer.

## Remaining implementation — required before cluster handoff

1. Compile a generation-scoped final-owner dependency plan. The Qperp stencil
   needs PV at neighboring edges; producing that PV from prognostic values
   requires a larger dependency neighborhood than merely importing Qperp's
   nine nodes. Audit coverage explicitly rather than assuming the existing
   scalar ghost plan suffices. Cache geometry weights, addresses and per-level
   service lists, with dynamic validity reset each RK stage.
2. Implement one authoritative shared primitive producer (mass flux, PV,
   pressure/geopotential, KE/Bernoulli and Exner) and complete its boundary and
   fine-to-coarse phases. Retained mass compatibility must consume this work
   rather than trigger duplicate producers. B/Exner gradient inputs are the
   final restricted fields, not Stage 174's dynamic Exner reconstruction.
3. Execute source kernels at source-time primitive phases; route native child
   sources and preserve mask/absent-child ordering. Execute the gradient only
   after final B/Exner restriction. Keep completed Domain values reference-only.
4. Replace the raw velocity tendency import and Domain velocity-source/gradient
   calls in production. Complete velocity RK boundary/scaffold values natively;
   `RK_sub_step_compatibility` must no longer read Domain velocity tendencies.
   Preserve Stage 174's accepted add/subtract rounding where needed for
   numerical continuity, using the native result rather than a Domain trend.
5. Prove the cut with production source/gradient/import counts of zero and
   poisoned or denied completed Domain velocity reads. Run checked RK3/RK4,
   optimized restart/checkpoint comparison, then full 83-rank oracles and paired
   performance tests. No tolerance relaxation or completed-source fallback.

## Local validation

Standalone tests:

```sh
sh test/parallel_block_velocity/run.sh
```

The local test artifacts are isolated under
`/private/tmp/wavetrisk-stage175.oCDe5s`. They do not modify the user's standard
restart test directory or `bin/climateJ5`. RK3 uses a separate source overlay;
the repository's `src/shared.f90` remains RK4.

Test results are recorded below after completion. These short oracle runs
exercise one timestep from checkpoint 3; they do not substitute for the full
restart/remap acceptance required for the eventual production replacement.

- Standalone checked `-O0` and `-O2` tests: PASS.
- Fresh RK3 and RK4 checked application builds: PASS with the normal 132-column
  limit, warnings as errors, bounds checks and invalid/zero/overflow traps.
- Fresh optimized RK4 application build: PASS. GNU Fortran reports a
  maybe-uninitialized array-descriptor warning for `air_temperature` in unchanged
  `src/parallel_block.f90`; there are no warnings in the new velocity module.
- Four-rank RK3 oracle, final source: PASS, exit 0, three RK substages. Last
  printed state: 0.3016 d, dt 63.4 s, Jmax 7, DOF 62816, relative mass 0.95627.
  Log: `rk3-oracle/RK3_debug_oracle.log` under the isolated artifact directory.
- First four-rank RK4 oracle: PASS, exit 0, four RK substages. Last printed
  state: 0.3021 d, dt 103.6 s, Jmax 7, DOF 62782, relative mass 0.95592. The
  final revision adds explicit nonfinite rejection and node diagnostics; a
  rerun of that final binary is recorded separately below.
- Final four-rank RK4 oracle: PASS, exit 0, four RK substages, with the same
  printed numerical state. Log: `rk4-oracle/RK4_final_debug_oracle.log`.

Final tested checked-binary SHA-256:

```text
20625283819521d8b306ee421c5251713e404fa2ee6bffb56e54138722c99a88  bin/stage175-rk3-check/climate
1f5e8963c20dc5f61333f44cbc58550491c53096b13740855dcc280259cbb51e  bin/stage175-rk4-check/climate
```

The checked tests require no new oracle flag. They used
`WAVETRISK_VALIDATE_BLOCK_DYNAMICS=1`,
`WAVETRISK_VALIDATE_BLOCK_ADAPTATION=1`,
`WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS=0`, and
`WAVETRISK_PROFILE_PARALLEL_BLOCKS=0`. These are local internal tests, not a
request for another cluster testing cycle.

No production-path removal, 83-rank test, paired performance result or complete
Stage 175 restart/checkpoint acceptance is claimed for this internal gate.
