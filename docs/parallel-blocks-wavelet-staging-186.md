# Stage 186: oracle scratch ordering and RK wavelet initialization

This stage continues local development with the frozen Stage 185 J4/J6
four-rank fixture. It tested two bounded changes in parallel_block_mpi.f90. Only the oracle
correction is retained; the zero-import experiment was not promoted.
No checkpoint input, numerical tolerance, physics parameter or integration
formula changes. The original legacy-screen and checked-versus-optimized
failures remain separate unresolved results.

## Oracle noninterference correction

The final RK wavelet transform first reconstructs a Domain-layout scratch
image and installs native velocity restrictions on the fixed coarse scaffold.
That scaffold is below the compact block roots. The complete vector wavelet
producer subsequently consumes those native coarse values.

Previously, the scalar-wavelet input oracle ran between restriction and vector
production. Its assert_block_domain_field_family_match call reconstructs the
same scratch arrays: it zeros the arrays and reseeds the fixed scaffold from
the independently accumulated Domain solution. This replaced the native coarse
restrictions with values differing in their last bits. Enabling validation
therefore changed a later producer's inputs.

An isolated checked diagnostic copied the scratch vector immediately before
that assertion and compared fixed coarse patches immediately after it, then
stopped intentionally. In the first final RK stage, the retained log reports
29 changed values on rank 3, 33 on rank 2, and 36 on rank 1, with maximum
3.5527136788005009e-15. Rank 0's diagnostic line was not retained in the merged
abort output, so no four-rank total is claimed. This direct intervention supports
the causal path suggested by Stage 185's level-4-only checkpoint mismatch.

The final-stage input oracle now runs as part of the initial reconstruction,
before native restriction writes the coarse stage. Provisional-stage input
checks stay in their existing location. Both use the same independent Domain
reference, integrated-field coverage and existing tolerance. No scratch copy,
new persistent storage, replay or relaxed checkpoint gate is introduced.

## Unpromoted experiment: reduce repeated zero-image imports

The native WT caller clears atmospheric Domain wavelets at each RK substage.
Previously, every substage transported that image back to final-owner blocks,
including unchanged soil/surface wavelets. The first substage still performs
that full import so each timestep starts with current non-atmospheric state.
Later production substages clear atmospheric block patch interiors directly.
The rest of the import lifecycle and post-production compression import remain.

The oracle retains the original full import at every substage, independently
checking the production initialization through the exact output comparison.
The new kernel does not touch non-atmospheric coefficients, solution arrays,
boundary/ghost payloads or rollback state. It uses the current field layout and
has no cache that could survive topology/restart changes. Layout tests exercise
mixed soil/atmosphere extents, variables and multiplicities, idempotence, and
rejection of malformed storage under both O0 and O2 with runtime checks.

For RK4 this removes three of the four zero-image imports per timestep; for
RK3 it removes two of three. These counts refer only to zero-image imports,
not all wavelet communication. Profiling adds rank-stage and cleared-value
counts to an existing reduction, with no new collective. End-to-end speedup
still requires repeatable, pressure-qualified timing evidence.

## Short numerical regression fixture

`test/parallel_block_profile/local_j4_wavelet_check.in` uses the same evolved
seed but time_end=0.3060 and dt_write=0.3050. The solver suppresses checkpoint
writes on its first step; this window instead crosses a save boundary after
that suppression (three RK3 steps, or two RK4 steps with the current seed).
It is a distinct numerical fixture, not a neutral truncation for timing:
time_end changes integer-clock scaling. Compared executables must use identical
inputs. Its checkpoint/reload comparison does not establish post-reload
stepping or remap coverage; the full RK4 fixture supplies that evidence.

## Validation and evidence

Optimized and checked RK4 builds succeeded and all 55 helper tests passed.
Full RK4 restart acceptance passes: optimized and checked production match
their respective Stage 185 references exactly, and checked oracle-on/off
checkpoints now match exactly. Each full run completes 12 steps, checkpoint 3,
its reload and a subsequent remap. The previous 208 differences are all gone.
The short RK3 regression also passes: three checked steps and checkpoint/reload
complete in both modes, with exact oracle-on/off checkpoint equality. This
short test has no post-reload timestep/remap coverage. Repeated timing
measurements are recorded separately; no timing gain is inferred from a
single validation run.

- Original-source diagnostic build: `/private/tmp/wavetrisk-stage186-oracle-diag-01/`
- Diagnostic run: `/private/tmp/wavetrisk-stage186-oracle-diag-run-01/`
- Diagnostic preparation: `/private/tmp/wavetrisk-stage186-diagnostic-build.py`
- Oracle-only fix build: `/private/tmp/wavetrisk-stage186-oracle-fix-01/`
- Combined candidate build: `/private/tmp/wavetrisk-stage186-wavelet-opt-01/`
- Acceptance runs: `/private/tmp/wavetrisk-stage186-acceptance-01/`
- Acceptance driver: `/private/tmp/wavetrisk-stage186-acceptance.py`
- Helper tests: `/private/tmp/wavetrisk-stage186-helper-tests.log`

At this stage, cluster acceptance was still pending. Subsequent Stage 187
cluster acceptance is recorded separately with job 54452.

## Timing outcome and promotion decision

| Measured run | Elapsed seconds | Median ordinary step seconds | VM screen |
| --- | ---: | ---: | --- |
| Stage 185 baseline 1 | 121.69 | 9.31 | pass |
| Zero-import candidate 1 | 125.04 | 9.76 | fail |
| Zero-import candidate 2 | 126.43 | 9.81 | fail |
| Stage 185 baseline 2 | 123.12 | 9.37 | fail |

A separate 120.66-second baseline warmup is excluded. All numerical gates
pass. The profile confirms 144 avoided rank-stage zero imports (3 per step,
12 steps, 4 ranks), but the experiment does not demonstrate an elapsed-time
gain. Its performance code and focused test are preserved in the temporary
candidate archive and removed from the proposed production changes. The oracle
ordering correction remains. The recorded source hashes describe the measured
experimental candidate, not the subsequent source after removing that experiment.

The profile identifies scalar restriction capture as a larger next target: about
13.22 seconds maximum-rank inclusive time over the 12 profiled steps, with
870,420 capture calls per rank in the two report windows. This timer includes
patch physics capture and boundary capture; it does not isolate patch writes. Avoided transfer
counts alone are not substituted for repeatable performance evidence.
