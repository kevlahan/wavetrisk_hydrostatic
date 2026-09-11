# Stage 184: full local checked restart acceptance

The cluster is offline. This local acceptance set tests the unchanged Stage 184
column-divergence candidate on the evolved J6 checkpoint-4 fixture with four
MPI ranks, 30 atmospheric and 10 soil levels plus surface. Full J7/83-task
cluster acceptance and matched performance measurements remain pending.

## Coverage and gates

Run checked RK4 production, checked RK4 all-oracle, checked RK3 production, and
checked RK3 all-oracle sequentially to time_end 0.3350. Require checkpoints 5
and 6, a reload of each, and a remap after each reload before the next reload or
termination. Bounds checks, signalling-NaN initialization, FPE traps, strict
132-column and warning-as-error solver builds remain enabled. The external
physics package retains its normal build flags; this is checked solver coverage,
not new instrumentation of the external physics library.

RK4 reuses the verified Stage 184 checked binary. RK3 is a fresh isolated
checked build of the same candidate, changing only the timeint_type selection.
The working source remains RK4. No copied parent executable is used for RK3.

Require exact checkpoint comparison for checked RK4 against the completed
optimized Stage 184 run, and for oracle-on against oracle-off within each RK
scheme. The optimized Stage 184 RK4 checkpoints already match Stage 183 exactly.
These comparisons do not establish independent legacy equivalence or a new
RK3-versus-Stage183 baseline comparison. The original unmodified legacy fixed
screen still fails; the Stage 182 causal diagnosis is unchanged and replay is
not linked into production.

## Why retain time_end 0.3350

A shorter endpoint reduces steps but cannot remove required restart/remap
coverage. On the existing RK4 trajectory, checkpoint 5 is at 0.32004 days,
checkpoint 6 at 0.33202, and the step following the second reload/remap ends at
0.3344. Ending just after checkpoint 6 would save only the final additional
step compared with the existing 0.3350 input (which finishes at 0.3368).

Moreover, time_end is not purely a loop-stop setting: main_mod::init_structures
uses it to set time_mult, and checkpoint records include itime. Restart resets
itime from the stored physical time and the new multiplier. Thus changing
time_end changes integer-clock scaling and can change headers and rounding;
it is not safe to assume exact equivalence to existing reference checkpoints.
Keep the established endpoint for this full acceptance matrix. Short first-step
smoke checks remain useful for future edits, but do not replace this coverage.

## Memory interpretation

Stage 183 reduced scalar workspace allocated capacity from 6.517 to 2.329 GiB
and the four-module counted total from 10.487 to 6.300 GiB at the aligned phase.
That reduction is established capacity evidence, not a whole-process RSS peak.
Stage 184 changes access/computation and adds small automatic column buffers;
it is not another substantial persistent-memory reduction.

Historical full-run observations:

| Run | System-wide compressed traffic | System-wide decompressed traffic |
| --- | ---: | ---: |
| Corrected block baseline, Stage 180 | 507.62 GiB | 467.57 GiB |
| Stage 183 compact production | 99.76 GiB | 87.26 GiB |

Those are cumulative traffic counters, not resident compressed footprints.
Different machine conditions prevent a controlled attribution, and Stage 183
still failed the existing local timing-quality screen for heavy compression.
Stage 184's earlier validation had no memory monitor, so it supplied no evidence
of further compression improvement.

This set records existing per-run VM counter deltas and sampled simultaneous
child RSS. A separate five-second VM sampler records logical pages stored in
the compressor and physical pages occupied by it. These are system-wide:
other applications contribute, and child RSS excludes compressed footprint.
The additional sampler began after the first run launched, so its sampled peaks
are not guaranteed whole-run peaks. No background application is closed or
modified to manufacture a pressure-free result.

## Evidence locations

- Full-run results: `/private/tmp/wavetrisk-stage184-full-checked-01/results.json`
- VM samples: `/private/tmp/wavetrisk-stage184-full-checked-01/system-memory.json`
- Reproduction driver: `/private/tmp/wavetrisk-stage184-full-checked.py`
- VM sampler: `/private/tmp/wavetrisk-stage184-vm-sampler.py`
- RK3 build: `/private/tmp/wavetrisk-stage184-column-rk3-check-01/`
- RK4 build: `/private/tmp/wavetrisk-stage184-column-opt-01/`

Drivers require fresh output directories on reuse. Source/binary hashes and
final validation status are recorded in the adjacent results JSON when the
matrix completes. No production source or numerical tolerance changes are
part of this acceptance work.

## Initial checked-mode finding

The first full checked RK4 production run completes all 11 steps, both reloads
and both following remaps without runtime-check failures. However, exact
comparison against the optimized executable fails at both checkpoints. Headers,
topology and thresholds remain exact, but checkpoint-6 atmospheric wavelet
maxima differ by velocity 6.62949e-4, mass 6.03806e-4 and mass-weighted temperature
0.1962742; soil temperature also differs. This cross-build-mode gate is retained
as a failure. It is not the previously diagnosed legacy-versus-block discrepancy
and has not been explained by that diagnosis.

Before continuing the matrix, a full Stage 183 checked RK4 control is run with
the identical fixture, endpoint and oracle-off settings. Its source differs
from the checked candidate only in the two Stage 184 production files (apart
from nonexecuted tracked symlink placeholders absent from the old archive).
The control will determine whether the candidate introduced the checked-path
change. Passing a same-build-mode control will not relabel the failed
checked-versus-optimized comparison as passing.

Continuation evidence is under
`/private/tmp/wavetrisk-stage184-full-checked-02/`, with driver
`/private/tmp/wavetrisk-stage184-full-checked-continue.py`.

The matched Stage 183 checked control passes exact candidate comparison at
checkpoints 5 and 6, including all compared fields, soil, headers, thresholds
and topology. Thus the checked-versus-optimized difference is already present
before Stage 184 on this fixture; its numerical cause is not yet diagnosed.
All 23 physics object files are byte-identical between those builds (the static
archive container hashes differ). The remaining full oracle/RK3 cases continue
against the matching checked-mode production references. The initial cross-mode
failure remains in the results rather than being discarded.

## Superseded by the local J4/J6 redesign

At the user's request, the expensive matrix was stopped during the RK4 oracle
run after two completed steps. Its MPI descendants and VM sampler were stopped;
no RK3 full run was started. Full-matrix acceptance is not claimed. The completed
RK4 production coverage and same-mode Stage 183 equality remain valid evidence.
A fresh param_J4 / max_level6 fixture is prepared separately in Stage 185,
retaining the original cluster fixture unchanged.
