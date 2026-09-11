# Stage 185: local J4/J6 fixture and acceptance protocol

The user clarified that the cluster fixture is intended to restart from
checkpoint 3 near t=0.3 with min_level5/max_level7, take about eight steps,
exercise one checkpoint/restart and end at 0.3225. Its dt_write=0.01 belongs
to that short restart test, not to a general local spin-up configuration.
The cluster fixture remains unchanged. This stage redesigns local testing;
it does not change the Stage 184 production sources or numerical tolerances.

## Selected local configuration

| Setting | Value |
| --- | --- |
| Build parameter | param_J4 (min_level=4, DOMAIN_LEVEL=1, 40 Domains) |
| Maximum adaptive level | 6 |
| Local MPI ranks | 4, one thread per rank/library |
| Physics / atmospheric / soil levels | Simple / 30 / 10 plus surface |
| Tolerance | 0.025, retained from the existing local fixture |
| Genuine evolved seed | localJ4J6_checkpoint_0002.bin.zst, t approximately 0.30205 days |
| Restart selection | resume=2 |
| End time | 0.33 days |
| Output interval | dt_write=0.02 days, CP_EVERY=1 |
| Required new restart | checkpoint 3, followed by a remap before completion |

The proposed input is `test/parallel_block_profile/local_j4_restart.in`.
The immutable local seed directory is
`/private/tmp/wavetrisk-stage185-j4-fixture-02/seed`.
Seed SHA256:
`0272c3717e3ed520d58a542707cc886e85a9412e2b5d61f20ddcb96dd007f3f4`.
Do not substitute or rename the original 160-Domain J5 checkpoint: the J4
layout has 40 Domains and needs its independently evolved seed.

The endpoint and save interval are now frozen across compared binaries. As
previously noted, time_end also controls integer-clock scaling, so changing
it is a fixture change, not a guaranteed neutral truncation of an old run.
The J4 comparison helper requires identical inputs, seed, grids and rank count.

## Refinement evidence and why this seed was selected

The one-time optimized fresh pilot used `local_j4_spinup.in`: resume=-1,
time_end=0.35 and dt_write=0.15. It produced checkpoints 1 and 2 after the
initial checkpoint-0 bootstrap. It first reported Jmax5 near 0.163 days and
Jmax6 near 0.218 days. The bootstrap checkpoint is not counted as a timestep
restart cycle. No unmodified legacy executable or replay intervention is used.

Neither Jmax, raw patch capacity, finest coefficient magnitude nor complete
triangles alone adequately measures this edge-heavy adaptive workload. The
checkpoint has 1,436 stored level-6 patches, yet its exported triangles are
levels 4/5 and its finest wavelet coefficients are below threshold. Direct
existing `cpt_dt` active-mask counters resolve the ambiguity:

| Level | Active nodes | Active edges |
| --- | ---: | ---: |
| 4 | 2,570 | 7,710 |
| 5 | 4,702 | 18,584 |
| 6 | 844 | 7,886 |

Thus 8,730 of 42,296 active DOFs, **20.64%**, are already at level 6. Counts
use the solver's existing mask>=ADJZONE definition, shared horizontally over
layers. They include the adaptive stencil/adjacency work counted by the
standard DOF report; they are not claimed to be only TRSK physics points.

A later candidate checkpoint near 0.3524 has 14,836 finest DOFs (26.58%) but
55,809 total DOFs. The earlier seed gives substantial finest-level computation
at lower cost, so it was selected. The later exploratory extension is recorded
by `local_j4_seed_extension.in` and its raw results but is not required to
regenerate the selected seed.

The counter census is an isolated read-only source copy: it only prints the
already-computed local arrays, without new numerical calls or collectives.
Its executable is diagnostic, never a timing/production binary. The production
source stays unchanged. VTK/threshold summaries are auxiliary diagnostics,
not substitutes for these counters.

## Running and comparing

Use fresh isolated build directories, with serial make and PARAM=param_J4.
RK4 checked and optimized solver objects must use separate BUILD_DIR/BIN_DIR.
RK3 uses an isolated source copy with only the integrator selection changed;
record actual source and executable hashes. The external physics dependency
Makefile.inc must exist (as in the verified local archives). Fresh source-only
copies lack that generated file; restore it from verified unchanged physics
sources or regenerate the package dependencies before building.

Example using the verified local optimized candidate archive:

```sh
python3 test/parallel_block_profile/local_j4.py \
  --build /private/tmp/wavetrisk-stage185-j4-opt-01 \
  --fixture /private/tmp/wavetrisk-stage185-j4-fixture-02/seed \
  --out /private/tmp/wavetrisk-j4-new-production-run
```

For checked runs add `--binary bin-checked/climate`. Add `--oracles` to enable
all three independent validation oracles. Production clears all inherited
WAVETRISK switches; all three validation switches and profiling/diagnostics
are off by default. `--profile` selects a separate detailed profile. Keep all
outputs fresh siblings outside the immutable seed. The runner copies the
selected executable as climateJ4 and verifies min_level4/DOMAIN_LEVEL1 in the
actual run log. It records input, grid, seed, source and binary identities.

`compare_windows(..., expected_new=(3,))` requires exactly the one newly written
checkpoint, a reload and following remap, completion, and exact compared fields,
headers, topology and thresholds. This one-cycle gate is explicit to the new
local fixture; the existing J5 two-cycle helper and all numerical thresholds
remain unchanged. J4 decoding requires explicit expected_domains=40, preventing
accidental acceptance of the wrong checkpoint layout.

The bounded acceptance set is two matched optimized builds (Stage 183 and 184),
then checked RK4/RK3 each off/on. Within each checked integrator, compare oracle
on/off exactly. Do not demand or imply RK3/RK4 equality. The earlier checked
versus optimized failure is retained separately, and Stage 184's exact checked
Stage 183 control on the old J5/J6 fixture remains documented. Neither that
control nor the new local tests establishes independent legacy equivalence.

## Memory and timing interpretation

The purpose is a smaller, sufficiently refined local workload, not a forecast
of 83-task scaling. Compare identical seeds/inputs and compatible builds. Keep
checked and oracle runs out of production timing comparisons. Report per-step
solver time separately from startup/output/restart elapsed time.

The runner records system-wide compression/decompression and swap deltas plus
sampled simultaneous simulation-child RSS. RSS excludes compressed pages;
system counters include other applications. A lower total alone does not prove
a pressure-free comparison. Preserve the existing local timing-quality screen;
if compression remains substantial, report it and do not claim a clean speedup.
Authoritative main is unchanged, diagnostic binaries are excluded, and the
original unmodified legacy fixed-screen failure is not waived.

## Measured acceptance results

All four completed runs take 12 RK4 steps, write exactly checkpoint 3, reload
it, perform a subsequent vertical remap and finish near t=0.3308. The optimized
Stage 184 candidate matches the Stage 183 control exactly at checkpoint 3.
Checked production and checked oracle runs both complete without a runtime
check or internal oracle failure. The checked scope covers solver objects;
the external physics dependency retains its normal compilation flags.

**Full checked acceptance fails:** oracle-on versus oracle-off differs in 208
atmospheric velocity wavelet values, maximum absolute difference
7.105427357601002e-15. All other compared fields, headers, topology and thresholds
match exactly. This small magnitude does not satisfy the exact gate and is not
waived. Its cause remains undiagnosed. RK3 production and oracle cases were not
run after this failure; their checked executable is built and available.

Checkpoint equality tests the saved state before reload. Completion and the
following remap establish execution coverage after reload; this one-cycle
fixture does not compare a second post-reload checkpoint. Coverage self-checks
are not independent numerical controls.

| Run | Elapsed (s) | Median non-checkpoint step (s) | Compression (GiB) | Decompression (GiB) |
| --- | ---: | ---: | ---: | ---: |
| Optimized Stage 183 | 123.8 | 9.43 | 4.72 | 2.05 |
| Optimized Stage 184 | 119.7 | 9.17 | 3.71 | 2.10 |
| Checked RK4 production | 407.2 | 31.9 | 5.57 | 4.08 |
| Checked RK4 oracle | 1108.7 | 87.5 | 412.18 | 387.47 |

The optimized candidate peaks at 6.76 GiB sampled simultaneous child RSS.
Ordinary local testing is substantially more manageable, but all four runs
fail the retained timing-quality screen. In particular, independent oracle
storage still causes heavy compression. The 123.8 versus 119.7 second result
is an observation, not a clean speedup claim. No cluster scaling conclusion
is drawn. When bbserv returns, retain its standard 83-task fixture for a basic
matched optimized timing comparison with oracles disabled; this will not close
the failed numerical acceptance gates.

The 54 helper tests passed. The explicit-window helper's focused tests also
pass after requiring callers to supply checkpoint indices, avoiding a stale
default from the discarded later-seed candidate. Production sources are
unchanged from Stage 184. This protocol was validated before the retained changes were committed.

## Evidence locations

- Spin-up: `/private/tmp/wavetrisk-stage185-j4-pilot-01/`
- Later exploratory seed: `/private/tmp/wavetrisk-stage185-j4-refinement-01/`
- Selected-seed active counters: `/private/tmp/wavetrisk-stage185-j4-early-census-01/active-levels.json`
- Later-seed active counters: `/private/tmp/wavetrisk-stage185-j4-census-run-01/active-levels.json`
- Acceptance: `/private/tmp/wavetrisk-stage185-j4-acceptance-01/`
- Reproduction driver: `/private/tmp/wavetrisk-stage185-j4-acceptance.py`

Measured results and remaining limitations are recorded in the adjacent Stage 185
results JSON. Large checkpoints and binaries remain local temporary artifacts;
the inputs, helpers and recorded hashes provide the durable reproduction path.
