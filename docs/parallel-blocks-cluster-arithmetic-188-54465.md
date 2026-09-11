# Stage 188 cluster acceptance: job 54465

The user supplied the complete campaign results.json on 2026-09-11. The
[original file](parallel-blocks-cluster-arithmetic-188-54465-campaign.json) is
archived unchanged, SHA256
`c5d8117b925bcdd1faa95541fe1dba9fc72a1f51018ea10267f98271f2f07105`.

All ten runs passed on 83 ranks, bb[01-02], with 42/41 tasks, GNU Fortran
13.2.0, PARAM=param_J5, max_level=7, RK4, time_end=0.315 and dt_write=0.01.
Each run completed twelve steps, reloaded checkpoints 3 and 4, and reported
a remap after the new checkpoint-4 reload. The archived manifest matches all
112 checked source files and the original arithmetic-option Makefile. Every
run's executable hash matches its declared build and the recorded build hash.

The mandatory gates include eight independent comparisons and two reference
self-checks. All compared checkpoint headers, topology, thresholds and sixteen
field categories are exact, with no nonfinite values. Optimized contraction-off
runs compare against checked-off. Default optimized repeats compare against
their own warmup. Checked oracle-on compares against checked oracle-off.

The separate default-versus-off diagnostic is false, as recorded. That is a
comparison between different arithmetic policies, not a failed mandatory gate.
There is no new cluster legacy comparison in this campaign.

| Median metric | Compiler default | Contraction off | Off increase |
| --- | ---: | ---: | ---: |
| Elapsed time | 30.541876 s | 30.573341 s | 0.1030% |
| Ordinary timestep | 1.53 s | 1.53 s | 0.0% |

Three measured runs per mode used alternating order; warmups and checked runs
are excluded. The elapsed median difference is 0.031465 s, smaller than the
observed run-to-run ranges (0.099921 s default, 0.056871 s off). This supports
no material performance penalty observed on this short case. It does not
establish zero cost on every workload. Distributed memory pressure was not
automatically qualified by the harness.

## Default change

The Makefile now selects `FP_CONTRACT=off` when the option is omitted, for both
optimized and checked builds. This chooses the explicit setting already tested
locally and on bbserv. `FP_CONTRACT=default` retains the compiler's own behavior
for historical comparisons; `fast` is also available. No Fortran numerical
sources, physics-package flags, tolerances or checkpoint formats change.

Use fresh build and binary directories when moving an existing checkout to
the new setting: make does not infer that existing objects need rebuilding
just because compiler flags changed. Existing executable files remain as built.
The Stage 188 comparison helper still passes its settings explicitly, so its
default/off comparison retains the same meaning after this default change.
The original kit patch and tested Makefile hash remain frozen; the helper's
manifest additionally recognizes the promoted Makefile hash.

The Stage 188 commit records the default change, diagnostics, harness and
cluster evidence together. The earlier job 54452
Stage 187 speedup remains evidence for its original compiler-default builds;
this new comparison measures the cost of turning contraction off in Stage 187.

Validation of the changes: all 67 helper tests pass, including actual
Makefile flag selection for omitted and explicit arithmetic settings in both
optimized and checked modes. All 112 recorded Fortran source hashes remain
unchanged. The promoted Makefile hash matches the helper manifest.
