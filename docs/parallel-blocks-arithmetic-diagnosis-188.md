# Stage 188: checked/optimized arithmetic and legacy equivalence

The local J4 checked-versus-optimized discrepancy is eliminated by compiling
the optimized solver with `-ffp-contract=off`. Solver sources, inputs, ranks,
external physics object files and exact comparison criteria are unchanged.
The two-step and twelve-step J4 checkpoints then match their checked references
exactly. The longer run reloads checkpoint 3 and performs a subsequent remap.

GNU Fortran 15.2.0 reports `-ffp-contract=fast` for the original O2 build.
[GCC documents](https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html) that
`off` disables expression contraction, while `fast` permits fused operations
such as multiply-add. The local instruction audit finds fused instructions in
the optimized geometry, wavelet, physics interface and block MPI objects, and
none in those four checked or O2/no-contraction objects. The external physics
objects are byte-identical between the block arithmetic-control builds.

Fused and separate operations can round differently. This establishes a
compiler-arithmetic cause for the observed J4 cross-mode differences; it does
not establish that the checked trajectory is physically more accurate. It is
also distinct from the Stage 186 oracle-state correction and the cluster
checkpoint I/O corruption. No tolerance is loosened and no replay is used.

## J4 evidence

The original fixture is retained: min_level=4, max_level=6, four ranks, RK4,
30 atmospheric layers, 10 soil layers, restart from checkpoint 2. The short
window ends at 0.3060 and writes/reloads checkpoint 3 after two steps. The full
window ends at 0.33, takes twelve steps, and includes a remap after the reload.
Comparisons use identical inputs within each window, including time_end.

| Comparison | Exact checkpoint fields | Existing numerical screen |
| --- | --- | --- |
| Default optimized vs checked, short/full | Fail | Kept as the original cross-mode failure |
| Optimized with contraction off vs checked, short/full | Pass | Pass |
| Default legacy vs default block, full | Fail | Fail |
| Legacy vs block, both with contraction off, full | Fail | Pass |

For the full J4 legacy comparison, disabling contraction reduces maximum
atmospheric wavelet differences from 1.59481e-9 to 8.33111e-12 for velocity,
2.47872e-11 to 3.29692e-12 for mass, and 2.46045e-8 to 1.62254e-9 for
mass-weighted temperature. Headers, thresholds, topology and soil are exact.
The pre-existing screening bounds are 1e-10, 1e-10 and 1e-7 respectively.
Passing that screen is not bitwise legacy equivalence.

Both legacy executables are built from an archive whose tracked files verify
against authoritative main `b82647b36eb825aed0c0aae4e92adcf4ab46d608`.
The arithmetic control adds a separate, untracked make include containing only
`FFLAGS += -ffp-contract=off`; it does not edit legacy numerical sources.

## Original J5 failure fixture

The original min_level=5/max_level=6 fixture was also checked with four
ranks, resume=4 and time_end=0.3350. It requires checkpoints 5 and 6, both
reloads, and a following remap for each. The checked reference is the completed
Stage 184 production run. Both new J5 runs complete all eleven steps, both
checkpoint/reload cycles and both following remaps.

The optimized block solver with contraction off matches the checked reference
exactly at checkpoints 5 and 6. This confirms the compiler-arithmetic cause on
the original failure fixture, where the default optimized versus checked run
previously differed by up to 6.62949e-4 in velocity and 0.1962742 in
mass-weighted temperature at checkpoint 6.

Legacy versus block, both with contraction off, passes the unchanged numerical
screen at both checkpoints. It still fails exact equality. Checkpoint-6
atmospheric wavelet differences are:

| Field | Original default legacy/block comparison | Both with contraction off | Existing bound |
| --- | ---: | ---: | ---: |
| Velocity | 5.99723e-9 | 1.22107e-11 | 1e-10 |
| Mass | 5.70420e-9 | 4.32010e-12 | 1e-10 |
| Mass-weighted temperature | 2.78759e-6 | 3.31420e-9 | 1e-7 |

Headers, thresholds, topology, soil and poles are exact, and no compared values
are nonfinite. The 23 external physics objects are byte-identical between
legacy and block. Legacy tracked sources still verify against authoritative
main after both J4 and J5 builds. The original failed default-build comparison
remains a failure; the passing result belongs to the explicitly controlled
arithmetic configuration. The Stage 182 midpoint-crossing diagnosis remains
valid for the original trajectory. Residual nonzero legacy differences are
bounded by the existing screen on these fixtures, not proved harmless for all
future trajectories or traced to a unique source expression here.

## Reproduction and build option

`FP_CONTRACT=off` is now the default following bbserv job 54465.
`FP_CONTRACT=default` explicitly selects the original compiler behavior for
comparisons. Use `FP_CONTRACT=off` explicitly when documenting matched-arithmetic
diagnostics. `fast` is
also accepted; unknown values fail at make parsing. The option affects solver
code, including its physics interface; the external physics package keeps its
own flags. Preserve identical external physics objects in a comparison.

Use fresh build/output directories when changing compiler settings:

```sh
make -j1 PARAM=param_J4 DEBUG=false FP_CONTRACT=off \
  BUILD_DIR=build-j4-fp-off BIN_DIR=bin-j4-fp-off
make -j1 PARAM=param_J4 DEBUG=check FP_CONTRACT=off \
  BUILD_DIR=build-j4-check-fp-off BIN_DIR=bin-j4-check-fp-off
```

The existing `local_j4.py` runner and exact checkpoint reader reproduce the J4
comparisons. For the original J5 window, build with PARAM=param_J5 and use:

```sh
python3 test/parallel_block_profile/run_arithmetic_diagnostic.py \
  --block /absolute/path/to/block-fp-off/climate \
  --legacy /absolute/path/to/legacy-fp-off/climate \
  --checked-reference /absolute/path/to/completed-checked-run \
  --fixture /absolute/path/to/original-j5-seed \
  --out /absolute/path/to/new-diagnostic-output
```

The driver runs sequentially, verifies fixture identity and both post-restart
remaps, and records exact and existing-screen comparisons separately. It rejects
unmatched reference inputs or rank counts before running. The local diagnostics
do not by themselves establish cluster equivalence or performance. Subsequent
job 54465 provides that evidence for the 83-task short fixture. Job 54452 timings
remain evidence for the original compiler-default builds.

Raw runs, compiler logs and instruction counts are under
`/private/tmp/wavetrisk-stage188-crossmode-01`. The optimized arithmetic-control
build is `/private/tmp/wavetrisk-stage188-nofma-01`; the verified legacy copy
is `/private/tmp/wavetrisk-stage188-legacy-j4-01`. All 64 helper tests pass,
including three new tests enforcing the original two-restart coverage.

The compact source/binary identities, comparisons and coverage are archived in
[the Stage 188 result record](parallel-blocks-arithmetic-diagnosis-188-results.json).
The Stage 188 commit includes the build option, drivers, tests and evidence.
[Cluster job 54465](parallel-blocks-cluster-arithmetic-188-54465.md) passed all
mandatory gates with a 0.103% observed elapsed increase and unchanged median
ordinary-step time. On that evidence the Makefile default is now
contraction off; numerical Fortran sources remain unchanged.
