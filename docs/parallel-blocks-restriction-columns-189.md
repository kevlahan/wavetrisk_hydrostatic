# Stage 189: column restriction interpolation

Prepared against Stage 188 (`61509429`) with `FP_CONTRACT=off`. Cluster
acceptance on the tested single-node `nk` configuration is recorded below. This stage implements the first bounded portion of the
[performance plan](parallel-blocks-performance-targets-189.md).

## Production change

Compact production scalar restriction now resolves interpolation inputs across
physical columns, instead of locating and decoding every scalar record at each
vertical layer. Parent edge masks are read once per horizontal parent node,
using the compact storage contract that shares geometry across physical layers.
Inactive parents skip interpolation. All consumed fluxes and weights retain
finite/poison validation through the checked column-storage interface.

The interpolation formulas and ordered three-term accumulations are preserved.
For each active parent node with 30 physical layers, interpolation makes 60
column location calls instead of 1,800 scalar location calls. Construction of
its two horizontal index sets is also performed once per column. This is an
operation-count reduction, not a 30-fold speedup prediction.

Only parent direct flux is written; the inputs belong to the finalized finer
level. Node/layer loop interchange therefore does not move reads across writes
to the same fine-level flux. Level traversal, MPI exchanges, boundary closure,
partial/coarse-flux formulas, soil exclusions, field precision and tolerances
are unchanged. Mass keeps its native production path. Oracle mode retains the
original scalar restriction path and independent full records.

This first implementation batches interpolation; partial/coarse closure and
parent writes still use the original scalar routines. It does not complete the
entire restriction-column roadmap. Invocation-local column arrays add bounded
scratch storage, with no persistent numeric cache or new collective. Checked
Fortran reports a four-value array temporary when passing a strided row of the
interpolation result to the unchanged coarse-flux routine; no checks are
suppressed. That small temporary is included in solver execution costs.

The existing profiling reduction carries three additional counters:
`restriction interpolation columns: columns locations replaced scalar reads`.
The last counter is the equivalent scalar read count, not measured memory
traffic. Original scalar access and oracle routines remain available.

## Local validation and measurement

Use the evolved J4/max-level-6 fixture with four local ranks first. This catches
indexing, arithmetic and restart regressions without relying on cluster access.
Use separate fresh optimized and checked build directories. External physics
objects are unchanged; checked solver flags do not instrument that separately
built package.

The new component test compiles the actual old and new interpolation functions,
actual restriction-index generation and scalar storage module. A synthetic
location resolver exercises patch, boundary and ghost classes. At O0 and O2,
with contraction off, it requires bitwise equality across physical layers,
variables, node positions and both compact/full storage. Full-record geometry
varies by layer to check that the column reader does not assume oracle geometry
is shared. Negative cases reject poisoned, nonfinite and out-of-range columns.

The associated benchmark uses this synthetic routing and actual interpolation
and storage routines. It alternates eight trials, 30,000 column evaluations per
trial, requiring matching checksums. It excludes actual MPI and block route
lookup, and does not represent solver memory pressure or whole-step speed.

Run helpers:

```sh
python3 -m unittest discover -s test/parallel_block_profile -p 'test_*.py'
python3 test/parallel_block_profile/benchmark_restriction_columns.py \
  --out /absolute/path/to/new-component-results
```

Run the component benchmark when solver runs and builds are finished. The
benchmark records compiler arguments and system VM/RSS evidence. Solver
correctness runs, profiling and unprofiled timing must remain distinct.

## Local results

All 73 helper tests pass (including the subsequent cluster build/submission helpers). Clean optimized and checked RK4 builds and the
isolated checked RK3 build succeed. All 23 external physics objects match the
frozen Stage 188 build. Only `parallel_block_mpi.f90` changes among the 112
tracked numerical source files; the isolated RK3 build additionally selects RK3
in its own `shared.f90`, without modifying the production checkout.

- Full optimized RK4: 12 steps, checkpoint 3/reload and following remap; saved
  checkpoint fields match Stage 188 contraction-off output exactly.
- Full checked RK4: the same 12-step coverage; checkpoint fields match the
  archived independent full oracle result exactly. No new full oracle run was
  needed because its scalar path is unchanged.
- Short checked RK3: three steps and checkpoint/reload; matches the archived
  independent RK3 oracle checkpoint exactly. This short case does not cover a
  subsequent remap.
- Separate optimized two-step profile: checkpoint matches the Stage 188 short
  reference exactly. Across four ranks, 22,800 interpolation columns use
  1,368,000 location calls, replacing 41,040,000 equivalent scalar reads.
- The unchanged local legacy comparison screen passes. This is not bitwise
  legacy equality or new cluster legacy acceptance.

The final eight-trial component benchmark has median scalar/column wall times
of 0.386638/0.167393 seconds, a 56.7% reduction in this isolated operation. All
checksums agree. Its VM screen passes: 2.67 MiB system-wide decompression, no
compression or swapping, and 3.66 MiB sampled process RSS. Earlier benchmark
trials agree in direction; one sandboxed monitoring attempt could not sample
RSS and is not the final qualified measurement.

No whole-solver timing gain is claimed locally. The optimized correctness run
exceeded the conservative VM screen (31.66 GiB system-wide decompression) and
overlapped a checked build; its elapsed time is validation evidence only.
Checked and profiling runs are also excluded from performance comparisons.

Durable results, gates and source/binary identities are in
[the result record](parallel-blocks-restriction-columns-189-results.json).
Raw outputs are under `/private/tmp/wavetrisk-stage189-*`; the final component
benchmark is `wavetrisk-stage189-interpolation-benchmark-04`. Builds and raw
outputs are temporary local artifacts, not committed solver data.


## Cluster acceptance

After local checks pass, build the Stage 188 baseline and Stage 189 candidate
with `PARAM=param_J5`, `FP_CONTRACT=off`, fresh build directories and the same
compiler/MPI/external physics settings. Build a checked Stage 189 executable
separately. Preserve the baseline source and executable identities before
applying the candidate patch. Do not reuse Stage 187 default-contraction
binaries as the baseline for this test.

Keep the existing J5/max-level-7 checkpoint-3 fixture, 83 ranks, `dt_write=0.01`
and `time_end=0.315`. Do not overwrite the frozen fixture. If those parameters
need changing, prepare a new fixture directory and record its identity.

The v2 kit adds `build_restriction_189.py`, which verifies the Stage 188 source
and Makefile hashes in a fresh checkout, prepares physics dependencies, builds
an optimized baseline, applies the supplied patch, then builds optimized and
checked candidates. Every build uses explicit J5 and contraction-off flags.
It verifies unchanged sources during each build and identical 23-object external
physics inputs, and saves `restriction-189-builds.json`. Existing build paths
are rejected. No clean/reset/commit operation is used.

`cluster_189.sbatch` now takes four arguments: kit directory, prepared source
tree, immutable fixture, and new output directory. Use the v2 kit consistently;
the original seven-argument launcher is superseded. The submission helper
verifies the three recorded executable hashes before starting the campaign and
passes the complete build manifest as build notes.

The campaign runs checked production/oracles, warmups, three timing pairs and
two separate profiles. Every run must reload checkpoint 4 and subsequently
remap. It retains exact baseline/candidate and oracle comparisons and adds
mandatory checked/optimized comparisons. A failed additional gate makes the
final result failed. Three new tests cover build order/flags, changed physics,
stale build directories and changed submission binaries. The physics dependency
setup also compiled successfully from a fresh local package checkout.

Compare ordinary steps and elapsed time separately; exclude warmups, checked
runs and profiles from timing. Cluster results are needed to establish a
repeatable 83-task gain. Further inverse/ghost and boundary-buffer changes
remain separate future candidates.

## Cluster follow-up, 14 September 2026

Job 54549 on `nk` passed the supplied exact checkpoint gates and showed 1.24%
lower median elapsed time and 1.96% lower ordinary-step median. The supplied compilation lines confirm that both optimized builds used
`-march=native -ffp-contract=off`; the checked build used `-O0 -ffp-contract=off`.
The earlier flag discrepancy is resolved.
See [the campaign analysis](parallel-blocks-cluster-189-54549.md) for results,
profile attribution and the precise acceptance limits.
