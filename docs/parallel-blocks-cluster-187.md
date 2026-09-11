# Stage 187 cluster test kit

The reproducible pre-187 baseline is committed as
`6d5c8a3a14cbba71eb53324e9a6f2c7526ec087f`. It adds Stage 184 column
divergence and the Stage 186 oracle correction to Stage 183, excluding Stage
187 capture and the rejected zero-import experiment. The Stage 187 candidate
is `073713bcdd972cc926b1b86247d2f8c9491d1db6`.

Both comparison versions must also receive the identical checkpoint I/O fix
from `checkpoint_io_fix.patch` before building, as described in
[the workaround instructions](parallel-blocks-checkpoint-io-workaround.md).
The source commits alone precede that common fix. The latest validation branch
contains the candidate, common fix, harness and evidence together.

After fetching these commits, create comparison trees directly:

```sh
git worktree add --detach ../wt-pre187 6d5c8a3a14cbba71eb53324e9a6f2c7526ec087f
git worktree add --detach ../wt-stage187 073713bcdd972cc926b1b86247d2f8c9491d1db6
git -C ../wt-pre187 apply "$kit/checkpoint_io_fix.patch"
git -C ../wt-stage187 apply "$kit/checkpoint_io_fix.patch"
```

Set `kit` to the extracted kit directory containing the common fix and runner.
Use new tree names if these directories already exist. The original patch-based
setup below is an alternative for repositories that have only Stage 183; do not
apply its baseline/candidate patches again to the commits above.

The kit contains two independent patches against that same Stage 183 commit:
`pre187-baseline.patch` and `stage187-candidate.patch`. Apply each to a separate
clean source tree; do not apply the candidate patch on top of the baseline.
The accompanying JSON manifests record the base commit, patch hash and resulting
Fortran source hashes. Both patches were reapplied locally and reproduce those
hashes exactly. The runner does not commit, build, submit jobs or change sources.

## Prepare the executables

For the original patch-based setup, transfer the entire extracted kit to bbserv.
Set `kit` to its absolute path and run the following from an existing repository
containing Stage 183:

```sh
kit=/absolute/path/to/wavetrisk-stage187-cluster-kit
git worktree add --detach ../wt-pre187 99571e6129261c36804f41d98ebf1dc41bd4f626
git worktree add --detach ../wt-stage187 99571e6129261c36804f41d98ebf1dc41bd4f626
git -C ../wt-pre187 apply "$kit/pre187-baseline.patch"
git -C ../wt-stage187 apply "$kit/stage187-candidate.patch"
```

Use the usual cluster build environment and physics setup in both trees.
`src/physics/Makefile.inc` is generated/ignored and is not supplied by a clean
Git worktree; the kit's `physics-Makefile.inc` is the portable source-list include
used for these builds. Install it in both trees before compiling:

```sh
cp "$kit/physics-Makefile.inc" ../wt-pre187/src/physics/Makefile.inc
cp "$kit/physics-Makefile.inc" ../wt-stage187/src/physics/Makefile.inc
bash "$kit/prepare_physics_dependencies.sh" ../wt-pre187 ../wt-stage187
make -C ../wt-pre187 -j1 PARAM=param_J5 DEBUG=false BUILD_DIR=build-opt BIN_DIR=bin-opt
make -C ../wt-stage187 -j1 PARAM=param_J5 DEBUG=false BUILD_DIR=build-opt BIN_DIR=bin-opt
make -C ../wt-stage187 -j1 PARAM=param_J5 DEBUG=check BUILD_DIR=build-check BIN_DIR=bin-check
```

There are two different generated includes: `src/physics/Makefile.inc` selects
physics source files, while `src/physics/simple_physics/phyparam/Makefile.inc`
orders the physics modules. The preparation script generates the latter in both
trees. Without it, a clean serial build can fail with missing
`read_param_mod.mod`. The script avoids the package’s `clean` target, which also
reformats source files. This setup omission was corrected after the first
cluster build attempt. A fresh local physics library and driver build passes
with the generated dependency include.

The resulting executables are `wt-pre187/bin-opt/climate`,
`wt-stage187/bin-opt/climate` and `wt-stage187/bin-check/climate`. Keep the same
RK scheme, compiler, MPI and physics options. The current candidate uses RK4;
run RK3 separately if required. The Makefile's optimized flags include
`-march=native`: build on a matching compute-node CPU or use the site's known
compatible target flags consistently in both optimized builds. Record the
actual flags and versions in the build notes. External physics objects use
their own flags; the checked solver does not imply checked external physics.

## Fixture and job

Use **J5/J7 with 83 tasks**. The current `param_J4` documents a 40-rank limit.
Prepare a fixture directory containing:

- The original `simple.in`: `max_level=7`, `resume=3`, `time_end=0.315`,
  `dt_write=0.01`, `CP_EVERY=1`, Simple physics, 30 atmospheric and 10 soil levels,
  `NCAR_topo=F`, `sso=F`. Keep its existing `run_id` and other parameters.
- Its genuine `<run_id>_checkpoint_0003.bin.zst` near t=0.3, with 160 Domains.
- The matching `grids` directory (a directory symlink is accepted).

Python 3.9+, NumPy, `zstd`, and `srun` must be available in the job environment.
The output directory must not already exist. With absolute paths:

```sh
sbatch "$kit/cluster_187.sbatch" \
  "$kit" \
  /absolute/path/wt-pre187/bin-opt/climate \
  /absolute/path/wt-stage187/bin-opt/climate \
  /absolute/path/wt-stage187/bin-check/climate \
  /absolute/path/cluster-fixture \
  /absolute/path/new-stage187-results \
  'RK4; pre187-baseline/stage187-candidate manifests; record actual compiler/MPI versions and flags here'
```

The supplied job requests 83 tasks in partition `bb`, one thread per task,
exclusive nodes and one hour. Adjust site-specific allocation settings if
needed. The Python driver runs once in the batch shell and launches sequential
83-task steps; do not launch 83 copies of the Python driver.

## What the script does

1. Checks J5 fixture inputs and seed header, snapshots executables/input/seed/grids,
   records their hashes, grid identity, build notes and Slurm allocation details.
2. Runs candidate checked oracle-off and oracle-on. Requires successful
   completion, checkpoint 4, its reload and subsequent stepping, and exact
   semantic checkpoint equality. It stops before timing if this gate fails.
3. Runs one excluded warmup of each optimized executable.
4. Runs three measured pairs, reversing alternate pairs: baseline/candidate,
   candidate/baseline, baseline/candidate. Every run starts from a fresh copy of
   checkpoint 3 and writes to its own directory. All optimized checkpoints must
   match the baseline warmup exactly. There is no checked-versus-optimized gate.
5. Reports all elapsed times, ordinary-step medians, group medians and the
   observed elapsed-time reduction in `results.json`. Full logs, per-run analysis
   and exact checkpoint comparison evidence remain in the output directory.

All inherited `WAVETRISK_*` settings are cleared before selecting the explicit
oracle/profile mode; thread counts are set to one. Detailed diagnostics and
profiling are off during measured runs. The first checked and optimized
references are self-compared only to establish finite readable checkpoints;
subsequent gates compare independent runs.

To append separate baseline/candidate detail profiles, run `cluster_187.py`
directly inside an existing allocation with the same arguments and `--profiles`.
The selected short window now ends at 0.315 and expects checkpoint 4 only.
The original 0.3225 window wrote checkpoints 4 and 5, as demonstrated by cluster
job 54447 (19 steps and both reloads). Its first harness failure was an incorrect
expected checkpoint list, before numerical comparison. The short checked-off
run later passed in 60.64 seconds; checked-on stopped on missing grid data,
so that attempt did not establish cluster oracle equality or timings.
Subsequent job 54452 passed all gates and measured a 4.71% elapsed-time
reduction; see [the full result](parallel-blocks-cluster-187-54452.md).

The runner now snapshots the grid directory and materializes all symlinks, then
makes and hashes a real grid copy for each run before timing starts. The required
J5 file `HR95JT_WT_004` must be present. It no longer depends on the original
grid directory after staging. The original fixture and output directory must
initially be accessible from the allocated nodes; on bbserv use shared home
storage rather than an uninitialized `/mnt/beeond` path.

For a deliberately different window, select its expected checkpoint list
explicitly and use identical inputs for every compared run. The default short
fixture does **not**
require a post-restart remap: it records whether one occurred and requires
continued stepping. Add `--require-post-restart-remap` when the selected window
is intended to provide that coverage. Changing `time_end` changes integer-clock
scaling; do not alter it between compared binaries.

Inspect placement, node load, memory/swap evidence and run-to-run spread before
claiming a speedup. The runner does not pretend to qualify distributed memory
pressure from the batch process alone. Slurm accounting is disabled on bbserv, so retain the batch and per-run logs;
`scontrol show job JOBID` can supply allocation state while Slurm retains it.
A checkpoint match checks semantic numerical fields, thresholds and topology;
compressed bytes and directory load weights are not required to match.

## Local verification

All 60 helper tests pass, including four new harness tests. These exercise a complete synthetic campaign with a
fake `srun`, exact rejection of a 1e-15 field change, failure before timing on a
zero-exit runtime error, environment cleanup, immutable seed handling and
restart-coverage checks. The synthetic campaign now moves the original grids
after the first run and verifies subsequent runs still read their own copies. No real cluster test has been run while bbserv is down.
The original legacy and checked-versus-optimized discrepancies remain separate
unresolved gates. No numerical tolerance is relaxed by this script.
