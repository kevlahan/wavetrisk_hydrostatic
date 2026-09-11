# Stage 188: bbserv arithmetic-control test

Use the existing successful J5/J7 fixture and 83 tasks. This compares three
builds of the same Stage 187 source with the common checkpoint I/O fix:

| Runner role | Solver build | FP_CONTRACT |
| --- | --- | --- |
| baseline | Optimized Stage 187 | default |
| candidate | Optimized Stage 187 | off |
| checked | Checked Stage 187 | off |

The baseline here is Stage 187 itself, not the pre-187 source. Checked-off and
checked-on mean oracle validation disabled/enabled; both have contraction off.

The build helper verifies source hashes against the tested Stage 187 sources,
applies only the Makefile option patch, generates missing physics dependency
includes, and builds in three fresh directories. It verifies identical external
physics object hashes across the three builds and records commands and binary
hashes. Existing bin-opt/bin-check executables are preserved. It refuses unknown
source changes or existing Stage 188 build directories. On a build failure,
inspect its printed build log before retrying; do not reuse partial builds.

## Transfer and build

The prepared archive is `/private/tmp/wavetrisk-stage188-cluster-kit.tar.gz` on
the Mac. Transfer it using the usual bbserv SSH connection:

```sh
scp /private/tmp/wavetrisk-stage188-cluster-kit.tar.gz bbserv:~/
```

On bbserv, load the same GNU Fortran/MPI/NetCDF/Python environment that produced
job 54452. Then:

```sh
cd "$HOME"
tar -xzf wavetrisk-stage188-cluster-kit.tar.gz
kit="$HOME/wavetrisk-stage188-cluster-kit"
tree="$HOME/wt-stage187"
fixture="$HOME/stage187-fixture-20260911-145307"

python3 "$kit/build_arithmetic_188.py" "$kit" "$tree"
```

This uses the already tested Stage 187 worktree. Alternatively, a clean checkout
of `b270c0f7` is accepted and receives the same Makefile-only patch. No numerical
source patch is applied. The helper creates:

* `bin-188-default/climate`
* `bin-188-off/climate`
* `bin-188-check-off/climate`
* `arithmetic-188-builds.json`

Do not continue unless the helper ends with `PASS: three fresh builds`.

## Submit

Keep the home-directory fixture that passed job 54452: max_level=7,
time_end=0.315, dt_write=0.01, resume=3. It contains the checkpoint-3 seed and
real grid assets, including `grids/HR95JT_WT_004`. The seed has 160 Domains.
The runner snapshots these inputs before starting.

```sh
grep -E '^(max_level|time_end|dt_write|resume)' "$fixture/simple.in"
out="$HOME/stage188-results-$(date +%Y%m%d-%H%M%S)"
sbatch "$kit/cluster_188.sbatch" "$kit" "$tree" "$fixture" "$out"
```

Submit from the home directory so the batch log appears there. Each submission
requires a new output directory. Monitor the job ID returned by sbatch:

```sh
tail -f "$HOME/wavetrisk-188-JOBID.log"
```

## Acceptance and interpretation

The runner requires oracle-on/off equality for the checked build and exact
optimized-off/checked-off checkpoint equality. All runs must reload checkpoint
4 and perform a subsequent remap. Default optimized repeats must reproduce their
own warmup checkpoint exactly. It excludes correctness runs and one warmup per
optimized build from timings, then runs three alternating measured pairs.

Default versus contraction-off checkpoint equality is recorded separately as
`default_vs_off_diagnostic`; either true or false is allowed because the builds
use different arithmetic settings. A difference between optimized-off and
checked-off is a mandatory gate failure and stops the campaign before timing.
No tolerance is loosened to accept such a failure. There is no new legacy test
in this campaign; its purpose is cross-mode arithmetic validation and the cost
of the setting on bbserv.

Success means all ten runs report PASS and results.json has status `passed`.
The summary includes `off_elapsed_increase_percent` and
`off_ordinary_step_increase_percent`: positive means contraction off is slower,
negative means faster. These are observed timings, with the same placement/load
and distributed-memory limitations as the previous campaign. The new setting
does not become a production default automatically.

Return the new output directory's `results.json` and the complete
`wavetrisk-188-JOBID.log`. If a run fails, also return that run's `run.log`.

## Completed campaign and default

Job 54465 passed the full campaign. See the [validated result](parallel-blocks-cluster-arithmetic-188-54465.md).
The Makefile now defaults to contraction off. The build helper continues
to request `FP_CONTRACT=default` and `FP_CONTRACT=off` explicitly, preserving the
comparison described here. The transferred kit remains a frozen reproduction
of the original campaign.
