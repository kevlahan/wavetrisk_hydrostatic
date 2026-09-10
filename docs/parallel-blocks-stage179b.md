# Stage 179b — complete the profiling evidence, not another numerical change

Builds on the validated, uncommitted Stage 179a instrumentation. No commit is
created by this stage. Production arithmetic, communication ordering, routing
and numerical tolerances are unchanged. `main` is not edited or checked out.

## Standalone srun workflow (launcher update)

If the usual standalone `srun` works but batch/interactive allocations do not,
use `--independent-srun --partition=bb`. **No recompilation is needed.** Update
`test/parallel_block_profile/experiment.py` and `compare.py` on the cluster.

From a normal shell, not under `srun`, `sbatch` or `salloc`:

```sh
python3 ~/wav/test/parallel_block_profile/experiment.py \
  --independent-srun --partition=bb \
  --legacy ~/wav-legacy-179b/bin/climate \
  --legacy-build-log ~/wav-legacy-179b/build.log \
  --legacy-identity ~/wav-legacy-179b/legacy-build.json \
  --block ~/wav/bin/stage179b-opt/climate \
  --block-build-log ~/wav/build-stage179b-opt.log \
  --fixture ~/hydro/climate/profile-seed \
  --out "$PWD/profile179b-srun-01"
```

For each run the driver copies the appropriate executable into its fresh run
directory as `climateJ5`, verifies its hash, and executes exactly:

```sh
srun -n 83 --partition=bb ./climateJ5 simple.in
```

Combined stdout/stderr is displayed live and saved, equivalent to `2>&1 | tee`,
but preserving `srun`'s failure status. Logs have `RK4_production.log` or
`RK4_profile.log` aliases alongside `run.log`. Each subprocess finishes before
the next starts. There is no encompassing allocation, no Python/GNU-time rank
wrapper and no added CPU-binding options in this mode. All-oracle-off timing
and application-detail switches are still selected automatically. `--partition`
can be changed or omitted to match a working srun command; optional `--nodelist`
is forwarded only when explicitly supplied. No resource restrictions are bypassed.

Numerical/input/hash checks and per-step application profiles are retained.
Node/core identity and per-rank RSS are **not measured** in this direct pattern.
The comparator therefore reports `placement_verified: false` and
`valid_for_comparison: false`, while still returning usable, **provisional**
timing ratios if the numerical/run checks pass (`measurements_usable: true`).
Do not present these as a verified same-hardware speedup. The original
shared-allocation mode below retains strict placement verification. Optional
`--sample` explicitly adds the external perf shell wrapper to either mode;
omit it for the exact direct launch shown above.

The seed must contain inputs only, not `climateJ5`, prior run logs or profiling
output directories. Choose a new `--out` directory for each experiment. Avoid
overlap with a queued batch experiment that may start later.

This launcher update has local command-construction, exit-status/tee,
placement-policy and parser tests. No Fortran source or executable changed.

## Questions these measurements can answer

1. How much slower is this exact block binary than the unchanged legacy binary
   for matching ordinary timesteps on the same ranks/cores?
2. Which previously unclassified dynamics costs are RK construction, stage
   retention/publication, native tendency assembly, state refresh, wavelet
   driving, grid-change synchronization or scalar/velocity preparation?
3. Which callers account for boundary work, and how much is repeated route
   scanning/packing, receive planning, allocation/clearing, local copying,
   MPI posting, waiting or installation?
4. Which ranks do the work versus spend time in MPI? Are aggregate window
   impressions also true on each ordinary step, excluding restart/remap?
5. Do external samples in *both unmodified numerical executables* implicate
   allocation/copy/lookup dispatch, native kernels or MPI progress? Do extra
   instructions or reduced IPC corroborate a computational/layout hypothesis?

We cannot label all MPI time as network cost, all extra instructions as block
overhead, or the whole measured legacy gap as irreducible overhead. The legacy
executable deliberately has **no new source-level timers**. Its external samples
provide call-path evidence, not a matching set of exclusive phase timers.
No conclusion that the complete gap has been attributed is justified while
significant block self time or sampled helper callers remain unexplained.

## Added instrumentation

Existing region IDs 1–42 are unchanged. New IDs are:

| IDs | Meaning |
|---|---|
| 43–47 | Native wavelets, compression, RK assembly, retention/publication, tendency assembly residual |
| 48–53 | Grid-change synchronization, trend/prognostic refresh, scalar capture, velocity preparation, RK compatibility |
| 54–57 | Split physics, remap, restart, oracle work |
| 58–63 | Boundary scan/pack, receive-route scan, MPI posting, scan/install, local copies, allocation/zeroing |
| 64–65 | Boundary caller outside a timed scope; wavelet driver residual |

Boundary caller attribution includes its children (including waits) and is a
**cross-cut**, not additional self time. The nearest timed caller is reported;
start and finish can have different callers. Counts are entered start/finish
routines, not distinct messages. Already-up-to-date early returns are excluded.

Counters 9–16 record boundary send/receive doubles, posted requests, candidates
examined by pack/route/install scans, actual buffer capacity zeroed, and fields
per entered start. Scan candidates are counted from loop lengths outside the
node loops; no per-node clocks or instrumentation collectives are added.
These counters quantify traversal amplification and retained buffer clearing,
not hardware memory traffic. Field counts are field slots, not unique columns.

Per-step rank-local self-wall/self-CPU vectors and total wall/CPU are buffered
until the existing report. Sequence IDs persist across report resets; restart
and remap entries are recorded. There is no per-step reporting collective or
barrier. Up to 256 steps per window are retained. Overflow is printed explicitly
and the comparison tool rejects incomplete detail attribution. Use shorter
report windows for longer experiments. Additional gathers/broadcast occur only
at reporting. Added detail overhead still must be measured empirically.

Per-step sums must equal each rank's timestep total within printed precision.
Window self sums can also include instrumented work between timesteps. Neither
independent phase maxima nor inclusive rows may be summed into a critical path.

## Build and tests

Python helpers require Python 3.9 or newer. Build both executables on the same
partition/CPU family used for execution. Both Makefiles default to
`-march=native`; do not move an EPYC-specific binary to an incompatible chip.
The Slurm harness also requires GNU `/usr/bin/time` on the compute nodes.

First run the usual RK3 and RK4 checked oracle tests, with the integrator chosen
at build time as before. All five switches for each debug run:

```sh
export WAVETRISK_VALIDATE_BLOCK_DYNAMICS=1
export WAVETRISK_VALIDATE_BLOCK_ADAPTATION=1
export WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS=1
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=1
export WAVETRISK_PROFILE_BLOCK_DETAIL=1
```

Production timing uses all five switches set to **0**. Detailed optimized
profiling uses the first three **0** and the last two **1**. There are no new
profiling environment switches. The experiment driver sets these itself and
refuses unclassified inherited `WAVETRISK_*` switches.

Local instrumentation/unit checks:

```sh
sh test/parallel_block_profile/run.sh
```

The runner enforces 132-column source width and checks both `-O0` and `-O2`
with bounds and floating-point traps. New tests cover boundary caller
conservation, step snapshots/reset/overflow, parsing, alternating run order,
placement rejection, checkpoint exclusion and paired comparison arithmetic.

## Compile authoritative legacy main without editing it

The requested main tip is `b82647b36eb825aed0c0aae4e92adcf4ab46d608`
(`Minor tidying.`). The helper archives exactly this revision into a **new
directory outside the repository**, compiles its own Makefile, and checks all
319 tracked file contents against Git afterward. It never switches branches,
commits, modifies a tracked legacy file, or adds legacy instrumentation.

The archive lacks the generated physics `Makefile.inc`; the helper creates
that untracked dependency file using the package's supplied generation scripts.
It does not call the physics `clean/nice` target, which formats source files.

On the cluster, choose a new build directory:

```sh
python3 ~/wav/test/parallel_block_profile/build_legacy.py \
  --repo ~/wav --out ~/wav-legacy-179b
```

Output: `~/wav-legacy-179b/bin/climate`, `build.log`, `legacy-build.json`.
The identity records source revision, unchanged-source verification, executable
SHA-256 and build host. A build failure is reported without patching main.

Build the optimized block executable separately and retain its complete log:

```sh
cd ~/wav
set -o pipefail
make -j1 DEBUG=false PARAM=param_J5 TEST_CASE=climate \
  BUILD_DIR=build/stage179b-opt BIN_DIR=bin/stage179b-opt \
  2>&1 | tee build-stage179b-opt.log
```

Ensure the block build is **RK4**, matching main's default RK4. Do not run this
comparison with an RK3 block binary. Do not overwrite `climateJ5` mid-experiment.

## Matched experiment — one existing Slurm allocation

Create a small, frozen seed directory, e.g. `~/hydro/climate/profile-seed`,
containing the current `simple.in`, the original checkpoint 3, and a `grids`
link. Use the checkpoint filename implied by `run_id` in that input. Do not
point the driver at the whole live test/output directory. Directory-linked
assets are identified by resolved path, not recursively content-hashed: keep
grids immutable. Regular files and file-linked checkpoints are content-hashed.
Every execution gets a fresh copy; file-linked inputs are materialized so a
checkpoint write cannot overwrite the seed. Results directories must be new.

Run inside the same 83-task allocation, with both build logs available:

```sh
python3 ~/wav/test/parallel_block_profile/experiment.py \
  --legacy ~/wav-legacy-179b/bin/climate \
  --legacy-build-log ~/wav-legacy-179b/build.log \
  --legacy-identity ~/wav-legacy-179b/legacy-build.json \
  --block ~/wav/bin/stage179b-opt/climate \
  --block-build-log ~/wav/build-stage179b-opt.log \
  --fixture ~/hydro/climate/profile-seed \
  --out "$PWD/profile179b-timing-01"
```

Default: six alternating-order legacy/block pairs (12 runs), followed by two
alternating off/detail pairs (4 runs). `srun` uses 83 tasks, core binding and
kill-on-bad-exit, with one CPU per task and OpenMP/BLAS thread limits set to one
for both binaries. These settings and selected MPI progress variables are
recorded. It does not allocate nodes or submit a job. Each run records
per-rank hostname, CPU model, CPU affinity and allowed memory nodes before exec.
GNU time records each rank's process CPU, peak RSS, faults and context switches
without adding numerical source instrumentation. These are whole-run resource
measurements; summed rank RSS peaks are not simultaneous resident memory.
The driver stops on failure or printed-state disagreement, rather than timing
a different workload. It verifies the legacy identity and checks executable
hashes and frozen regular inputs before each run. Do not change shared libraries
or grid assets during the experiment. Check the saved compiler/link commands
for matching compiler/MPI/physics libraries and appropriate ISA flags.

```sh
python3 ~/wav/test/parallel_block_profile/compare.py profile179b-timing-01 \
  > profile179b-comparison.json
```

The comparison rejects changed placement, incomplete runs or different printed
workloads. It excludes checkpoint-containing timesteps and treats a **run**,
not each timestep, as the replication unit. It reports paired ratios, their
range and a descriptive paired-run bootstrap interval (at least three pairs).
This is not protection against systematic bias or changing cluster contention;
increase repetitions if intervals are broad or run order matters. Detail
on/off pairs quantify enabled instrumentation perturbation, not the residual
cost of disabled calls compared with a pre-instrumentation executable.

## External CPU evidence from both executables

Repeat the same `experiment.py` command with a **new** `--out` directory and
`--sample record`. This runs one legacy and one block executable through
`perf`, default selected ranks `0,30,52,82`; override with `--sample-ranks` after
reading the new profile. All application oracles/profilers are off. Use another
new directory with `--sample stat` for instructions/cycles/cache counters.
Sampling runs are separate from timing observations. No code is injected into
the legacy source, no MPI interposition is used, and cluster security settings
are not changed. Missing/denied perf access is an error, not a silent fallback.
Unsupported/multiplexed counters are flagged explicitly; retain those warnings
and any errors if the cluster disallows sampling.

```sh
python3 ~/wav/test/parallel_block_profile/samples.py profile179b-record-01 \
  > profile179b-sample-summary.json
python3 ~/wav/test/parallel_block_profile/samples.py profile179b-stat-01 \
  > profile179b-counter-summary.json
```

The record command also creates per-rank self-symbol and caller/callee reports
beside the raw recordings; it refuses to overwrite existing reports. Keep the
exact executables available for symbol resolution. Counters report IPC and
multiplexing/unavailable warnings, not invented DRAM bandwidth. Whole-run
samples include startup, output, compression children and MPI spinning; do not
use their percentages as ordinary-step shares. Attribute hot helpers by call
path and use the application timers to bound their relevant phase costs.

If startup/output dominates samples, collect a separate, identical *input-only*
longer benchmark for both binaries with output deferred beyond its time_end.
Keep the normal restart test as the correctness check; do not alter main source
to isolate timing. Record the changed input hash and compare it only with its
matching counterpart.

## Turning evidence into an optimization decision

The primary result is measured excess `block_time - legacy_time` for matching
ordinary steps, plus the fraction of block runtime that must disappear to
reach `1.2 * legacy_time`. For a measured ratio 2.5, this fraction is 52%.

Optional explicit what-if calculation (remove the three serialization proofs):

```sh
python3 ~/wav/test/parallel_block_profile/compare.py profile179b-timing-01 \
  --remove 6:1 --remove 7:1 --remove 9:1
```

This removes only the specified fraction of **exclusive rank-average self
cost**, assumes no replacement cost and leaves MPI waits unchanged. It is not
a predicted job speedup. Do not include both an inclusive phase and its nested
children, add independent slow-rank maxima, or count boundary caller totals
again. No automated removal scenario labels the retained time irreducible.

The next implementation should be selected using this table of evidence:

| Observed evidence | Interpretation to test |
|---|---|
| Large boundary scan/local-copy time and repeat candidates per transferred value | Repeated horizontal traversal/packing; batch columns and reuse schedules |
| Large RK retention/grid-change/self-copy costs with allocation/copy samples | Repeated representation transfer or snapshot overhead; keep data resident |
| Native replay arithmetic dominates, with excess instructions or low IPC | Block kernel/layout/dispatch issue; inspect same physics call paths in legacy |
| Wait-heavy ranks lag compute-heavy owners | Dependency/imbalance or too many phases; not automatically a bandwidth problem |
| Small proof/geometry-expansion self time | Cleanup may help, but cannot explain a multi-fold gap by itself |

Return `experiment.json`, the comparison JSON, all run logs and rank-placement
JSONs/resource reports, both build logs/legacy identity, and the sample reports/counter summary
(or the perf access error). With these we can report measured/likely removable
costs, explicit speedup scenarios and the *unexplained remainder*, separately
from any genuinely necessary block bookkeeping. An inherent overhead floor
requires evidence of what must remain after a proposed design, not just today's
timing difference.

## Local verification completed

- GNU Fortran 15.2: checked RK3/RK4 builds and optimized RK4 build completed.
  The checked builds enforce `-ffree-line-length-132`, `-Werror`, bounds checks
  and floating-point traps. The optimized compiler retains an existing
  `air_temperature` maybe-uninitialized warning in `parallel_block.f90`;
  this stage does not modify that numerical code.
- Four-rank RK3 and RK4 one-step oracle tests passed, including all RK stages
  and the final inverse/adaptation. Printed states match Stage 179a.
- Four-rank optimized ten-step detail run passed remaps and checkpoint/restart.
  All ten per-step rank sets are complete; sequence numbers survive reset;
  no records were dropped. Self sums reconcile within printed precision.
- That run's printed states match the completed Stage 179a serial run.
  Checkpoint 4 is byte-identical, SHA-256:
  `6086e1fab875369095b1dea9800af7a39261ab751d6f10916718ed3c081179f3`.
- Four-rank optimized one-step test with all oracles/profilers off passed.
  Its printed state matches an unchanged-main four-rank smoke test.
- Main `b82647b3` compiled and all 319 tracked contents were verified unchanged
  after building. Its one-step restart smoke test passed.
- Instrumentation tests pass at `-O0`/`-O2`; seven Python unit tests pass.
  Strict source-width checks and `git diff --check` pass.

MPI tests ran serially in isolated temporary fixtures; the standard test
directory and `bin/climateJ5` were not modified. Local durations are not an
83-rank performance comparison. Slurm placement collection and Linux perf
execution cannot be exercised on this Mac: the supplied harness has local
syntax/unit coverage, but its cluster integration must still be verified.
