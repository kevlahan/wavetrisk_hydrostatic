# Stage223: transfer, build, run and retrieve

Reserve one exclusive nk node (two 128-core sockets, 256 physical cores).
The campaign uses both sockets and runs the 40/41-rank J4 gates before the
160/167-rank J5 gates. All runs include checkpoint/restart correctness checks;
this stage does not measure or accept speedups. Stop on the first failed gate.

Use the usual bbserv compiler/MPI/Python environment used for the previous
campaigns. Required commands include `mpif90`, `make`, Python with NumPy,
`zstd`, `srun`, `numactl`, `lscpu` and the existing NetCDF/LAPACK dependencies.
The kit includes its J4 fixture. It expects the existing J5 fixture at
`$HOME/stage187-fixture-20260911-145307`.

## 1. Transfer from murier

Run on **murier**:

```bash
scp "$HOME/Desktop/wavetrisk-stage223-cluster-kit.tar.gz" kevlahan@bbserv.math.mcmaster.ca:~/
```

## 2. Extract and submit the build on bbserv

Run on **bbserv**. These are scripts, not commands to source into the login
shell; failures will not log you out. No Git clone/bundle is required.

```bash
tar -xzf "$HOME/wavetrisk-stage223-cluster-kit.tar.gz" -C "$HOME" && bash "$HOME/wavetrisk-stage223-cluster-kit/bootstrap_223.sh"
```

Bootstrap verifies and installs the source snapshot into a fresh
`$HOME/wt-stage223/reference`, writes `$HOME/stage223-state.sh` and submits the
build. It refuses an existing source tree and preserves any previous state.

Watch the build:

```bash
source "$HOME/stage223-state.sh"
tail -F "$HOME/wavetrisk-223-build-${build_job}.log"
```

Wait for `PASS: six fresh J4/J5 builds...`. Press **Ctrl-C** to stop following
the log; this does not cancel the job. If the build fails, send the build log
and the reported `build-223-*.log`; do not submit the campaign.

## 3. Submit the campaign

Run on **bbserv**, after the build passes:

```bash
bash "$HOME/wavetrisk-stage223-cluster-kit/launch_223.sh"
source "$HOME/stage223-state.sh"
tail -F "$HOME/wavetrisk-223-${job}.log"
```

The launcher verifies all six builds and submits onto the same nk node used for
compilation. It records the actual job ID and output directory in the state
file. It refuses duplicate submission. The batch reserves all 256 physical
cores, while individual runs use 40, 41, 160 or 167 pinned ranks.

Successful completion ends with `PASS: 41/40 and 167/160 capability gates...`.
A failure stops the remaining run matrix. Retrieve the review archive for
either outcome; no `sacct` command is required.

## 4. Create, verify and transfer the review archive FROM bbserv TO murier

After the campaign ends, run this **on bbserv**:

```bash
bash "$HOME/wavetrisk-stage223-cluster-kit/collect_223.sh"
```

That single command checks the user's active queue, creates the archive only
if absent, verifies every archived file hash, prints its checksum and executes
`scp` to `kevlahan@murier.math.mcmaster.ca:~/Desktop/`. It avoids querying a
purged job ID and works with disabled Slurm accounting. If transfer is
interrupted, rerun the same command; the verified archive is reused.

For clarity, the archive/verification/transfer portion is:

```bash
(
set -e
source "$HOME/stage223-state.sh"
archive="$HOME/stage223-review-${job}.tar.gz"
if [ ! -f "$archive" ]; then
    python3 "$kit/collect_cluster_logs_223.py" --stage 223 --job "$job" --results "$out" --tree "$tree"
fi
python3 "$kit/collect_cluster_logs_223.py" --verify "$archive"
sha256sum "$archive"
scp "$archive" kevlahan@murier.math.mcmaster.ca:~/Desktop/
)
```

Prefer `collect_223.sh`, which additionally refuses collection while the job
is still queued or running. The returned file is
`~/Desktop/stage223-review-JOBID.tar.gz` on murier. Attach it for review.

## Placement repair after job 54724

The original kit requested `--cpus-per-task=1` for each srun step. Slurm documents
that an explicit request implies `--exact`, restricting the step allocation and
potentially making planned cores on the other socket unavailable. The original
probe discarded stderr, so the exact remote error is unconfirmed. Job 54724
stopped at the first 40-rank probe before any solver execution. The initial
`tail -F` missing-file message is normal while the job starts.

The repair removes the explicit step CPU request and batch CPU/thread hints,
records and clears inherited step restrictions, and keeps exclusive whole-node
allocation with manual pinning to one hardware thread per physical core. It
saves probe stdout/stderr and return codes even on failure. Sources, the source
manifest, candidate patch and six built binaries remain unchanged.

Transfer the small repair from murier:

```bash
scp "$HOME/Desktop/wavetrisk-stage223-launcher-repair.tar.gz" kevlahan@bbserv.math.mcmaster.ca:~/
```

On bbserv, press Ctrl-C to stop the old tail, then:

```bash
tar -xzf "$HOME/wavetrisk-stage223-launcher-repair.tar.gz" -C "$HOME" && bash "$HOME/wavetrisk-stage223-launcher-repair/retry_223.sh"
source "$HOME/stage223-state.sh"
tail -F "$HOME/wavetrisk-223-${job}.log"
```

The script refuses an active job 54724, checks its saved failure preceded any
solver run, backs up the original helpers and state, verifies repaired helper
hashes, and resubmits using the existing six verified builds on their build
host. It preserves the failed results/log and refuses duplicate submission if
a different campaign job is already recorded. Do not rerun bootstrap or extract
the original kit over the repaired kit. Retrieve the new campaign with
`bash "$HOME/wavetrisk-stage223-cluster-kit/collect_223.sh"` after it ends;
this creates/verifies the archive and transfers it from bbserv to murier.

Validation: 13 Stage223 helper tests pass, including failed-probe diagnostics;
the isolated repair workflow passes active-job, identity, backup, retry,
idempotence, duplicate submission and archive/transfer checks. This local
workflow does not reproduce Slurm cgroups; the repaired cluster probes remain
the actual placement gate.

Reference: <https://slurm.schedmd.com/srun.html#OPT_cpus-per-task>.
