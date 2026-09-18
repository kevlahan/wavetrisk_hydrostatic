# Stage215: transfer, build, run and retrieve

Candidate: compile the remaining generic Domain halo integer routes once per
topology/level interval and pack directly. Baseline is accepted Stage202,
Stage203 off. Local validation is included in `LOCAL_VALIDATION.json`.
The100ms target is unproven; the runner stops early for small/negative gains.
No profiled timing campaign or automatic promotion is included. Main is untouched.

## 1. Transfer from murier to bbserv

On murier:

```bash
scp "$HOME/Desktop/wavetrisk-stage215-cluster-kit.tar.gz" kevlahan@bbserv.math.mcmaster.ca:~/
```

## 2. Extract and build on bbserv

Use the same compiler/MPI/Python/numpy/zstd environment as the previous campaigns.
Run this block on bbserv. It uses a subshell, so failure cannot log you out.

```bash
(
set -e
cd "$HOME"
tar -xzf wavetrisk-stage215-cluster-kit.tar.gz
bash "$HOME/wavetrisk-stage215-cluster-kit/bootstrap_215.sh"
)
source "$HOME/stage215-state.sh"
tail -F "$HOME/wavetrisk-215-build-${build_job}.log"
```

Wait for `PASS: three fresh J5 builds...`, then press **Ctrl-C** to stop tail.
If the build fails, send the build log; do not launch the campaign. Bootstrap
extracts a complete source archive, verifies it and refuses an existing
`$HOME/wt-stage215`. It preserves any previous state file. It never uses an
incremental Git bundle or checks out a branch.

## 3. Run on bbserv

```bash
bash "$HOME/wavetrisk-stage215-cluster-kit/launch_215.sh"
source "$HOME/stage215-state.sh"
tail -F "$HOME/wavetrisk-215-${job}.log"
```

The runner validates at16 tasks, screens one16-task A/B pair, checks a promising
candidate at83, and confirms a sufficiently large gain with fresh A–B–B–A.
It uses `$HOME/stage187-fixture-20260911-145307`. The source, executable, input and
placement checks are automatic. A no-gain or below-target decision is a completed
experiment, not a reason to rerun. The checked validation run is deliberately
slower because it audits every route use.

Wait for `COMPLETE: ...`, then press **Ctrl-C**. `squeue -u "$USER"` shows active
jobs; Slurm accounting is not required.

## 4. Create, verify and transfer the archive FROM bbserv TO murier

After the job ends, run **all of this block on bbserv**:

```bash
(
set -e
source "$HOME/stage215-state.sh"
active=$(squeue -h -u "$USER" -o '%A')
if printf '%s\n' "$active" | awk -v wanted="$job" '$1==wanted {found=1} END {exit !found}'; then
    echo "Job $job is still queued/running; collect after it ends."
    exit 1
fi
archive="$HOME/stage215-review-${job}.tar.gz"
if [ ! -f "$archive" ]; then
    python3 "$kit/collect_cluster_logs_215.py" \
        --stage 215 --job "$job" --results "$out" --tree "$tree"
fi
python3 "$kit/collect_cluster_logs_215.py" --verify "$archive"
sha256sum "$archive"
scp "$archive" kevlahan@murier.math.mcmaster.ca:~/Desktop/
)
```

The identical operation is also available as:

```bash
bash "$HOME/wavetrisk-stage215-cluster-kit/collect_215.sh"
```

Attach `stage215-review-JOBID.tar.gz` from murier's Desktop. If scp is interrupted,
repeat step4: the existing archive is verified and retransferred, not recreated.
A job removed from Slurm's active records does not cause the earlier
`Invalid job id specified` error. A genuine scheduler-query failure stops safely.
If archive verification fails, preserve it for inspection.
