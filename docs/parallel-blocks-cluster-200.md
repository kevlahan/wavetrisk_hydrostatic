# Stage 200: controlled task-count and socket-placement tests

This kit compares unchanged Stage 197 block numerics with legacy at 21, 42 and
83 tasks on one socket, plus 83 tasks across two sockets. Five builds and 42
short solver runs retain correctness/restart gates while replacing the large
Stage 199 trace campaign with small per-rank phase/MPI/work records. There are
only two bare timings per implementation/configuration. Build and run on nk.
Do not edit `simple.in`, Makefiles or CPU maps between runs.

## Transfer from murier

```sh
scp /private/tmp/wavetrisk-stage200-cluster-kit.tar.gz kevlahan@bbserv:~/
```

## Prepare and build on bbserv

```sh
cd "$HOME"
tar -xzf wavetrisk-stage200-cluster-kit.tar.gz
kit="$HOME/wavetrisk-stage200-cluster-kit"
tree="$HOME/wt-stage200"
fixture="$HOME/stage187-fixture-20260911-145307"
state="$HOME/stage200-state.sh"
printf 'kit=%q\ntree=%q\nfixture=%q\n' "$kit" "$tree" "$fixture" > "$state"
git -C "$HOME/wav" bundle verify "$kit/stage196-baseline.bundle"
git -C "$HOME/wav" fetch "$kit/stage196-baseline.bundle" HEAD
git -C "$HOME/wav" worktree add --detach "$tree" 9bc3a19556dab79b3cba3944cf5df53517389da3
build_job=$(sbatch --parsable "$kit/build_200.sbatch" "$kit" "$tree")
printf 'build_job=%q\n' "$build_job" >> "$state"
tail -F "$HOME/wavetrisk-200-build-${build_job}.log"
```

The starting commit is intentional. The kit applies the frozen validated
candidate before building. `-march=native` is compiled on the same bb13 hardware
used for testing; contraction is off for both implementations. Wait for
`PASS: five fresh J5 builds ...`. Ctrl-C stops `tail`, not the job.

If the log has not appeared, check the existing job before submitting again:

```sh
squeue -j "$build_job"
scontrol show job "$build_job"
```

## Run with a fresh output directory

```sh
cd "$HOME"
source "$HOME/stage200-state.sh"
out="$HOME/stage200-nk-results-$(date +%Y%m%d-%H%M%S)"
job=$(sbatch --parsable "$kit/cluster_200.sbatch" "$kit" "$tree" "$fixture" "$out")
printf 'out=%q\njob=%q\n' "$out" "$job" >> "$state"
echo "Test job: $job"
echo "Results: $out"
tail -F "$HOME/wavetrisk-200-${job}.log"
```

The batch allocation requests the full 256-core node exclusively, with all node
memory available. **Do not change its `--ntasks=256` to 21, 42 or 83.** The
Python driver launches the four intended rank configurations inside that one
allocation and verifies the discovered core maps. Solvers run sequentially.

`numactl`, `lscpu`, `srun`, `zstd`, Python/numpy and the build toolchain must be
available. Placement is checked before the first solver starts. The fixture
must retain `time_end=0.315`, `dt_write=0.01`, `resume=3`, J5/J7 and 30 layers.
The driver copies actual grid files, so run-directory symlinks do not depend on
node-local paths.

Expected progress lines name the configuration, implementation and mode, such
as `Starting n21-s1-candidate-bare-1`. Completion is:
`PASS: Stage 200 scaling campaign; results.json, scaling.md and scaling.csv ready for review.`
A passed campaign means all execution, placement and correctness checks passed;
it does not automatically identify an optimization or guarantee useful timing.

For failure diagnostics:

```sh
squeue -j "$job"
scontrol show job "$job"
tail -n 80 "$HOME/wavetrisk-200-${job}.log"
python3 - "$out/results.json" <<'PY'
import json, pathlib, sys
p = pathlib.Path(sys.argv[1])
if p.exists():
    r = json.loads(p.read_text())
    for key in ('status', 'active', 'error'):
        print(f'{key}: {r.get(key)}')
else:
    print('Setup failed before results.json; inspect the Slurm log.')
PY
```

## Collect the single review archive

After the job has stopped, whether successful or failed:

```sh
cd "$HOME"
state="$HOME/stage200-state.sh"
source "$state"
python3 "$kit/collect_cluster_logs_200.py" --stage 200 --job "$job" --results "$out" --tree "$tree"
scp "$HOME/stage200-review-${job}.tar.gz" kevlahan@murier.math.mcmaster.ca:~/Desktop/
```

Return that tarball. It includes the scaling tables, results and all gates,
per-rank raw phase/MPI/work records, observed CPU maps and topology, resource
snapshots, run/Slurm logs and build provenance. It excludes checkpoints and
executables. The collector refuses to overwrite an existing archive.

The decision criteria and measurement limitations are in `EXPERIMENT.md`.
