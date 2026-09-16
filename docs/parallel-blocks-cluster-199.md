# Stage 199 cluster diagnosis

Seven fresh J5/J7 builds and 27 short runs compare the unchanged Stage 197
numerical candidate, Stage 196 historical control, and pinned legacy Domain code.
The four additions split basic operators and adaptation construction, record
ordered producer/exchange batches, and count storage accesses without per-access
clocks. This campaign diagnoses causes; it does not test a new optimization.

Build and run on exclusive nk/bb13 with 83 ranks, contraction off and common physics
objects. Keep the existing fixture. All previous correctness, restart, workload,
placement and build-identity gates remain enabled. Observer overhead is checked
separately for timers, ordered records and the access census.

## Transfer from murier

```sh
scp /private/tmp/wavetrisk-stage199-cluster-kit.tar.gz kevlahan@bbserv:~/
```

## Prepare and build on bbserv

```sh
cd "$HOME"
tar -xzf wavetrisk-stage199-cluster-kit.tar.gz
kit="$HOME/wavetrisk-stage199-cluster-kit"
tree="$HOME/wt-stage199"
fixture="$HOME/stage187-fixture-20260911-145307"
state="$HOME/stage199-state.sh"
printf 'kit=%q\ntree=%q\nfixture=%q\n' "$kit" "$tree" "$fixture" > "$state"
git -C "$HOME/wav" bundle verify "$kit/stage196-baseline.bundle"
git -C "$HOME/wav" fetch "$kit/stage196-baseline.bundle" HEAD
git -C "$HOME/wav" worktree add --detach "$tree" 9bc3a19556dab79b3cba3944cf5df53517389da3
build_job=$(sbatch --parsable "$kit/build_199.sbatch" "$kit" "$tree")
printf 'build_job=%q\n' "$build_job" >> "$state"
tail -F "$HOME/wavetrisk-199-build-${build_job}.log"
```

The old revision is intentional: the kit first builds the historical control,
then applies the frozen numerical candidate and observer patches. Do not edit
Makefiles after the build. The C clock helper uses `mpicc`.

Wait for `PASS: seven fresh J5 builds ...`. A pending job may not have a log yet;
check `squeue -j "$build_job"` before submitting another build. Ctrl-C stops tail,
not the job.

## Submit tests with a fresh output directory

```sh
cd "$HOME"
state="$HOME/stage199-state.sh"
source "$state"
out="$HOME/stage199-nk-results-$(date +%Y%m%d-%H%M%S)"
job=$(sbatch --parsable "$kit/cluster_199.sbatch" "$kit" "$tree" "$fixture" "$out")
printf 'out=%q\njob=%q\n' "$out" "$job" >> "$state"
echo "Test job: $job"
echo "Results: $out"
tail -F "$HOME/wavetrisk-199-${job}.log"
```

To monitor or diagnose a job that disappears from the queue:

```sh
squeue -j "$job"
scontrol show job "$job"
tail -n 80 "$HOME/wavetrisk-199-${job}.log"
```

## Collect and transfer after completion or failure

```sh
cd "$HOME"
state="$HOME/stage199-state.sh"
source "$state"
python3 "$kit/collect_cluster_logs_199.py" --stage 199 --job "$job" --results "$out" --tree "$tree"
scp "$HOME/stage199-review-${job}.tar.gz" kevlahan@murier.math.mcmaster.ca:~/Desktop/
```

Return that single archive. It includes `results.json`, signed operation rankings,
`batches-summary.json`, each run's detailed batch summary and raw per-rank records,
logs, build identities and placement/memory evidence. Checkpoints and executables
are excluded. The collector refuses to overwrite an existing review archive.
