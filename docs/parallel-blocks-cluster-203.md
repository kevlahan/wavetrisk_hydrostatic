# Stage203 bb13 cluster protocol

The kit is self-contained and reconstructs validated Stage202-on and Stage203
from the supplied Stage196 bundle plus frozen patches. It also includes pinned
legacy sources and the common checkpoint I/O correction. No modification of
`~/wav` or an older test checkout is required. Build and run on bb13/nk because
both optimized builds use `-march=native`. Use fresh directories throughout.

## Transfer from murier

```bash
scp /private/tmp/wavetrisk-stage203-cluster-kit.tar.gz kevlahan@bbserv:~/
```

## Prepare and build on bbserv

```bash
cd "$HOME"
tar -xzf wavetrisk-stage203-cluster-kit.tar.gz
kit="$HOME/wavetrisk-stage203-cluster-kit"
tree="$HOME/wt-stage203"
fixture="$HOME/stage187-fixture-20260911-145307"
state="$HOME/stage203-state.sh"
printf 'kit=%q\ntree=%q\nfixture=%q\n' "$kit" "$tree" "$fixture" > "$state"
git -C "$HOME/wav" bundle verify "$kit/stage196-baseline.bundle"
git -C "$HOME/wav" fetch "$kit/stage196-baseline.bundle" HEAD
git -C "$HOME/wav" worktree add --detach "$tree" 9bc3a19556dab79b3cba3944cf5df53517389da3
build_job=$(sbatch --parsable "$kit/build_203.sbatch" "$kit" "$tree")
printf 'build_job=%q\n' "$build_job" >> "$state"
echo "Build job: $build_job"
tail -F "$HOME/wavetrisk-203-build-${build_job}.log"
```

Wait for `PASS: seven fresh J5 builds`. A missing log immediately after submission
is normal while queued. Check `squeue -j "$build_job"` rather than submitting a
duplicate. Ctrl-C stops tail, not the job. Do not edit the Makefile after building.
Build failures have individual `build-203-*.log` files in the source checkout.

## Run after the build passes

```bash
cd "$HOME"
source "$HOME/stage203-state.sh"
out="$HOME/stage203-nk-results-$(date +%Y%m%d-%H%M%S)"
job=$(sbatch --parsable "$kit/cluster_203.sbatch" "$kit" "$tree" "$fixture" "$out")
printf 'out=%q\njob=%q\n' "$out" "$job" >> "$state"
echo "Test job: $job"
echo "Results: $out"
tail -F "$HOME/wavetrisk-203-${job}.log"
```

The fixture is still J5/J7, resume 3, `time_end=0.315`, `dt_write=0.01`, with
checkpoint 4 reloaded and a later remap. The harness snapshots input and grid
files and verifies placement before solver runs. There are 18 runs. Checked-on
is deliberately slower. The final PASS means correctness completed; the separate
`prototype-summary.json` determines whether the substantial-gain target was met.

Do not time checked binaries or observers as production performance. Accept only
both reversed bare pairs recovering >=0.10 s/step, with clean numerical and
observer gates. `results.json` preserves an active run/error if something fails.

## Collect once the job stops; transfer from bbserv to murier

```bash
cd "$HOME"
source "$HOME/stage203-state.sh"
python3 "$kit/collect_cluster_logs_203.py" --stage 203 --job "$job" --results "$out" --tree "$tree"
scp "$HOME/stage203-review-${job}.tar.gz" kevlahan@murier.math.mcmaster.ca:~/Desktop/
```

Return that one archive. It includes results, raw phase/PMU evidence, work ledgers
in run logs, build identities/logs and Slurm output; large executable/checkpoint
payloads are excluded. Keep the results directory until review is complete.
