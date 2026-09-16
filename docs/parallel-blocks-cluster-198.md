# Stage198: diagnosis-only cluster measurements

This kit preserves validated Stage197 numerical code and the pinned legacy.
Stage196 is an additional untouched control. No new optimization is included.
The campaign measures the full step, including adaptation/reconstruction, and
outputs a signed operation ranking plus calibration and unresolved limitations.

Seven fresh J5/RK4 builds on bb13 share physics objects and FP_CONTRACT=off.
The seventh binary has optional storage-access probes, absent from the ordinary
observer build. Default campaign: 27 short executions (two checked, three
warmups, six bare timing, four minimal-observer, four operation-timing, four
work-count, two storage-probe controls and two storage-sampled runs). The two
observation rounds reverse implementation order. This replaces the old 58-run
campaign; these additional passes check observation effects, not precise CPU
statistics for a small optimization.

Transfer from the local machine:

```sh
scp /private/tmp/wavetrisk-stage198-cluster-kit.tar.gz kevlahan@bbserv:~/
```

The observer's small monotonic-clock C helper is compiled with `mpicc`, ensuring
that bbserv uses the configured MVAPICH/GCC header search paths.

On bbserv, create a fresh worktree and build:

```sh
cd "$HOME"
tar -xzf wavetrisk-stage198-cluster-kit.tar.gz
kit="$HOME/wavetrisk-stage198-cluster-kit"
tree="$HOME/wt-stage198"
fixture="$HOME/stage187-fixture-20260911-145307"
state="$HOME/stage198-state.sh"
printf 'kit=%q\ntree=%q\nfixture=%q\n' "$kit" "$tree" "$fixture" > "$state"
git -C "$HOME/wav" bundle verify "$kit/stage196-baseline.bundle"
git -C "$HOME/wav" fetch "$kit/stage196-baseline.bundle" HEAD
git -C "$HOME/wav" worktree add --detach "$tree" 9bc3a19556dab79b3cba3944cf5df53517389da3
build_job=$(sbatch --parsable "$kit/build_198.sbatch" "$kit" "$tree")
printf 'build_job=%q\n' "$build_job" >> "$state"
tail -F "$HOME/wavetrisk-198-build-${build_job}.log"
```

Wait for the seven-build PASS. Pending jobs may not have a log yet; inspect
`squeue -j "$build_job"` before resubmitting. Stop tail with Ctrl-C, then:

```sh
cd "$HOME"
state="$HOME/stage198-state.sh"
source "$state"
out="$HOME/stage198-nk-results-$(date +%Y%m%d-%H%M%S)"
job=$(sbatch --parsable "$kit/cluster_198.sbatch" "$kit" "$tree" "$fixture" "$out")
printf 'out=%q\njob=%q\n' "$out" "$job" >> "$state"
tail -F "$HOME/wavetrisk-198-${job}.log"
```

After completion or failure, collect and transfer one archive:

```sh
cd "$HOME"
state="$HOME/stage198-state.sh"
source "$state"
python3 "$kit/collect_cluster_logs_198.py" --stage 198 --job "$job" --results "$out" --tree "$tree"
scp "$HOME/stage198-review-${job}.tar.gz" kevlahan@murier.math.mcmaster.ca:~/Desktop/
```

The archive includes ranking.md, ranking.csv, results.json, provenance, run logs
and raw common/operation/dependency/storage records. Checkpoints and executables
are excluded. A completed simulation campaign does not automatically establish
causality or optimization readiness. Review calibration flags, fused-scope
limitations and unclassified time before using the ranked locations as causes.
