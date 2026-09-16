# Stage 201 cluster instructions

Use the new kit and a fresh source tree. Keep the Stage 200 builds/results for
comparison. This diagnostic uses exclusive bb13/nk, 83 tasks across two sockets,
unchanged J5 numerics, and 20 solver runs. Do not edit Makefiles or profiler files
after preparing the kit/builds. No sudo or kernel configuration change is needed
or performed; inaccessible counters cause an early stop with collectable logs.

## Transfer from murier

Run on the local machine:

```bash
scp /private/tmp/wavetrisk-stage201-cluster-kit.tar.gz kevlahan@bbserv:~/
```

## Prepare and build on bbserv

```bash
cd "$HOME"
tar -xzf wavetrisk-stage201-cluster-kit.tar.gz
kit="$HOME/wavetrisk-stage201-cluster-kit"
tree="$HOME/wt-stage201"
fixture="$HOME/stage187-fixture-20260911-145307"
state="$HOME/stage201-state.sh"
printf 'kit=%q\ntree=%q\nfixture=%q\n' "$kit" "$tree" "$fixture" > "$state"

git -C "$HOME/wav" bundle verify "$kit/stage196-baseline.bundle"
git -C "$HOME/wav" fetch "$kit/stage196-baseline.bundle" HEAD
git -C "$HOME/wav" worktree add --detach "$tree" 9bc3a19556dab79b3cba3944cf5df53517389da3

build_job=$(sbatch --parsable "$kit/build_201.sbatch" "$kit" "$tree")
printf 'build_job=%q\n' "$build_job" >> "$state"
echo "Build job: $build_job"
tail -F "$HOME/wavetrisk-201-build-${build_job}.log"
```

A missing log immediately after submission is normal while the job is pending.
Use `squeue -j "$build_job"` if necessary. Wait for the five-build PASS message;
Ctrl-C stops tail, not the build. Do not submit a second build into the same tree.
The detached baseline is only the reconstruction base; the kit applies the
unchanged numerical candidate and then the diagnostic observer patches.

## Submit tests

```bash
cd "$HOME"
source "$HOME/stage201-state.sh"
out="$HOME/stage201-nk-results-$(date +%Y%m%d-%H%M%S)"
job=$(sbatch --parsable "$kit/cluster_201.sbatch" "$kit" "$tree" "$fixture" "$out")
printf 'out=%q\njob=%q\n' "$out" "$job" >> "$state"
echo "Test job: $job"
echo "Results: $out"
tail -F "$HOME/wavetrisk-201-${job}.log"
```

The first step checks PMU access and sampling capability, before running climate.
If it fails, return the review archive; it contains `preflight/capabilities.json`,
compiler/perf logs, event availability and scheduling evidence. Do not weaken the
counter checks or replace missing measurements with zeros. The necessary support
request will be specific to the observed failure (for example perf installation
or permission for user-space per-process hardware events).

After successful preflight, monitor the same `wavetrisk-201-JOB.log`. Individual
solver output is in `"$out"/n83-s2-*/run.log`. Successful completion reports
Stage 201 hardware campaign PASS and writes `results.json`,
`hardware-summary.json` and `review.md`. Numerical PASS does not certify that
all observer overhead or sparse-sample review flags are acceptable.

## Collect and return one archive

After the job stops (success or failure):

```bash
cd "$HOME"
source "$HOME/stage201-state.sh"
python3 "$kit/collect_cluster_logs_201.py" \
  --stage 201 --job "$job" --results "$out" --tree "$tree"
scp "$HOME/stage201-review-${job}.tar.gz" \
  kevlahan@murier.math.mcmaster.ca:~/Desktop/
```

Return that archive. It includes raw perf recordings, on-cluster symbolized
stacks, per-rank phase counters/timelines, calibration/quality flags, placement,
provenance, build logs and job output. It excludes solver binaries and checkpoints.
Keep the source/build tree and binaries until review finishes: re-symbolization
of raw perf data can require the exact binary/shared-library versions. The
exported text supports the normal review without installing cluster libraries
on the local machine.

A repeat must use a new `out` directory. Reuse successful builds only while the
kit and source hashes remain unchanged. No second rank-count sweep is required.
