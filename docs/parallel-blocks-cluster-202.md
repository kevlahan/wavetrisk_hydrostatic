# Stage202 cluster prototype test

Use the supplied `wavetrisk-stage202-cluster-kit.tar.gz`. Build and run on exclusive
`nk`/bb13, using the existing J5 fixture. Keep `-march=native`: compilation and
execution occur on the same verified CPU. Do not edit either Makefile or input
between building and testing.

## Transfer from murier

```bash
scp /private/tmp/wavetrisk-stage202-cluster-kit.tar.gz kevlahan@bbserv:~/
```

## Prepare and build on bbserv

```bash
cd "$HOME"
tar -xzf wavetrisk-stage202-cluster-kit.tar.gz
kit="$HOME/wavetrisk-stage202-cluster-kit"
tree="$HOME/wt-stage202"
fixture="$HOME/stage187-fixture-20260911-145307"
state="$HOME/stage202-state.sh"
printf 'kit=%q\ntree=%q\nfixture=%q\n' "$kit" "$tree" "$fixture" > "$state"
git -C "$HOME/wav" bundle verify "$kit/stage196-baseline.bundle"
git -C "$HOME/wav" fetch "$kit/stage196-baseline.bundle" HEAD
git -C "$HOME/wav" worktree add --detach "$tree" 9bc3a19556dab79b3cba3944cf5df53517389da3
build_job=$(sbatch --parsable "$kit/build_202.sbatch" "$kit" "$tree")
printf 'build_job=%q\n' "$build_job" >> "$state"
tail -F "$HOME/wavetrisk-202-build-${build_job}.log"
```

The builder applies the frozen baseline, prototype and matched observer patches.
Wait for `PASS: six fresh J5 builds...`; Ctrl-C only exits `tail`. A missing log
before Slurm starts the job is normal. Inspect `squeue -j "$build_job"` if waiting.
If the script reports an error, retain its build log and use a fresh tree after
resolving it; do not bypass identity checks.

## Submit tests

```bash
cd "$HOME"
source "$HOME/stage202-state.sh"
out="$HOME/stage202-nk-results-$(date +%Y%m%d-%H%M%S)"
job=$(sbatch --parsable "$kit/cluster_202.sbatch" "$kit" "$tree" "$fixture" "$out")
printf 'out=%q\njob=%q\n' "$out" "$job" >> "$HOME/stage202-state.sh"
tail -F "$HOME/wavetrisk-202-${job}.log"
```

Watch the main Slurm log above. It prints each run and PASS. On completion,
`results.json`, `prototype-summary.json` and `REVIEW.md` distinguish correctness
from the substantial-gain decision. Numerical PASS alone is not approval to merge.

## Collect and return one archive after completion (or failure)

```bash
cd "$HOME"
source "$HOME/stage202-state.sh"
python3 "$kit/collect_cluster_logs_202.py" --stage 202 --job "$job" --results "$out" --tree "$tree"
scp "$HOME/stage202-review-${job}.tar.gz" kevlahan@murier.math.mcmaster.ca:~/Desktop/
```

Return that archive. It contains build identities, Slurm/run logs, numerical gates,
raw rank counters, placement and memory evidence; it excludes executables and
checkpoint payloads. The collector refuses to overwrite an existing archive.
