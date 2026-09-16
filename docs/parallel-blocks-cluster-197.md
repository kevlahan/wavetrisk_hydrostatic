# Stage 197 cluster campaign

Use the frozen kit with a new detached worktree at Stage196 commit 9bc3a195.
The six builds comprise untouched Stage196, optimized candidate, checked
candidate, candidate observer, untouched legacy plus common I/O fix, and
legacy observer. The same candidate observer with resident selector off
provides the Stage196 diagnostic control and must match its untouched binary.

Both Slurm scripts use nk, one exclusive node and 83 tasks. Build on bb13 to
keep march=native consistent. J5/RK4, contraction off, common physics objects.

From the local machine:

```sh
scp /private/tmp/wavetrisk-stage197-cluster-kit.tar.gz kevlahan@bbserv:~/
```

On bbserv:

```sh
cd "$HOME"
tar -xzf wavetrisk-stage197-cluster-kit.tar.gz
kit="$HOME/wavetrisk-stage197-cluster-kit"
tree="$HOME/wt-stage197"
fixture="$HOME/stage187-fixture-20260911-145307"
state="$HOME/stage197-state.sh"
printf 'kit=%q\ntree=%q\nfixture=%q\n' "$kit" "$tree" "$fixture" > "$state"
git -C "$HOME/wav" bundle verify "$kit/stage196-baseline.bundle"
git -C "$HOME/wav" fetch "$kit/stage196-baseline.bundle" HEAD
git -C "$HOME/wav" worktree add --detach "$tree" 9bc3a19556dab79b3cba3944cf5df53517389da3
build_job=$(sbatch --parsable "$kit/build_197.sbatch" "$kit" "$tree")
printf 'build_job=%q\n' "$build_job" >> "$state"
tail -F "$HOME/wavetrisk-197-build-${build_job}.log"
```

Wait for the six-build PASS. A pending job may not have a log yet; check
`squeue -j "$build_job"` before submitting again. Stop tail with Ctrl-C.
Then submit the 14-execution campaign:

```sh
cd "$HOME"
state="$HOME/stage197-state.sh"
source "$state"
out="$HOME/stage197-nk-results-$(date +%Y%m%d-%H%M%S)"
job=$(sbatch --parsable "$kit/cluster_197.sbatch" "$kit" "$tree" "$fixture" "$out")
printf 'out=%q\njob=%q\n' "$out" "$job" >> "$state"
tail -F "$HOME/wavetrisk-197-${job}.log"
```

The campaign has checked-off/on, three warmups, two alternating timing rounds
across Stage196/candidate/legacy, and one focused profile for each. Profiles
are excluded from timing medians. All exact/workload/legacy screens remain.
No CPU-time confidence interval is claimed from two rounds; add rounds only
if results are small or inconsistent. Compare profile overhead against bare
runs before interpreting the additional scalar counters as precise timing.

Collect after completion or failure, then transfer the single archive:

```sh
cd "$HOME"
state="$HOME/stage197-state.sh"
source "$state"
python3 "$kit/collect_cluster_logs.py" --stage 197 --job "$job" --results "$out" --tree "$tree"
scp "$HOME/stage197-review-${job}.tar.gz" kevlahan@murier.math.mcmaster.ca:~/Desktop/
```

Return that tarball, including results, build records, run logs, common
workflow records, operation records and dependency records. Checkpoint and
binary payloads remain excluded as in previous review kits.
