# Stage 192 cluster protocol

Compare the accepted Stage 191 configuration (`30c83114`: installation on, mass
overlap off), Stage 192 sparse outer-vector interior installation, and pinned
legacy Domain `b82647b36eb825aed0c0aae4e92adcf4ab46d608` with the common checkpoint
I/O fix. Stage 192 sends only written interior and boundary vector edges in outer scatter;
it preserves scaffold resets, alias writes and arithmetic order.

## Transfer and build

From the local machine, transfer the complete archive:

```sh
scp /private/tmp/wavetrisk-stage192-cluster-kit.tar.gz kevlahan@bbserv:~/
```

On bbserv, load the same GCC/MVAPICH/netCDF/LAPACK/Python environment as Stage 191.
Python 3.11.9+, NumPy and zstd are required. Use new kit/tree names; preserve
previous results and the existing fixture. Run from HOME so Slurm logs land there.

```sh
cd "$HOME"
tar -xzf wavetrisk-stage192-cluster-kit.tar.gz
kit="$HOME/wavetrisk-stage192-cluster-kit"
tree="$HOME/wt-stage192"
fixture="$HOME/stage187-fixture-20260911-145307"
state="$HOME/stage192-state.sh"
printf 'kit=%q\ntree=%q\nfixture=%q\n' "$kit" "$tree" "$fixture" > "$state"
git -C "$HOME/wav" bundle verify "$kit/stage191-baseline.bundle"
git -C "$HOME/wav" fetch "$kit/stage191-baseline.bundle" HEAD
git -C "$HOME/wav" worktree add --detach "$tree" 30c831149653205210531427f81b184663a3acbf
build_job=$(sbatch --parsable "$kit/build_192.sbatch" "$kit" "$tree")
printf 'build_job=%q\n' "$build_job" >> "$state"
echo "Build job: $build_job"
tail -F "$HOME/wavetrisk-192-build-${build_job}.log"
```

The exclusive nk job builds on bb13. Keep `-march=native`: both compilation and
execution use this CPU. All optimized executables use O2 and FP_CONTRACT=off;
the checked executable uses O0 and runtime checks. Five fresh J5 builds include
baseline, candidate, checked candidate, legacy and separately instrumented legacy.
All 23 external physics objects must match. Legacy directories are
`$tree-legacy192` and `$tree-legacy192-profile`.

Wait for `PASS: five fresh J5 builds ...`. Stop tail with Ctrl-C. A missing log
while the job is pending is normal; inspect `squeue -j "$build_job"` before
submitting another build. Build failure details remain in the Slurm and compiler logs.

## Run

Use the frozen J5/J7 fixture with resume=3, time_end=0.315, dt_write=0.01,
CP_EVERY=1, 30 atmospheric and 10 soil levels, and materializable optimized grids.
The kit requires exactly checkpoint 4 to be written/reloaded, with subsequent
remapping. Do not modify an input inside an existing results directory.

```sh
cd "$HOME"
source "$HOME/stage192-state.sh"
out="$HOME/stage192-nk-results-$(date +%Y%m%d-%H%M%S)"
job=$(sbatch --parsable "$kit/cluster_192.sbatch" "$kit" "$tree" "$fixture" "$out")
printf 'out=%q\njob=%q\n' "$out" "$job" >> "$HOME/stage192-state.sh"
echo "Test job: $job"
echo "Results: $out"
tail -F "$HOME/wavetrisk-192-${job}.log"
```

There are 27 runs: two checked modes (oracles off/on), three warmups, six rotating
Stage 191/192/legacy timing triples, a candidate with Stage 192 disabled, and
three separate profiles. Only six timing runs per implementation enter the main
medians. All block runs explicitly use installation=1 and mass overlap=0.
`WAVETRISK_STAGE192_SPARSE_OUTER=0` restores the Stage 191 scatter selection;
1 selects the new route, including in full oracle validation.

The campaign stops on failed exact block/checkpoint gates, unequal workload,
missing restart/remap, or failed fixed legacy screen. Legacy velocity/mass limits
remain 1e-10, temperature 1e-7, with exact soil/topology/thresholds. Instrumented
legacy must match uninstrumented legacy exactly. Build identities, core affinity,
rank resources and node memory/swap evidence are retained.

`inverse-comparison.json` and `results.json` report full inclusive WT_after_step
and gather/alias/scatter bytes/messages. Whole-window profile times include the
checkpoint step and do not enter headline medians. Transfer helpers are contained
in inverse time, not added to it. Compare the current baseline and candidate,
not the elapsed seconds from older jobs. A traffic reduction is not automatically
a performance gain; acceptance requires repeatable ordinary-step improvement.

## Collect one archive, then transfer it

After the job stops (including on failure):

```sh
cd "$HOME"
source "$HOME/stage192-state.sh"
python3 "$kit/collect_cluster_logs.py" --stage 192 --job "$job" --results "$out" --tree "$tree"
scp "$HOME/stage192-review-${job}.tar.gz" kevlahan@murier.math.mcmaster.ca:~/Desktop/
```

Return that single archive. It includes results, all run/profile logs, raw legacy
rank profiles, rank resource/affinity records, five compiler logs, build manifest
and main Slurm log, with a hash inventory. Checkpoints/binaries stay on the
cluster until review. To repeat, reuse builds/fixture and choose a new `out`.
