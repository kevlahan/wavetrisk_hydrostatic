# Stage 193 cluster protocol

Compare accepted Stage 191 production (`30c83114`: installation on, mass
overlap off), Stage 193 scalar inverse columns, and pinned legacy Domain
`b82647b36eb825aed0c0aae4e92adcf4ab46d608` with the common checkpoint I/O fix.
The candidate hoists stencil address and geometry preparation across physical
layers and scalar fields. Stage 192 sparse scatter remains OFF throughout this
campaign: the primary comparison isolates computational work, using the accepted
production baseline rather than the inconclusive Stage 192 experiment.

## Transfer and build

From the local machine, transfer the complete archive:

```sh
scp /private/tmp/wavetrisk-stage193-cluster-kit.tar.gz kevlahan@bbserv:~/
```

On bbserv, load the same GCC/MVAPICH/netCDF/LAPACK/Python environment as Stage 191.
Python 3.11.9+, NumPy and zstd are required. Use new kit/tree names; preserve
previous results and the existing fixture. Run from HOME so Slurm logs land there.

```sh
cd "$HOME"
tar -xzf wavetrisk-stage193-cluster-kit.tar.gz
kit="$HOME/wavetrisk-stage193-cluster-kit"
tree="$HOME/wt-stage193"
fixture="$HOME/stage187-fixture-20260911-145307"
state="$HOME/stage193-state.sh"
printf 'kit=%q\ntree=%q\nfixture=%q\n' "$kit" "$tree" "$fixture" > "$state"
git -C "$HOME/wav" bundle verify "$kit/stage191-baseline.bundle"
git -C "$HOME/wav" fetch "$kit/stage191-baseline.bundle" HEAD
git -C "$HOME/wav" worktree add --detach "$tree" 30c831149653205210531427f81b184663a3acbf
build_job=$(sbatch --parsable "$kit/build_193.sbatch" "$kit" "$tree")
printf 'build_job=%q\n' "$build_job" >> "$state"
echo "Build job: $build_job"
tail -F "$HOME/wavetrisk-193-build-${build_job}.log"
```

The exclusive nk job builds on bb13. Keep `-march=native`: both compilation and
execution use this CPU. All optimized executables use O2 and FP_CONTRACT=off;
the checked executable uses O0 and runtime checks. Five fresh J5 builds include
baseline, candidate, checked candidate, legacy and separately instrumented legacy.
All 23 external physics objects must match. Legacy directories are
`$tree-legacy193` and `$tree-legacy193-profile`.

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
source "$HOME/stage193-state.sh"
out="$HOME/stage193-nk-results-$(date +%Y%m%d-%H%M%S)"
job=$(sbatch --parsable "$kit/cluster_193.sbatch" "$kit" "$tree" "$fixture" "$out")
printf 'out=%q\njob=%q\n' "$out" "$job" >> "$HOME/stage193-state.sh"
echo "Test job: $job"
echo "Results: $out"
tail -F "$HOME/wavetrisk-193-${job}.log"
```

There are 27 runs: two checked modes (oracles off/on), three warmups, six rotating
Stage 191/193/legacy timing triples, a candidate with the new kernel disabled, and
three separate profiles. Only six timing runs per implementation enter the main
medians. All block runs explicitly use installation=1 and mass overlap=0.
`WAVETRISK_STAGE193_SCALAR_COLUMNS=1` enables the candidate, including in full
oracle validation. The candidate-control run selects 0. Both production defaults
remain conservative: the new kernel and Stage 192 sparse scatter are off unless
explicitly enabled. All candidate/checked logs must identify the requested mode.

The campaign stops on failed exact block/checkpoint gates, unequal workload,
missing restart/remap, or failed fixed legacy screen. Legacy velocity/mass limits
remain 1e-10, temperature 1e-7, with exact soil/topology/thresholds. Instrumented
legacy must match uninstrumented legacy exactly. Build identities, core affinity,
rank resources and node memory/swap evidence are retained.

`inverse-comparison.json` and `results.json` report full inclusive WT_after_step
and gather/alias/scatter bytes/messages. The campaign requires those transfer
counts to remain identical. It also records combined scalar/inner-vector kernel
time and scalar-inverse work counts (address resolution and geometry reads),
including reference-equivalent counts for the unchanged output workload. Whole-window profile times include the
checkpoint step and do not enter headline medians. Transfer helpers are contained
in inverse time, not added to it. Compare the current baseline and candidate,
not the elapsed seconds from older jobs. Fewer geometry reads are not automatically
a whole-step performance gain; acceptance requires repeatable ordinary-step improvement.

## Collect one archive, then transfer it

After the job stops (including on failure):

```sh
cd "$HOME"
source "$HOME/stage193-state.sh"
python3 "$kit/collect_cluster_logs.py" --stage 193 --job "$job" --results "$out" --tree "$tree"
scp "$HOME/stage193-review-${job}.tar.gz" kevlahan@murier.math.mcmaster.ca:~/Desktop/
```

Return that single archive. It includes results, all run/profile logs, raw legacy
rank profiles, rank resource/affinity records, five compiler logs, build manifest
and main Slurm log, with a hash inventory. Checkpoints/binaries stay on the
cluster until review. To repeat, reuse builds/fixture and choose a new `out`.
