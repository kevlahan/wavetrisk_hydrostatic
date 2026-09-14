# Stage 194 cluster protocol

Compare the validated Stage 193 scalar-column kernel, Stage 194 with reused
inner-vector geometry, and pinned legacy Domain code with the common checkpoint
I/O fix. All block cases keep Stage 193 scalar columns on and Stage 192 sparse
scatter off. The Stage 194 kernel is explicit and reversible.

The fresh source tree starts at the committed Stage 191 bootstrap. The build
script first applies the exact frozen Stage 193 patch, builds the previous-stage
baseline, then applies the Stage 194 delta and builds the candidate and checked
candidate. No prior uncommitted state or merge is required on bbserv.

## Transfer and build

From the local machine, transfer the complete archive:

```sh
scp /private/tmp/wavetrisk-stage194-cluster-kit.tar.gz kevlahan@bbserv:~/
```

On bbserv, load the same GCC/MVAPICH/netCDF/LAPACK/Python environment as Stage 191.
Python 3.11.9+, NumPy and zstd are required. Use new kit/tree names; preserve
previous results and the existing fixture. Run from HOME so Slurm logs land there.

```sh
cd "$HOME"
tar -xzf wavetrisk-stage194-cluster-kit.tar.gz
kit="$HOME/wavetrisk-stage194-cluster-kit"
tree="$HOME/wt-stage194"
fixture="$HOME/stage187-fixture-20260911-145307"
state="$HOME/stage194-state.sh"
printf 'kit=%q\ntree=%q\nfixture=%q\n' "$kit" "$tree" "$fixture" > "$state"
git -C "$HOME/wav" bundle verify "$kit/stage191-baseline.bundle"
git -C "$HOME/wav" fetch "$kit/stage191-baseline.bundle" HEAD
git -C "$HOME/wav" worktree add --detach "$tree" 30c831149653205210531427f81b184663a3acbf
build_job=$(sbatch --parsable "$kit/build_194.sbatch" "$kit" "$tree")
printf 'build_job=%q\n' "$build_job" >> "$state"
echo "Build job: $build_job"
tail -F "$HOME/wavetrisk-194-build-${build_job}.log"
```

The exclusive nk job builds on bb13. Keep `-march=native`: both compilation and
execution use this CPU. All optimized executables use O2 and FP_CONTRACT=off;
the checked executable uses O0 and runtime checks. Five fresh J5 builds include
baseline, candidate, checked candidate, legacy and separately instrumented legacy.
All 23 external physics objects must match. Legacy directories are
`$tree-legacy194` and `$tree-legacy194-profile`.

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
source "$HOME/stage194-state.sh"
out="$HOME/stage194-nk-results-$(date +%Y%m%d-%H%M%S)"
job=$(sbatch --parsable "$kit/cluster_194.sbatch" "$kit" "$tree" "$fixture" "$out")
printf 'out=%q\njob=%q\n' "$out" "$job" >> "$HOME/stage194-state.sh"
echo "Test job: $job"
echo "Results: $out"
tail -F "$HOME/wavetrisk-194-${job}.log"
```

There are 28 runs: two checked modes (oracles off/on), three warmups, six rotating
Stage 193/194/legacy timing triples, a candidate with the new kernel disabled, and
four separate profiles: the original Stage 193 baseline, Stage 194 with the
geometry cache disabled, Stage 194 with it enabled, and instrumented legacy. Only six timing runs per implementation enter the main
medians. All block runs explicitly use installation=1 and mass overlap=0.
`WAVETRISK_STAGE193_SCALAR_COLUMNS=1` remains enabled for every block run.
`WAVETRISK_STAGE194_INNER_GEOMETRY=1` enables the new cache, including during
full oracle validation. Both candidate-control and candidate-control-profile
select 0. The new cache defaults to off pending validation. No manual environment
exports are needed; the campaign sets and checks these options.

The campaign stops on failed exact block/checkpoint gates, unequal workload,
missing restart/remap, or failed fixed legacy screen. Legacy velocity/mass limits
remain 1e-10, temperature 1e-7, with exact soil/topology/thresholds. Instrumented
legacy must match uninstrumented legacy exactly. Build identities, core affinity,
rank resources and node memory/swap evidence are retained.

`inverse-comparison.json` and `results.json` report complete WT_after_step,
inverse transfer counts, scalar work and the new split. `scalar`, `inner_vector`
and the existing `outer_vector` child scopes are separate; the combined parent
is not added to its children. An explicit conservation check verifies those
children and the parent overhead account for the same combined scope.

The old combined timer also includes the outer-vector kernel, despite its old
label. Stage 193 retains its original 65-region profile. The Stage 194
candidate-control-profile supplies the 67-region split with caching disabled;
candidate-profile supplies the identical split with caching enabled. Thus the
optimization and added timer overhead can be examined independently, while the
uninstrumented timing comparison remains a direct Stage 193 versus Stage 194
comparison. Four profile runs are excluded from headline medians.

The report checks identical inverse bytes/messages and physical cell work, and
measures geometry reads against the actual cache-disabled control. Profiles
include the checkpoint window; only separately reported ordinary-step diagnostics
exclude it. Repeated uninstrumented ordinary-step gain determines acceptance.

## Collect one archive, then transfer it

After the job stops (including on failure):

```sh
cd "$HOME"
source "$HOME/stage194-state.sh"
python3 "$kit/collect_cluster_logs.py" --stage 194 --job "$job" --results "$out" --tree "$tree"
scp "$HOME/stage194-review-${job}.tar.gz" kevlahan@murier.math.mcmaster.ca:~/Desktop/
```

Return that single archive. It includes results, all run/profile logs, raw legacy
rank profiles, rank resource/affinity records, five compiler logs, build manifest
and main Slurm log, with a hash inventory. Checkpoints/binaries stay on the
cluster until review. To repeat, reuse builds/fixture and choose a new `out`.
