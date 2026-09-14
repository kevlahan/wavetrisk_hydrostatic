# Stage 191: installation ownership and mass exchange overlap

The candidate compares directly with committed Stage 190 (`cdd2eb64`) and pinned
legacy Domain (`b82647b36eb825aed0c0aae4e92adcf4ab46d608`) with the common checkpoint
I/O fix. It transfers temporary block allocations at installation, omits redundant
serialization proofs in production, and starts mass exchange before independent
velocity work. Full oracle mode retains copying/proofs and synchronous exchange.
It does not yet reuse unchanged subtrees or remove the Domain/block round trip.

## Prepare and build on bbserv

Transfer and extract the complete `wavetrisk-stage191-cluster-kit.tar.gz` under
`$HOME`. Load the same GCC/MVAPICH/netCDF/LAPACK/Python environment as Stage 190.
Python 3.11.9+, NumPy and zstd are required. Use fresh tree names; no old tree is
reset or cleaned. The bundle requires the already merged Stage 189 history.

```sh
cd "$HOME"
kit="$HOME/wavetrisk-stage191-cluster-kit"
tree="$HOME/wt-stage191"
fixture="$HOME/stage187-fixture-20260911-145307"
git -C "$HOME/wav" bundle verify "$kit/stage190-baseline.bundle"
git -C "$HOME/wav" fetch "$kit/stage190-baseline.bundle" HEAD
git -C "$HOME/wav" worktree add --detach "$tree" cdd2eb64
build_job=$(sbatch --parsable "$kit/build_191.sbatch" "$kit" "$tree")
echo "Build job: $build_job"
tail -F "$HOME/wavetrisk-191-build-${build_job}.log"
```

The exclusive `nk` allocation builds on bb13 so `-march=native` targets the test
CPU. All optimized builds use O2 and contraction off; the checked build uses O0
and runtime checks. Compilation is serial. The builder produces five fresh J5
executables: Stage 190, Stage 191 optimized/checked, authoritative legacy, and a
separate legacy copy with phase timers. The latter two live in sibling directories
`$tree-legacy191` and `$tree-legacy191-profile`; neither may already exist.
All 23 external physics object hashes must match. Source/patch/binary manifests
are verified. Wait for `PASS: five fresh J5 builds ...` and status `passed` in
`$tree/restriction-191-builds.json` before submitting the tests.

## Test after a successful build

Use the unchanged J5/J7 fixture: resume=3, time_end=0.315, dt_write=0.01,
CP_EVERY=1, 30 atmosphere/10 soil layers and a real `grids/HR95JT_WT_004`.
The runner snapshots the fixture and requires newly written/reloaded checkpoint
4 and subsequent remapping. It uses 83 ranks, one exclusive nk node, fixed
nonoverlapping core binding and one thread per rank.

```sh
cd "$HOME"
out="$HOME/stage191-nk-results-$(date +%Y%m%d-%H%M%S)"
job=$(sbatch --parsable "$kit/cluster_191.sbatch" "$kit" "$tree" "$fixture" "$out")
echo "Test job: $job"
echo "Results: $out"
tail -F "$HOME/wavetrisk-191-${job}.log"
```

There are 33 runs: two checked modes; three warmups; six rotating timing triples
(Stage 190/191/legacy); one candidate with both changes disabled; three alternating
pairs with installation-only/overlap-only; and three separate profiles. Only the
six headline runs per implementation enter the main medians. Ablation medians
are reported separately. Profiles, control, warmups and checked runs are excluded.
The runtime controls `WAVETRISK_STAGE191_INSTALL` and
`WAVETRISK_STAGE191_MASS_OVERLAP` are explicitly set per run and rank-consistent.

Stage 191, its ablations and checked/oracle results must match Stage 190 exactly
in the semantic checkpoint comparison. The legacy numerical screen stays fixed:
velocity/mass 1e-10, temperature 1e-7, soil/topology/thresholds exact. Instrumented
legacy must match the authoritative legacy checkpoint exactly. Every run must
match the printed time/dt/Jmax/dof trajectory. These gates stop the campaign on
failure; legacy tolerance is a screen, not a proof of mathematical equivalence.

Accept only final status `passed` in `results.json` and successful job exit.
Inspect rank affinity, per-rank RSS/page faults and node VM/swap records before
interpreting timing. Compare within this allocation, not against historical times.
The legacy profile embeds every rank's aggregate phase comparison in results;
region IDs are shared but some implementation boundaries differ. Raw rank files
remain useful for investigation. No profile time enters headline ratios.

Return:

- `$out/results.json` and `wavetrisk-191-JOBID.log`.
- `run.log` from baseline-profile, candidate-profile and legacy-profile, renamed
  to identify each run.
- `$tree/restriction-191-builds.json`, three `$tree/build-191-*.log` files, plus
  `$tree-legacy191/build-191-legacy.log` and
  `$tree-legacy191-profile/build-191-legacy-profile.log`.
- An archive of `$out/legacy-profile/legacy-detail-rank-*.txt` for detailed analysis.

Retain raw checkpoints until review. For another campaign reuse the binaries and
fixture with a new output directory. A longer adaptive-window campaign remains a
follow-up after the short case passes; its save interval and expected checkpoint
set must be chosen together from a pilot, not changed independently.

## Runner correction after job 54560

Job 54560 failed before the campaign started with `OSError: [Errno 7]
Argument list too long`. The submit helper serialized the entire five-build
provenance record into one `--build-notes` argument. The runner update uses
`--build-notes-file` to pass the path of `restriction-191-builds.json` instead.
The campaign reads and embeds the complete record in `results.json`; source and
binary identity checks, numerical gates and timing settings are unchanged.

Install both updated Python files from `wavetrisk-stage191-runner-fix.tar.gz`
into the existing kit, then resubmit with a new output directory. No executable
rebuild is needed. The complete kit archive also incorporates this correction.
A regression test spawns a real Python child with a build record larger than
512 KiB and verifies that it arrives intact through the file argument.

## Default evidence collection

After job 54561, use the single-archive collection process documented in
`parallel-blocks-cluster-log-collection.md` for future campaigns. Future kits
include `collect_cluster_logs.py`; transfer the resulting archive to
`kevlahan@murier.math.mcmaster.ca:~/Desktop/`.
