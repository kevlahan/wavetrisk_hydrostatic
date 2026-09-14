# Stage 190: cluster comparison with legacy Domain

This kit compares Stage 189, Stage 190 and legacy Domain numerics using the
same J5/J7 checkpoint-3 fixture, 83 ranks and one allocation. Stage 190 must
match Stage 189 and checked/oracle checkpoints exactly. Legacy uses the existing
fixed numerical screen (velocity/mass 1e-10, temperature 1e-7, soil exact), with
exact topology and thresholds. The screen is not a proof of legacy equivalence.
A failed numerical or printed-workload gate stops the campaign.

The legacy executable is explicitly **legacy Domain plus the common checkpoint
I/O fix**. Its numerical sources come from authoritative revision
`b82647b36eb825aed0c0aae4e92adcf4ab46d608`. Only `src/checkpoint.f90` changes in a
separate experimental archive: Stage 187's validated ROMIO hints, truncation
barrier, write synchronization, transfer-count and persisted-directory checks.
The patch changes no field packing, payload offsets, arithmetic or checkpoint
format. The original archive and its complete tracked-content manifest remain
in the kit. No existing legacy checkout is modified. This prevents the known
bbserv checkpoint corruption from confounding the comparison. Full elapsed
ratios apply to this I/O-normalized legacy build; ordinary-step ratios are the
primary measure of the solver gap.

## Host and placement

The user reports `nk`/`bb13` has two AMD EPYC 9754 sockets, each with 128 physical
cores and one thread per core. Each `bb` node has two EPYC 7543 sockets with 32
physical cores each and one thread per core. Thus 83 ranks fit on one `nk` node,
whereas the same rank count spans two `bb` nodes. Shared memory still has NUMA
placement effects. The protocol records actual rank affinity, allowed memory
nodes, CPU model, CPU/core/socket/NUMA topology, per-rank peak RSS and page faults,
and per-node VM counters before/after each run. Allowed memory nodes are not a
measurement of actual page placement. Affinity sets must not overlap.

Both build and test batch scripts default to `nk`, one exclusive node. The build
script runs compilation there so `-march=native` targets the compute CPU. All
optimized builds use `-O2`, `-march=native` and `-ffp-contract=off`; the checked
block build uses `-O0` with runtime checks and contraction off. The legacy
Makefile receives the contraction flag through a separate include. All 23
external physics object hashes must match across the four builds.

## Prepare on bbserv

Transfer and extract the entire `wavetrisk-stage190-cluster-kit.tar.gz` under
`$HOME`. Load the same GCC/MVAPICH/netCDF/LAPACK/Python environment used for the
passing Stage 189 campaign. Python 3.11.9+, NumPy and zstd are required. The corrected runner uses
Python `os.wait4` for child resource accounting; GNU time is not required. The kit contains its Python helpers and the original legacy
source archive, so no additional legacy checkout is needed.

```sh
kit="$HOME/wavetrisk-stage190-cluster-kit"
tree="$HOME/wt-stage190"
fixture="$HOME/stage187-fixture-20260911-145307"

# The bundle supplies the committed Stage 189 baseline, even before it is pushed.
# Its prerequisite is the already merged Stage 188 commit 61509429.
git -C "$HOME/wav" bundle verify "$kit/stage189-baseline.bundle"
git -C "$HOME/wav" fetch "$kit/stage189-baseline.bundle" HEAD
git -C "$HOME/wav" worktree add --detach "$tree" \
  a0c75adb4dc1741f56ecb572ab1809e8a70a32e4

sbatch "$kit/build_190.sbatch" "$kit" "$tree"
```

Use a new worktree name if this one already exists; do not reset or clean an
existing tree. The builder also creates the sibling `$tree-legacy190` and
rejects an existing directory with that name. It builds Stage 189 first,
applies the Stage 190 patch in the new worktree, builds optimized/checked
candidates, then prepares and builds the separate legacy archive. Compilation
is serial and all build directories are fresh. It rejects environment flag
overrides and changed source/Makefile identities.

Monitor `wavetrisk-190-build-JOBID.log`. Successful completion prints
`PASS: four fresh J5 builds ...` and writes `$tree/restriction-190-builds.json`
with status `passed`. Detailed compiler commands are in the three
`$tree/build-190-*.log` files and `$tree-legacy190/build-190-legacy.log`.

## Run after the build passes

The immutable fixture must retain `resume=3`, `max_level=7`, `time_end=0.315`,
`dt_write=0.01`, `CP_EVERY=1`, 30 atmospheric and 10 soil layers, and a real
`grids/HR95JT_WT_004`. The campaign copies inputs and grids into its own snapshot.
It requires checkpoint 4, its reload, and a subsequent vertical remap.

```sh
out="$HOME/stage190-nk-results-$(date +%Y%m%d-%H%M%S)"
sbatch "$kit/cluster_190.sbatch" "$kit" "$tree" "$fixture" "$out"
```

Monitor `wavetrisk-190-JOBID.log`; individual solver logs are
`$out/RUN_LABEL/run.log`. The 16 runs are checked-off/on, three warmups, three
rotating timing triples, and two separate block profiles. The triples are
baseline/candidate/legacy, candidate/legacy/baseline, legacy/baseline/candidate.
No checked, warmup or profile time enters the medians. A Python parent waits
for each solver process and records its child-specific RSS/fault counts through
`os.wait4`, preserving inherited MPI descriptors. The waiting parents add some
memory outside the measured solver-child RSS. Wrapper/startup costs remain
in elapsed time but not the solver's printed ordinary-step CPU timing. All
three timed variants use the same wrapper; node VM records include its overhead.

`results.json` is updated after each run and gate. Accept completion only when
its final status is `passed`, and the job exits successfully. It reports the
Stage 190/189 improvement and both block/legacy elapsed and ordinary-step ratios.
A ratio above 1 means block is slower. Inspect VM/swap evidence and rank placement
before accepting a performance claim; the script does not infer a universal
speedup or a two-node result from one campaign.

Return:

- `$out/results.json` and `wavetrisk-190-JOBID.log`.
- `$out/baseline-profile/run.log` and `$out/candidate-profile/run.log`.
- `$tree/restriction-190-builds.json` and all four build logs.

The result JSON embeds rank resource/placement and node VM records. Keep raw
checkpoints until the numerical gates have been reviewed. For another campaign,
reuse the binaries and fixture but choose a new `out` directory. Do not compare
absolute `nk` times directly with older two-node `bb` timings.

## Runner correction after job 54552

The initial kit required `/usr/bin/time`, which is absent on `bb13`. Job 54552
failed in `rank_launch.py` before starting the solver. Its original error was
then masked by a missing-resource-file read. The runner-only update replaces
this dependency with `rank_launch_190.py` and checks the launch exit status
before reading resource reports. No Fortran, executable, fixture or build
identity changes are needed. Apply the small `wavetrisk-stage190-runner-fix`
update to an already extracted original kit, then submit with a new output
directory. The separate legacy I/O fix and every numerical gate are unchanged.
