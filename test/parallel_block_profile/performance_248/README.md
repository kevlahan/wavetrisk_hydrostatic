# RK migration performance campaign on bbserv

This companion kit runs the qualified Stage 247 source in two matched builds. `feature-off` compares the complete opt-in RK migration against the same current source with the new RK, physics, soil, adaptation, lifecycle, coefficient and mirror-free features disabled. `mirrored` keeps the native migration enabled in both builds and isolates removal of the persistent Domain SOL/WAV arrays. The candidate always enables the qualified Stage 247 feature set. These comparisons do not reuse an old binary from another compiler or node.

The Slurm script requests one exclusive `nk` node, which is bb13 with 256 effective CPUs. It compiles both variants there with `make -j1`, pins one MPI rank to each physical core and checks the placement before solving. `nk` contains only bb13, so no node exclusion is needed. The local machine's four-core limit applies only to local preparation; no local timing is used as cluster performance evidence.

For J4, correctness runs use 16 and 41 ranks; the 41-rank case includes a compute rank without Domain geometry. Timing screens use 16 and 41 ranks. For J5, correctness uses 16 and 167 ranks, while timing screens use 16 and 83 ranks. J5 requires the previously established `stage187-fixture-20260911-145307` directory on bbserv. The script checks that fixture's input, seed checkpoint and grid hashes against the Stage 230 accepted campaign. J4 uses the source kit's sealed fixture.

Every full run must complete 12 steps, two remaps, the checkpoint reload, 62 exact VTK members and the accepted numerical checkpoint gate against a same-rank baseline. Step-state records must agree exactly. Candidate runs also prove legacy SOL/WAV storage is absent between bounded output operations. Scientific failure stops further runs. Timing uses the high-precision median of nine ordinary steps per **independent solver run**, excluding startup, checkpoint/output and first post-reload steps (indices 0, 7 and 8). A 16-rank pair and a larger-rank pair screen for a 3% regression; a clear regression stops the campaign. Otherwise an A-B-B-A confirmation at the larger rank count follows, with at most one repeat batch if the result straddles the guard. The two screening pairs are excluded from confirmation. No branch promotion occurs automatically.

## 1. Transfer the two archives from the local machine

```bash
scp /Users/kevlahan/.codex/worktrees/7173/wavetrisk_hydrostatic/.stage247-validation/bbserv-candidate2.tar.gz \
    /Users/kevlahan/.codex/worktrees/7173/wavetrisk_hydrostatic/.stage248-validation/rk248-performance-kit.tar.gz \
    kevlahan@bbserv.math.mcmaster.ca:~/
```

The source archive SHA-256 is `dc6daba22953f5ab54398f6085b9b39aa04a2161b233fade692fedaa9d5d7d63`. Find the companion archive SHA-256 in `docs/parallel-blocks-stage248-performance-kit.json`. On bbserv, load the same GNU Fortran, MPI, NetCDF, Python/numpy, `numactl`, `lscpu` and `zstd` environment used by recent campaigns. Check hashes and extract into fresh names:

```bash
cd "$HOME"
sha256sum bbserv-candidate2.tar.gz rk248-performance-kit.tar.gz
[ ! -e bbserv-candidate2 ] && [ ! -e rk248-performance-kit ]
tar -xzf bbserv-candidate2.tar.gz
tar -xzf rk248-performance-kit.tar.gz
```

## 2. Submit one comparison

Start with the J4 feature-off comparison. Use a new results directory for each campaign:

```bash
cd "$HOME"
kit="$HOME/bbserv-candidate2"
perf="$HOME/rk248-performance-kit"
out="$HOME/rk248-J4-feature-off"
[ ! -e "$out" ]
sbatch "$perf/performance.sbatch" "$kit" "$out" J4 feature-off
```

For the J4 storage-only comparison, extract a **fresh copy** of the source kit under another name before submission: the first campaign creates build files inside its source tree.

```bash
[ ! -e "$HOME/bbserv-candidate2-mirrored" ]
tar -C "$HOME" --transform='s/^bbserv-candidate2/bbserv-candidate2-mirrored/' -xzf "$HOME/bbserv-candidate2.tar.gz"
kit="$HOME/bbserv-candidate2-mirrored"
out="$HOME/rk248-J4-mirrored"
[ ! -e "$out" ]
sbatch "$perf/performance.sbatch" "$kit" "$out" J4 mirrored
```

For J5, use another fresh extraction and the historical fixture:

```bash
[ ! -e "$HOME/bbserv-candidate2-J5" ]
tar -C "$HOME" --transform='s/^bbserv-candidate2/bbserv-candidate2-J5/' -xzf "$HOME/bbserv-candidate2.tar.gz"
kit="$HOME/bbserv-candidate2-J5"
out="$HOME/rk248-J5-feature-off"
fixture="$HOME/stage187-fixture-20260911-145307"
[ ! -e "$out" ]
sbatch "$perf/performance.sbatch" "$kit" "$out" J5 feature-off "$fixture"
```

GNU tar on bbserv supports the name transform. Do not start another campaign inside an already built kit. Each job's `rk248-JOBID.log` appears in the directory from which `sbatch` was called. `results.json` keeps the active case, source and binary hashes, placement evidence, science gates, run medians and decision; `status=failed` retains evidence when a gate fails. A single run's times are never treated as independent repeated measurements.

## 3. Collect the result after the job ends

On bbserv, inspect `squeue -j JOBID` and the Slurm log. After the job exits, create and verify a small review archive:

```bash
python3 "$perf/collect.py" "$out" "$HOME/rk248-J4-feature-off-review.tar.gz"
python3 "$perf/collect.py" --verify "$HOME/rk248-J4-feature-off-review.tar.gz"
sha256sum "$HOME/rk248-J4-feature-off-review.tar.gz"
```

Change the archive name to match each campaign. The archive includes logs, placement and resource records, hashes, numerical comparisons and the final decision. It also includes the Slurm log when available and can collect an interrupted job after it leaves the queue. It omits the large checkpoints, VTK payloads and binaries; keep the full results directory on bbserv until review. Rank peak RSS is reported per process; their sum is not a simultaneous node-memory measurement. On the local machine, retrieve and verify the review archive:

```bash
scp kevlahan@bbserv.math.mcmaster.ca:~/rk248-J4-feature-off-review.tar.gz ~/Desktop/
python3 /Users/kevlahan/.codex/worktrees/7173/wavetrisk_hydrostatic/test/parallel_block_profile/performance_248/collect.py \
    --verify ~/Desktop/rk248-J4-feature-off-review.tar.gz
```

A completed scientific gate is a prerequisite for interpreting any timing result. J5 has not yet been validated with the new solver on bbserv; its correctness runs are the first gate of that campaign.
