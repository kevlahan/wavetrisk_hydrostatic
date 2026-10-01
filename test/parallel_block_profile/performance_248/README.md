# RK migration performance campaign on bbserv

This kit runs the qualified Stage247 source on the exclusive `nk` node (bb13, 256 effective CPUs). It builds the candidate and matched baseline there with `make -j1`, then checks correctness and independent timing runs. The kit's launcher follows the successful Stage227/229 pattern: it extracts a fresh source tree, submits with `sbatch --parsable`, records the job ID and prints the log and results paths. The collector verifies a review archive and transfers it from bbserv to murier. If the job fails before `results.json` exists, the collector still preserves the Slurm log and submission record.

The source archive SHA-256 is `dc6daba22953f5ab54398f6085b9b39aa04a2161b233fade692fedaa9d5d7d63`. The performance archive SHA-256 is recorded in `docs/parallel-blocks-stage248-performance-kit.json`. The source archive `bbserv-candidate2.tar.gz` must already be in `$HOME` on bbserv; its old extracted directory is not used.

## 1. Transfer from the local machine

```bash
scp /Users/kevlahan/.codex/worktrees/7173/wavetrisk_hydrostatic/.stage248-validation/rk248-performance-kit-v3.tar.gz kevlahan@bbserv.math.mcmaster.ca:~/
```

## 2. Extract and submit on bbserv

```bash
(
set -e
cd "$HOME"
if [ ! -d rk248-performance-kit-v3 ]; then tar -xzf rk248-performance-kit-v3.tar.gz; fi
bash "$HOME/rk248-performance-kit-v3/launch.sh"
)
```

The launcher defaults to J4 feature-off. It prints `Submitted job`, `Results`, and `Slurm log`. It creates a fresh source tree from the qualified archive for every submission; it never silently skips submission because an old results directory exists. To see the latest job's status and log:

```bash
source "$HOME/rk248-latest.state"
squeue -j "$job"
sacct -j "$job" -X -o JobID,JobName%25,State,ExitCode,Elapsed,NodeList
tail -n 80 "$HOME/rk248-${job}.log"
```

An empty queue means only that the job is no longer pending or running; inspect `sacct` and the log. A log file may not exist while a job is still pending.

## 3. Collect after the job ends on bbserv

```bash
bash "$HOME/rk248-performance-kit-v3/collect.sh"
```

The collector checks the queue, creates or reuses `$HOME/rk248-review-JOBID.tar.gz`, verifies its file hashes, and copies it to `kevlahan@murier.math.mcmaster.ca:~/Desktop/`. If that transfer fails, the verified archive remains on bbserv; pull it from the local machine with `scp` using the printed archive path. Keep the full results directory on bbserv until review. The review archive omits large checkpoints, VTK files and binaries.

For J4 storage-only or J5 feature-off, use `bash "$HOME/rk248-performance-kit-v3/launch.sh" J4 mirrored` or `bash "$HOME/rk248-performance-kit-v3/launch.sh" J5 feature-off`. J5 requires `$HOME/stage187-fixture-20260911-145307`.

The candidate enables the qualified Stage247 migration. `feature-off` compares it to the same source with the new RK, physics, soil, adaptation, lifecycle, coefficient and mirror-free features disabled. `mirrored` retains the native migration in both builds and isolates persistent Domain SOL/WAV storage. Every full run completes 12 steps, two remaps, checkpoint reload, exact VTK and numerical gates before timing. Timing uses the median of nine ordinary steps per independent run, followed by a bounded A-B-B-A confirmation when screening permits. No branch promotion occurs automatically. The local four-core machine is used only for preparation; cluster timing is measured on bb13.
