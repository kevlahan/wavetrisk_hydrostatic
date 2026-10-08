# Stage273 balanced contiguous partition candidate

This isolated candidate changes how blocks are assigned to ranks. It minimizes
the maximum contiguous load under the existing leaf-weight model, then chooses
feasible cuts near cumulative equal-work targets. Catalogue order, block
subdivision, weights, numerical kernels and payload proof checks are retained.

The evidence supports testing this change, not promising a solver speedup.
Stage272 found that 96 of 256 ranks received only 22.79% of block patch/divergence
work and none of the Domain producer work. The original assignment rule can
front-load early ranks even with equal weights. The new rule guarantees an
optimal peak for the existing contiguous weight model; actual runtime also
includes communication, fixed per-block costs and retained Domain work.

One exclusive nk/bb13 job makes four serial builds and twenty runs. Expect
**20–25 minutes after allocation**, based on Stage271's 21:47 for the same
build/run counts. Validation-only catalogue output may add some I/O time.
Queue and transfer time are additional. Actual Slurm and runner elapsed time
are collected automatically. There is no automatic repeat or promotion.

On the local machine:

```bash
scp /Users/kevlahan/.codex/worktrees/26d8/wavetrisk_hydrostatic/.stage273-validation/rk273-balanced-partition-kit-v1.tar.gz kevlahan@bbserv.math.mcmaster.ca:~/
```

On bbserv:

```bash
tar -xzf "$HOME/rk273-balanced-partition-kit-v1.tar.gz" -C "$HOME"
bash "$HOME/rk273-balanced-partition-kit-v1/launch.sh"
```

Check progress:

```bash
bash "$HOME/rk273-balanced-partition-kit-v1/status.sh"
```

After completion, collect and transfer the review to the local Desktop:

```bash
bash "$HOME/rk273-balanced-partition-kit-v1/collect.sh"
```

The kit needs only the established hash-checked J5 fixture on bbserv, not any
older campaign kit. Launch uses a fresh result directory and blocks accidental
duplicates. Collection also works for failed/purged jobs. No bb04 exclusion is
needed on nk.

The four builds are untouched Stage266, legacy main with established checkpoint
I/O hardening, optimized candidate, and checked candidate, sharing 23 physics
objects. Six validations cover main/baseline/candidate/checked/oracle at 16 ranks
and candidate oracle at 256. Fourteen timing runs cover main plus baseline–candidate–
candidate–baseline at 83 and 160 ranks, then the same block pair order at 256.
All use the established J5 trajectory; exclude steps 1, 8 and 9 from timings.

Candidate validations print the actual catalogue weights, roots and geometric/
block owners. Independent review verifies candidate cuts against a dynamic
programming reference and reconstructs baseline assignment on those same weights.
This offline baseline assignment is a counterfactual, not a baseline timing
measurement. Every timing run disables detailed catalogue output and has no
added phase profiler. Complete-step performance determines acceptance.
