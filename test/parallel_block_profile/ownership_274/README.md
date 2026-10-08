# Stage274 ownership-contract diagnostic

This begins the direct block producer design work on accepted Stage266. It does
not include Stage273's unqualified partition change. No numerical formula or
production ownership changes. It exports actual velocity source/restriction/
gradient plans and compiles a prospective ownership and data-movement budget.

The diagnostic tests whether moving velocity alone would introduce bulk input
transfers or same-level remote dependencies. Unmapped boundary/scaffold inputs
remain explicit. It helps choose a complete producer/consumer slice before
implementing it. It does not measure or promise a performance gain.

Four serial builds and eight full J5 runs on exclusive nk/bb13 retain main and
Stage266 science anchors, checked/oracle tests and 83/160/256-rank plan samples.
Estimated execution: **12–18 minutes after allocation**, with extra uncertainty
from the new export/postprocessing. Queue and transfer time are additional.
There are no timing comparisons, automatic repeats or automatic promotion.

On the local machine:

```bash
scp /Users/kevlahan/.codex/worktrees/26d8/wavetrisk_hydrostatic/.stage274-validation/rk274-ownership-contract-kit-v1.tar.gz kevlahan@bbserv.math.mcmaster.ca:~/
```

On bbserv:

```bash
tar -xzf "$HOME/rk274-ownership-contract-kit-v1.tar.gz" -C "$HOME"
bash "$HOME/rk274-ownership-contract-kit-v1/launch.sh"
```

Check progress:

```bash
bash "$HOME/rk274-ownership-contract-kit-v1/status.sh"
```

After completion, collect and transfer the review to the local Desktop:

```bash
bash "$HOME/rk274-ownership-contract-kit-v1/collect.sh"
```

The kit requires only the established hash-checked J5 fixture, not earlier kits.
The first three actual plan builds per Domain are sampled. These are structural
samples, not a synchronized census or a complete time-weighted RK trajectory.
Full kernel argument packing gives a conservative protocol budget. Hypothetical
upstream block ownership is labeled separately from retained Domain producers.
All source dependency blockers and adapter work remain visible in the report.
The source compiler and reviewer reject missing/duplicate/corrupted plan records.
