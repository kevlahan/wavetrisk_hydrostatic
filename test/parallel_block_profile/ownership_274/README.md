# Stage274 ownership-contract diagnostic: v2 recovery

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

Version 2 fixes the reader's incorrect assumption that Domain sample numbers
were global. They are local to each observing rank. It preserves observations
on multiple ranks, reports unique Domain counts separately from observation
counts, and retains corruption checks within each rank's sequence. No solver
source or numerical binary changes.

The commands below recover the reported v1 failure. All four sealed builds and
three completed solver runs are reused after verification; all scientific gates
are rerun on the copied outputs. Only the five remaining solver cases execute.
The original result directory is untouched. Both job attempts' accounting is
retained. Allow roughly **6–12 minutes after allocation**; recovery copy and
postprocessing costs make this estimate approximate. Do not resubmit the v1 kit.

On the local machine:

```bash
scp /Users/kevlahan/.codex/worktrees/26d8/wavetrisk_hydrostatic/.stage274-validation/rk274-ownership-contract-kit-v2.tar.gz kevlahan@bbserv.math.mcmaster.ca:~/
```

On bbserv:

```bash
tar -xzf "$HOME/rk274-ownership-contract-kit-v2.tar.gz" -C "$HOME"
bash "$HOME/rk274-ownership-contract-kit-v2/recover.sh"
```

Check progress:

```bash
bash "$HOME/rk274-ownership-contract-kit-v2/status.sh"
```

After completion, collect and transfer the review to the local Desktop:

```bash
bash "$HOME/rk274-ownership-contract-kit-v2/collect.sh"
```

The kit requires only the established hash-checked J5 fixture, not earlier kits.
The first three actual plan builds per Domain per observing rank are sampled. These are structural
samples, not a synchronized census or a complete time-weighted RK trajectory.
Full kernel argument packing gives a conservative protocol budget. Hypothetical
upstream block ownership is labeled separately from retained Domain producers.
All source dependency blockers and adapter work remain visible in the report.
The source compiler and reviewer reject missing/duplicate/corrupted plan records.

A fresh campaign can use `launch.sh` instead of `recover.sh`; it performs all
four builds and eight runs. Recovery deliberately accepts only the sealed v1
failure after the first three successful solver exits and refuses an active or
already-submitted recovery. It validates source and binary hashes before reuse.
Repeated observations across ranks are not assumed simultaneous: their totals
are observation budgets, not a census of distinct mesh work at one instant.
