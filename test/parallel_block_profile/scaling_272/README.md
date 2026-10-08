# Stage272 ownership and scaling diagnostic

Use the accepted Stage266 solver and a separate observer copy to explain the
160-to-256-rank slowdown. Production sources are unchanged. The kit is
self-contained apart from the existing, hash-checked J5 fixture on bbserv.
No Stage270/271 kit or results directory is needed there. In filenames and
runner logs, `candidate` means the optimized observer copy; `checked` means
the same observer source built with runtime checks.

One exclusive `nk`/bb13 allocation makes four serial builds (Stage266, legacy
main for correctness only, optimized observer, checked observer) and fourteen
solver runs. Budget **16–22 minutes after allocation**, based on Stage271's
21:47 for four builds/twenty runs and Stage270's 15:38 for five builds/fifteen
runs. Queue and transfer time are additional; the 256-rank oracle and filesystem
load can change this estimate. Actual Slurm and runner elapsed times are saved.

On the local machine:

```bash
scp /Users/kevlahan/.codex/worktrees/26d8/wavetrisk_hydrostatic/.stage272-validation/rk272-ownership-scaling-kit-v1.tar.gz kevlahan@bbserv.math.mcmaster.ca:~/
```

On bbserv:

```bash
tar -xzf "$HOME/rk272-ownership-scaling-kit-v1.tar.gz" -C "$HOME"
bash "$HOME/rk272-ownership-scaling-kit-v1/launch.sh"
```

Check progress on bbserv:

```bash
bash "$HOME/rk272-ownership-scaling-kit-v1/status.sh"
```

After the job ends, collect on bbserv. This also copies the review to the local
Desktop and works for failed jobs or jobs purged from the active queue:

```bash
bash "$HOME/rk272-ownership-scaling-kit-v1/collect.sh"
```

The launcher prints the job ID/log path, uses a fresh result directory, and
blocks accidental resubmission until the previous review has been collected.
The review name is `rk272-review-JOBID.tar.gz`; no filename editing is needed
for the commands above. Partition `nk` contains only bb13; no exclusion is used.

## What the evidence will distinguish

- Domain-owned producer work, including basic-operator patch-layer visits and
  primitive preparation, versus block restriction work and actual block/patch
  inventories. Groups use recorded ownership, never assumptions about rank IDs.
- Capture publication, restriction, inverse and migration communication scopes;
  selected remote logical payload counts are checked for send/receive conservation.
- Preparation, computation and MPI time distributions for Domain owners and
  ranks without Domains, separately for each independent observer run.
- Whole-step observer overhead at each rank count against an untouched build.
  Failed calibration remains visible and prevents qualified timing attribution.

The schedule uses six scientific validations, then these timing runs:
`160 bare, 160 observer, 256 bare, 256 observer, 256 observer, 256 bare,
160 observer, 160 bare`. Each count has its own A-B-B-A comparison. This revisits
160 ranks after 256; it reduces but does not eliminate run-order effects.
The same twelve-step J5 trajectory is used; ordinary timing excludes steps
1, 8 and 9. Checked/oracle runs never enter timing comparisons.

The diagnostic reuses the bounded Stage270 phase context and adds ownership,
producer work, restriction work and capture volumes. Local producer patch
batches and native-mass/Domain-boundary waits have separate scopes. It adds no MPI collective,
barrier, timed-step file output, or timer inside per-node physics callbacks. Records
are flushed per rank after the solver run. Basic-operator timing is at the
existing level/layer transaction, with integer work updates once per Domain.

MPI scope time includes both transfers and waiting for other ranks. CPU time
can include MPI spinning. Local durations are unsynchronized and do not prove
a global critical path. Inclusive phase times overlap. Inventories are not
executed work. Hardware bandwidth, NUMA and clock effects are not excluded.
The campaign does not promise a speedup, accept an implementation automatically,
or meet the eventual target of at least 640 MPI ranks for J5's 160 Domains.
