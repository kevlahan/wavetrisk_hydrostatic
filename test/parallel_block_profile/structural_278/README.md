# Stage 278 structural repair campaign — kit v1

Prepared 2026-10-14. Local delivery only: no transfer or submission has occurred.
Candidate is the final locally qualified v7 structural metadata/hash repair.
Previous is the exact attribution candidate measured by job 55503, with profiling
normally disabled. Stage266 remains accepted; this kit does not promote v7.

## Question and source boundary

Does removing expensive one-use structural plan construction and hash clustering
reduce ordinary complete-step time, close the 266 gap and improve 160→256 scaling,
without restoring retired Domain numerical paths or increasing memory costs?

Four sealed source trees produce five serial builds: accepted266, hardened legacy
main, exact55503 previous, v7 optimized and v7 checked. Candidate differs from
previous in exactly six files. Its 251-file source identity matches local v7
qualification; the full accepted266 overlay is 24 files. All controls are unchanged.
The established GCC13 physics dependency include repair and shared 23-object
Simple-physics build are retained. Do not merge another working checkout into this kit.

V7 substitutes one current-layout geometry directory and indexed, deduplicated
metadata pulls for sixteen general one-use plans. It fixes rank/slot hash
clustering in indexed and existing route deduplication. It preserves numerical
placement, request order, stage exchanges and existing Domain retirements. Layout
and epoch validity remain checked; no speculative cross-epoch cache is included.
No new communicator or persistent metadata cache is retained.

Local evidence: twelve complete optimized/checked/oracle replays pass; component
suites cover changing ownership/layout/epochs, duplicate requests, zero-geometry
consumers, errors and exchange stress. Coarse route/input/retention construction
costs fell about 50%/61%/90% locally. Bare run medians improved only 3.1%, below
15–21% observed repeat spreads; ordinary means improved 16.2% but are sensitive to
run variation. Neither statistic is a guaranteed cluster speedup.

## Fixed schedule and timing policy

One exclusive nk/bb13 allocation, all node memory, 256 allocated ranks, the proven
physical-core/NUMA placement policy. Five make -j1 builds and 48 serial solver runs:

- Eight separate science/ownership validations, including checked at16 and
  independent oracles at16 and256. No timing acceptance from these runs.
- At83 and160, main then266 controls; at256,266 only.
- At each count: previous, candidate, candidate, previous, candidate, previous,
  previous, candidate, all profiling/oracles/ownership audits disabled.
- Repeat the main/266 controls in reverse order after that group.
- One previous/candidate coarse-profile pair per count (order reversed at160).

This is eight validations, 24 repair/control bare runs, ten absolute-control bare
runs and six coarse runs. No fine runs, adaptive additions or automatic second job.
Four bare repetitions/version/count address the local timing uncertainty within
one predeclared job. Two main/266 repetitions bracket each group.

Primary report: mean of run-level ordinary medians, with all four individual
repair/control observations and run spread. Secondary: ordinary sums/means,
remap/non-remap medians and means. Exclude established indices0,7,8 for startup,
reload and saves; retain ordinary remaps4,9. Nine ordinary steps/run remain.
Steps are not independent repetitions. Short-fixture save frequency is for
lifecycle science, not a production performance model; rare production saves
must not outweigh ordinary-step performance. No ad-hoc frequency weighting.

Record all rank wait4 RSS/CPU/fault data, placement, sources/binaries, environments,
raw profiles and actual Slurm accounting. Sum of rank RSS peaks is not a concurrent
job peak. Coarse profiles include profiler perturbation and are a single pair.
New indexed metadata pulls sit in preparation parents whereas old generic pulls
also appeared in exchange region78: a fall in region78 alone is not transport
speedup. Never add inclusive profile regions as if disjoint.

## Review and scorecard

All established science gates must pass: exact accepted-family fields/topology,
main-relative tolerance, twelve states, two reloads/remaps and62 VTK members;
actual numerical work on zero-Domain ranks and existing retirement guards.
Review raw construction budgets against the six-file change, then bare full-step
medians and means, 266/main absolute gaps, repeat spread, remaps, RSS and scaling.
A construction gain without a repeatable complete-step benefit is insufficient.
An absolute-time gain with256 slower than160 does not fix scaling.
Unexpected regressions require attribution before another migration/campaign.
There is no invented numerical acceptance threshold or automatic promotion.
256 is an intermediate test; >=640 matched multi-node performance, whole-timestep
independence and legacy parity remain open. This repair preserves the existing
ownership migration while addressing its replacement cost; it retires no further
numerical Domain readers. Shared-memory transport is outside this campaign.

Budget 45–60 minutes allocated runtime, 90-minute wall limit; queue/collection
extra. Replaying the prior observed build/launch costs against this schedule gives
about46.5 minutes, not a measurement or guaranteed duration. Actual job55503 took
38:06 for five builds/37 runs. New Slurm start/end/elapsed remain pending.

## Run on bbserv

Use the established GCC13.2.0/MVAPICH4 environment and unmodified fixture
~/stage187-fixture-20260911-145307. The Mac archive must first be transferred.

```bash
scp "/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/rk278-structural-kit-v1.tar.gz" "/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/rk278-structural-kit-v1.tar.gz.sha256" kevlahan@bbserv.math.mcmaster.ca:~/
ssh kevlahan@bbserv.math.mcmaster.ca
```

On bbserv:

```bash
set -e
cd "$HOME"
sha256sum -c rk278-structural-kit-v1.tar.gz.sha256
test ! -e rk278-structural-kit-v1
tar -xzf rk278-structural-kit-v1.tar.gz
bash "$HOME/rk278-structural-kit-v1/launch.sh"
```

```bash
bash "$HOME/rk278-structural-kit-v1/status.sh"
```

After completion or failure:

```bash
bash "$HOME/rk278-structural-kit-v1/collect.sh"
```

Collection retains failed-job evidence, verifies the returned archive and copies
rk278-review-JOB.tar.gz to the established murier Desktop using the saved job ID.
The command needs no editable job number. Status/collection reject a latest-job
record belonging to another kit. Previous job collection is required before a new
launch. Binaries/checkpoint/VTK payloads are omitted from the return archive;
recorded hashes and numerical comparisons are retained, not a local rerun of them.
