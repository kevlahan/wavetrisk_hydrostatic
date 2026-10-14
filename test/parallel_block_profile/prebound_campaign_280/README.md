# Stage280 direct-provider and compact-publication campaign — kit v1

Prepared 2026-10-14. Ready for user transfer/submission; no campaign upload or submission has occurred. Stage266 remains accepted. The exact five source families are candidate (repair), Stage279 job55507, previous (v7 job55504), Stage266 and main with the established checkpoint-I/O hardening. The checked executable uses the candidate source. Do not overlay a working checkout.

## Repair and evidence

The two-file repair preserves Stage279's independent numerical ownership of three provisional RK4 transforms and paired retirement of their displaced Domain work. Sixteen ordinary routes bind freshly proved unique canonical providers directly, removing128 redundant preparation collectives per state. Ordered outer-vector alias proofs and fixed-input writer selection remain. Local publication reads canonical rows; remote duplicates share one import row. Communicator and scratch lifetimes, invalidation and state guards are retained. Payload volumes and runtime route counts are unchanged. Final transform/physics/adaptation, block granularity and load assignment are outside this repair.

All11 full local solver replays and48 optimized/checked exchange cases passed. Exact-source J5 GCC13/MVAPICH checked and optimized compilation passed in an isolated compiler copy. Four-rank J4/J6 bare comparisons (two runs each) showed4.79% improvement versus v7 and5.38% versus Stage279 in mean ordinary run medians, and29.64% less named extra-transform allocation. These are local observations, not a cluster forecast or total-memory qualification. The local coarse comparison showed a small whole-step benefit and increased RHS waiting; retain that limitation. IMPLEMENTATION.md and QUALIFICATION_SUMMARY.json contain quantitative evidence and caveats.

## Bounded campaign

Use the established GCC13.2.0/MVAPICH4 environment and fixture ~/stage187-fixture-20260911-145307. One exclusive nk/bb13 allocation,256 CPUs, all node memory, verified physical-core/NUMA placement. Six serial builds and exactly60 serial solver runs:

- Eleven science validations: main/266/v7/279/candidate/checked/oracle at16, production candidate at83/160/256, and candidate oracle at256. Only candidate/checked validations have ownership audits.
- Forty bare timing runs: four v7 and four candidate repetitions at83/160/256 in ABBABAAB order. Two main/266/279 controls bracket each83/160 group; two266/279 controls bracket256. Main cannot run above160 Domains.
- Nine separate coarse runs: one v7/279/candidate comparison at each count, alternating order at160. No fine profiling or adaptive runs.

Primary score: mean of independent bare-run medians over ordinary step indices1-6 and9-11, excluding startup/reload/checkpoint indices0,7,8, retaining remaps4,9. Report ordinary means/totals, remap/non-remap groups, independent-run observations/spread and profiles separately. Checkpoint-heavy validation is excluded from timing. Production checkpoints are rare; ordinary steps decide performance.

The revised scorecard requires both ownership/retirement and absolute performance/scaling. Preserve exact accepted-family science, main-relative tolerance,12 states,2 reloads,2 remaps,62 VTK members and global class-by-class kernel census. At256, require useful transform/RHS work on96 zero-Domain ranks. Audit compact-publication row counts for every epoch/rank. Compare rank RSS/CPU/faults and named route storage; individual peak sums are not a simultaneous allocation peak. MPI wait includes late producers and is not pure transport time. Descriptor counts are not calibrated work weights. Region85 separates transform plan construction;85-90 separate new transform costs. The Stage279 control has these regions; v7 does not.

Success requires repeatable absolute ordinary-time improvement with science/ownership preserved, separately acceptable memory and improved160-to256 scaling. Compare v7 and279 to establish recovery,266 for accepted-source parity and main for the ultimate legacy gap. No automatic promotion. Final lifecycle, whole-memory parity, finer granularity/load assignment and matched J5>=640 multi-node scaling remain open. J6/J8 and shared-memory windows are outside this campaign.

Planning allowance60-80 minutes within a90-minute limit, excluding queue/collection. Prior job55507 took43:52 for five builds49 runs. This is an estimate, not a measured duration.

## Transfer from the Mac

```bash
scp /Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage280-validation/rk280-prebound-kit-v1.tar.gz /Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage280-validation/rk280-prebound-kit-v1.tar.gz.sha256 kevlahan@bbserv.math.mcmaster.ca:~/
ssh kevlahan@bbserv.math.mcmaster.ca
```

## On bbserv in the established compiler/MPI environment

```bash
set -e
cd "$HOME"
sha256sum -c rk280-prebound-kit-v1.tar.gz.sha256
test ! -e rk280-prebound-kit-v1
tar -xzf rk280-prebound-kit-v1.tar.gz
python3 "$HOME/rk280-prebound-kit-v1/check_kit.py"
bash "$HOME/rk280-prebound-kit-v1/launch.sh"
```

## Status

```bash
bash "$HOME/rk280-prebound-kit-v1/status.sh"
```

## Collect after the job ends, including failures

```bash
bash "$HOME/rk280-prebound-kit-v1/collect.sh"
```

Collection preserves diagnostics, independently verifies archive member hashes and copies rk280-review-JOB.tar.gz to murier Desktop. Attach that returned archive for review. The launcher prevents an uncollected Stage280 relaunch; do not bypass its evidence gate. No collection or deletion of Stage279 evidence is required for this fresh Stage280 kit.
