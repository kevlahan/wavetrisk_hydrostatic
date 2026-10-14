# Stage 278: runnable attribution kit v1

Prepared 2026-10-13. The kit is sealed and ready for user transfer/submission.
No source was uploaded, no remote command was run, and no Slurm job was submitted.
Stage 266 remains accepted. This diagnostic campaign does not qualify a repair
or resolve the 256-rank scaling failure by itself.

## Deliverables and checks

- [Archive](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/rk278-attribution-kit-v1.tar.gz) (23,851,518 bytes).
- [SHA256 sidecar](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/rk278-attribution-kit-v1.tar.gz.sha256).
- [Full kit instructions](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/rk278-attribution-kit-v1/README.md), [schedule](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/rk278-attribution-kit-v1/campaign-spec.json).
- [Local qualification](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/rk278-attribution-kit-v1/LOCAL_VALIDATION.json), [26 helper tests](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/rk278-attribution-kit-v1/TEST_VALIDATION.json),
  [delivery verification](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/attribution-kit-v1-delivery.json).

Archive SHA256: `85e4960e80a82bb3373621e1803cfea4f7f0ad308725fec5b219decd1011fd23`.

The only source change from sealed v5 is the six-file attribution overlay.
The original v5 binary is built separately as a control; the canonical-fallback
repair is excluded. Packaged attribution source exactly matches its independent
local optimized off/coarse/fine, optimized-oracle and checked-fine replays.
All twelve states, exact accepted-family checkpoint science, 62 VTK members and
reload/remap checks pass. Five source/build families are represented by four
source trees, with the checked build sharing the attribution candidate source.
The delivered source includes 24 changed files relative to accepted 266.

The 26 helper tests exercise both complete 31/37-run schedules using mocked Slurm,
failed-run persistence and collection, duplicate submission guards, source-mutation
rejection, mode isolation, fine-run admission and raw-profile validation. These
workflow tests do not establish cluster speed or cluster science. The real local
solver tests use at most four MPI ranks with serial builds/replays. Exact packaged
source identities and analysis imports were verified from a fresh archive
extraction. Fresh GCC 13/MVAPICH compilation remains a first step inside the job;
the proven physics dependency-include fix is retained.

## Bounded campaign and decision

Use exclusive nk/bb13, 256 allocated tasks and all node memory, with the same
physical-core/NUMA placement policy as job 55502. Build main, 266, original v5,
attribution-only v5 optimized and checked serially, sharing the 23 external
physics objects. Keep the established GCC 13.2.0/MVAPICH 4 environment and J5 fixture.

There are eight separate science/ownership validation runs. At each of 83, 160
and 256 ranks, measure original-v5 off, attribution off, coarse, coarse, off,
original-v5 off in that order. Add matched main/266 controls (266 only at 256).
This gives 31 required runs. Record ordinary complete-step timing including
remaps, separate remap/non-remap summaries, and per-rank child-process memory.
Timing runs contain neither independent numerical oracles nor ownership audit
output. The repeated v5/off/coarse groups measure run spread, disabled-hook effects
and observed profiler perturbation separately.

Up to three additional fine/off pairs are admitted by the recorded coarse
results, making 37 runs maximum. Admission requires enough exchange time to cover
the positive 266 gap above the observed run-spread screen, or exchange growth
accompanying the 160-to-256 slowdown above that screen. This is a diagnostic
selection rule, not a performance threshold or proof of transport dominance.
There is no automatic second job or open-ended repeat loop.

Plan roughly 30–45 minutes of allocated runtime, with a 90-minute wall limit.
This is an estimate based on job 55502's actual 19:55 for four builds and 16 runs,
not a measured duration for this larger schedule. Queue and collection are extra.
Real start/end/elapsed accounting is collected for the runtime ledger after review.

Review must identify the cost budget for the next substantial repair: preparation,
compute/copy work, retained inverse/wavelet/compatibility paths, imbalance or
exchange. MPI waits can include late producers and do not establish network cost.
If attribution remains ambiguous, preserve that uncertainty rather than prescribe
batching. Absolute performance, paired Domain-path retirement, bounded memory
and useful scaling remain separate scorecard gates. 256 remains an intermediate
test; the >=640-rank goal is not evaluated here.

## Transfer from the Mac

```bash
scp "/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/rk278-attribution-kit-v1.tar.gz" "/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/rk278-attribution-kit-v1.tar.gz.sha256" kevlahan@bbserv.math.mcmaster.ca:~/
```

## Launch on bbserv

Use the same loaded compiler/MPI environment as job 55502. The existing fixture
is `~/stage187-fixture-20260911-145307`; no fixture edits are needed.

```bash
set -e
cd "$HOME"
sha256sum -c rk278-attribution-kit-v1.tar.gz.sha256
test ! -e rk278-attribution-kit-v1
tar -xzf rk278-attribution-kit-v1.tar.gz
bash "$HOME/rk278-attribution-kit-v1/launch.sh"
```

## Status and collection on bbserv

```bash
bash "$HOME/rk278-attribution-kit-v1/status.sh"
```

After the job ends, including after failure:

```bash
bash "$HOME/rk278-attribution-kit-v1/collect.sh"
```

Collection preserves failed build/run evidence, source and binary hashes,
profile rows, placement and memory records, numerical reports and actual Slurm
accounting. It verifies the review archive and copies it to the established
murier Desktop destination using the saved job number. Binaries and checkpoint/
VTK payloads are omitted; their cluster-computed hashes and numerical reports are
retained. The independent reviewer cannot rerun omitted payloads.

The existing latest-job guard requires collection of the preceding Stage 278 job
before another submission. Status and collection verify that the saved job belongs
to this kit. The archive, predecessor kit and previous campaign evidence are
preserved; no new performance acceptance claim accompanies this delivery.
