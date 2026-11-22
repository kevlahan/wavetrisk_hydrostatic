# RK migration: qualified bbserv candidate

## Resuming on another computer

Continue Stage 247 on `codex/performance-investigation`. Before editing, building or submitting a job, run `git status --short --branch` and `git remote -v`. Preserve local work. Set `mirror` below to the local name of the WAVETRISK GitHub or Bitbucket remote.

```sh
mirror=github
branch=codex/performance-investigation
git fetch "$mirror" "refs/heads/$branch"
git switch "$branch"
git merge --ff-only FETCH_HEAD
git merge-base --is-ancestor a4e821a07ce757709a4e4cc329957f602724ea31 HEAD
```

If the branch has no local checkout yet, use `git switch -c "$branch" FETCH_HEAD` for the switch step. Run the commands in order and stop if any fails. The final check exits successfully for the repair commit or a later descendant. If local work and the published branch have diverged, inspect both histories before proceeding.

Both mirrors received the checkpoint and restart repairs. Use the ancestor for the branch being resumed.

| Branch | Required repair ancestor | Role |
|---|---|---|
| `codex/performance-investigation` | `a4e821a0` | Active Stage 247 refactor |
| `codex/stage187-validated` | `bfdb68b8` | Stage 205 handoff |
| `parallel-blocks` | `441b6f64` | Earlier blocks baseline |

The sealed Stage 247 kit predates the checkpoint and restart repairs. Its manifest and validation results apply to the frozen source. Build from the verified checkout when testing the repaired code, and record the commit used for any new bbserv measurements.

## Stage 247 handoff

The implementation is complete for the restarted, compressible, unsplit climate RK4 case with Simple physics. Stage 247 closes the lifecycle, coefficient handoff and persistent legacy storage work listed in the Stage 244 handoff. Performance remains to be measured on bbserv.

Production stepping now retains native RK, physics, soil, adaptation and remap state. Vector coefficients pass directly from their native producer to the native transform workspace. Legacy Domain SOL/WAV numerical arrays are released outside initialization and bounded output/checkpoint adapters. Their per-field address capacities remain as metadata so adaptation can extend the topology without recreating numerical mirrors. Obsolete RK snapshots and persistent scalar writeback images are also absent.

Active vector restriction/wavelet scratch, native stencil and alias workspaces, fixed reference fields and communication buffers remain necessary production storage. This change does not claim to remove those buffers. Existing off and independent-oracle paths remain available; the candidate requires explicit build selectors and rejects incompatible runtime configurations.

The frozen candidate is `.stage247-validation/candidate6`; cases are under `.stage247-validation/final6`. The [source-bound Stage 247 manifest](parallel-blocks-stage247-results.json) records eight passing cases: a full bounds-checked 4-rank trajectory, optimized full trajectories on 4 and 41 ranks, an independent RK/physics oracle, an off-path control, and three deliberate failures for missing lifecycle state, missing coefficient handoff and forbidden persistent Domain storage. Each full trajectory has 12 exact reference states, two remaps, 62 exact VTK members and the existing checkpoint numerical gate. The independent RK oracle has 16 exact alias comparisons. Fifteen migrated transform kernels retain the audited legacy formulas. No numerical tolerance was relaxed.

All local builds used `make -j1`; solver runs were sequential with library threads capped at one. The oversubscribed 41-rank case verifies empty-geometry-rank behavior and communication correctness only. It provides no performance or scaling evidence.

The [sealed bbserv kit manifest](parallel-blocks-stage247-bbserv-kit.json) identifies the transferable archive and its SHA-256 hash. The kit contains the qualified sources, J4 restart fixture, exact build selectors, runtime environment, local validation manifest and file hashes. Its archive is unpacked and checked locally. Remote builds, allocations, rank placement, baseline comparisons and timing repetitions are the next step; avoid bb04. No bbserv job has been started.
