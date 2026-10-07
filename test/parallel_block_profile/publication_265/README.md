Stage265 tests direct canonical inverse SOL publication against exact Stage230
and legacy main. Both existing block installations remain. No other solver
changes are included. The kit bundles all three source trees and helpers and
requires only `$HOME/stage187-fixture-20260911-145307` on bbserv.

One exclusive nk/bb13 allocation builds four solvers serially, sharing the same
23 external physics objects. Seven correctness runs cover main/Stage230 at 16,
candidate optimized/checked/oracle at 16, and candidate production/oracle at
167 ranks, including ranks without geometric Domains. The independent oracle
compares all old-return destinations bitwise at each inverse exit, including
scaffold and soil. The production runtime checks the removed-return count.

Twelve unprofiled 83-rank runs form two batches:
`main, Stage230, candidate, candidate, Stage230, main`.
Each timing run also checks twelve states, both reloads/remaps, checkpoint
limits, 62 VTK members and exact block equivalence. Ordinary medians exclude
steps 1, 8 and 9. Candidate/Stage230 gaps use four order-balanced pairs. Main
remains a same-allocation anchor. No automatic promotion occurs.

On murier:

```bash
scp /Users/kevlahan/.codex/worktrees/26d8/wavetrisk_hydrostatic/.stage265-validation/rk265-direct-publication-kit-v1.tar.gz kevlahan@bbserv.math.mcmaster.ca:~/
```

On bbserv, with the usual compiler/MPI modules available:

```bash
tar -xzf "$HOME/rk265-direct-publication-kit-v1.tar.gz" -C "$HOME"
bash "$HOME/rk265-direct-publication-kit-v1/launch.sh"
```

Status:

```bash
bash "$HOME/rk265-direct-publication-kit-v1/status.sh"
```

After completion or failure:

```bash
bash "$HOME/rk265-direct-publication-kit-v1/collect.sh"
```

Collection transfers `rk265-review-JOBID.tar.gz` to murier Desktop. Submission
refuses duplicates until the previous review is collected. Status handles
purged job IDs; collection preserves early-failure diagnostics. Repeating
collection verifies and retransfers the existing review. Full results remain
at the path saved in `$HOME/rk265-latest.json`.

The transferred review omits binary/checkpoint/VTK payloads. Their cluster hashes
and numerical reports are reviewed, not rerun locally. The local reviewer also
verifies source/build identities, physical placement, storage, exact-return
proof coverage, raw states/timings, and recomputes the comparisons. A repeatable
regression blocks the next removal. Mixed/noisy evidence permits at most one
extra identical timing batch; no generic percentage allowance applies.
