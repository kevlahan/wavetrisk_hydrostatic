This tests one Stage266 change: defer the first full block installation on
provisional Stage203 production inverse exits until after the existing all-level
boundary refresh. Accepted exits, non-Stage203 execution and independent oracle
runs retain their earlier installation. Canonical tiled publication, finite and
extent preflight, boundary epochs, scatter and the final installation remain.

The kit bundles exact accepted Stage265t commit fab72486, legacy main with its
established checkpoint I/O hardening, and the new candidate. It requires only
`$HOME/stage187-fixture-20260911-145307` at runtime on bbserv.

One exclusive nk/bb13 allocation builds four solvers serially with the same 23
physics objects. Seven correctness runs cover main/baseline and candidate
optimized/checked/oracle at 16 ranks, then candidate production/oracle at 167
ranks. All existing checkpoint, VTK, restart/remap, storage, physical placement
and zero-Domain-rank checks remain. Each oracle proves all 48 inverse returns
and verifies the original installation counts (84 calls across 48 exits).
Production uses one installation per exit, enforced by an actual-call counter
on every rank. It adds no per-transaction proof logging to timing runs.

The initial timing batch has six unprofiled 83-rank runs:
`main, Stage265t, candidate, candidate, Stage265t, main`.
There are thirteen solver runs total. Compare each run's precise ordinary-step
median, excluding steps 1, 8 and 9. No automatic acceptance or promotion occurs.
A repeatable candidate/baseline regression blocks the next removal; a mixed
result remains inconclusive. Main remains the anchor for the inherited gap.

On murier:

```bash
scp /Users/kevlahan/.codex/worktrees/26d8/wavetrisk_hydrostatic/.stage266-validation/rk266-deferred-installation-kit-v1.tar.gz kevlahan@bbserv.math.mcmaster.ca:~/
```

On bbserv with the usual compiler/MPI modules:

```bash
tar -xzf "$HOME/rk266-deferred-installation-kit-v1.tar.gz" -C "$HOME"
bash "$HOME/rk266-deferred-installation-kit-v1/launch.sh"
```

Check status:

```bash
bash "$HOME/rk266-deferred-installation-kit-v1/status.sh"
```

After completion or failure, collect the review to murier Desktop:

```bash
bash "$HOME/rk266-deferred-installation-kit-v1/collect.sh"
```

The launcher rejects duplicates until the prior review is collected. Status
handles purged job IDs; collection includes early failure diagnostics. Results
remain at the path in `$HOME/rk266-latest.json`.

Only if review recommends one repeat for small/noisy/mixed evidence:

```bash
bash "$HOME/rk266-deferred-installation-kit-v1/repeat.sh"
```

This adds exactly six timing runs, reusing verified binaries, source seals,
physics objects, toolchain and physical placement without rebuilds or repeating
validation. Every additional timing run still passes science gates. The initial
review remains immutable and a second repeat is rejected. Collect with the same
collect.sh command. No repeat is automatic.

Before delivery, require real exact optimized, checked and oracle J4 replays
with no more than four local MPI ranks, serial builds and one library thread.
The independent reviewer recomputes timing, installation and scientific gates.
It verifies cluster-computed checkpoint/VTK hashes and numerical reports;
omitted binary/checkpoint/VTK payloads are not rerun locally.
