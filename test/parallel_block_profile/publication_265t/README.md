This refines Stage265's direct inverse publication by copying complete patch
ranges through one reusable physical-patch tile. Canonical value selection,
finite/extent preflight, both block installations, boundary epochs and the exact
old-return oracle remain. No further Domain pathway is removed.

Stage230 remains the accepted performance baseline because the original
Stage265 candidate was inconclusive. The kit is self-contained, bundles main,
Stage230 and the refined candidate, and requires only the established bbserv
J5 fixture `$HOME/stage187-fixture-20260911-145307`.

One exclusive nk/bb13 allocation builds four solvers serially using the same 23
physics objects. Seven correctness runs cover main/Stage230 and candidate
optimized/checked/oracle at 16 ranks, plus candidate production/oracle at 167
ranks. These preserve all existing checkpoint, VTK, restart/remap, storage and
zero-Domain-rank checks; each oracle proves all 48 inverse returns bitwise.

The initial timing batch has only six unprofiled 83-rank runs:
`main, Stage230, candidate, candidate, Stage230, main`.
This gives two order-balanced candidate/Stage230 pairs and two main anchors.
There are thirteen solver runs total. Ordinary medians exclude steps 1, 8 and 9.
No automatic performance acceptance or promotion occurs.

On murier:

```bash
scp /Users/kevlahan/.codex/worktrees/26d8/wavetrisk_hydrostatic/.stage265t-validation/rk265t-tiled-publication-kit-v1.tar.gz kevlahan@bbserv.math.mcmaster.ca:~/
```

On bbserv, with the usual compiler/MPI modules available:

```bash
tar -xzf "$HOME/rk265t-tiled-publication-kit-v1.tar.gz" -C "$HOME"
bash "$HOME/rk265t-tiled-publication-kit-v1/launch.sh"
```

Status:

```bash
bash "$HOME/rk265t-tiled-publication-kit-v1/status.sh"
```

After completion or failure:

```bash
bash "$HOME/rk265t-tiled-publication-kit-v1/collect.sh"
```

Collection transfers `rk265t-review-JOBID.tar.gz` to murier Desktop, including
early-failure diagnostics. Status handles purged job IDs. Submission refuses
another job until the previous review is collected. Full results remain at the
path saved in `$HOME/rk265t-latest.json`.

Only if review finds the initial timing difference small/noisy/mixed, one
optional repeat is available after collecting the initial review:

```bash
bash "$HOME/rk265t-tiled-publication-kit-v1/repeat.sh"
```

This submits exactly six more timing runs. It verifies and reuses all four
binaries, source seals, toolchain and physical placement, without rebuilding or
repeating the seven validation runs. The extra runs retain the science gates.
Initial results are preserved in `results-before-repeat.json`; the original
review archive remains unchanged. Collect the repeat with the same collect.sh
command. A second repeat is rejected. No repeat is submitted automatically.

The independent reviewer accepts thirteen runs initially or nineteen after the
single repeat. It checks source/build identities, raw timing and state records,
placement, storage and exact-return proofs, and recomputes the comparisons.
Binary/checkpoint/VTK payloads are omitted: their cluster-computed hashes and
comparison reports are reviewed, not rerun locally. New sources and tile use
must pass local exact optimized/checked/oracle tests before kit sealing.
