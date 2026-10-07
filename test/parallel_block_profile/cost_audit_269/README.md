Stage269 measures the largest remaining costs before selecting one optimization.
The production solver remains exactly accepted Stage266 ff740744; the read-only
legacy main anchor remains 6a245db7 with its established checkpoint I/O hardening.
Separate diagnostic source copies add bounded transaction observers only.

One exclusive nk/bb13 allocation makes five serial builds, sharing the same 23
physics objects: Stage266, main, their two observers, and checked Stage266 with
the observer. No prior results or old campaign kit is required on bbserv. The
only external prerequisite is the established J5 fixture at
`$HOME/stage187-fixture-20260911-145307`, plus the usual compiler/MPI modules.
Partition nk contains only bb13, so no bb04 exclusion is needed.

Seven correctness runs use 16 ranks for both controls, both optimized observers,
the checked observer and the block oracle; a 167-rank block oracle exercises
zero-Domain ranks. Checkpoint, exact same-family states/VTK, restart, remap,
storage, placement and independent oracle gates remain. All block configurations
match Stage266 exactly. The oracle retains all 48 inverse proofs and 84 actual
installations. No numerical or ownership change is proposed in this kit.

Eight 83-rank runs calibrate observer overhead in this fixed order:
`main, main-observer, main-observer, main, Stage266, block-observer, block-observer, Stage266`.
There are 15 solver runs total, two independent timing runs per binary. Use the
precise ordinary-step window excluding steps 1, 8 and 9. No automatic repeat or
promotion occurs. The Slurm request's three-hour limit is a ceiling, not an
estimated runtime.

The observers cover 27 nested transaction families: RK/tendency, wavelets,
Domain/native inverse, inverse setup/wait/publication/installation, provisional
and accepted refresh, adaptation, catalogue/source extraction/migration,
post-grid reconstruction, hydrostatic cache, native recomposition/velocity,
compression/forward transforms/ghosts, surface pressure and Domain boundaries.
They record per-rank/per-step self and inclusive wall/process CPU times, calls,
and named work counters. Actual work includes rebuilt cache columns/layers,
plan rebuilds, published/installed slots and inverse remote requests/payload.
There are no node/layer clocks, added synchronization or output during a timed
step. Fixed buffers are flushed once after the solver run.

Both observers must pass their own calibration before current region differences
are used to prioritize work. Prespecified limits are 3 ms or 1% for mean shift,
5 ms or 2% per pair, and 2% spread, using the larger absolute/relative limit.
Failures remain visible as `observer_costs_unqualified`; all raw evidence is
collected. An overhead pass does not establish that each tiny scope is unbiased.
Signed mean self times are additive; inclusive times overlap. The composition
of the longest local step uses one rank per step, not independent region maxima.
Local unsynchronized durations do not establish the global critical path, and
process CPU includes MPI spinning. Counters describe logical payloads and named
products; they do not measure all solver traffic or physical memory bytes.

This campaign locates costs; it does not establish recoverable speedup. Select
one demonstrably redundant operation after review, then qualify that intervention
with the short unprofiled main/Stage266/candidate comparison. Preserve the
Stage267 shared-reader and boundary requirements. Local four-rank validation
establishes correctness and complete conserved records, not bb13 performance.

On murier:

```bash
scp /Users/kevlahan/.codex/worktrees/26d8/wavetrisk_hydrostatic/.stage269-validation/rk269-cost-audit-kit-v1.tar.gz kevlahan@bbserv.math.mcmaster.ca:~/
```

On bbserv:

```bash
tar -xzf "$HOME/rk269-cost-audit-kit-v1.tar.gz" -C "$HOME"
bash "$HOME/rk269-cost-audit-kit-v1/launch.sh"
```

Check status:

```bash
bash "$HOME/rk269-cost-audit-kit-v1/status.sh"
```

After completion or failure, collect to murier Desktop:

```bash
bash "$HOME/rk269-cost-audit-kit-v1/collect.sh"
```

The launcher saves the submitted ID and fresh results path in
`$HOME/rk269-latest.json`. It rejects another submission until the prior review
is collected. Status handles purged IDs; collection preserves early failures
and raw observer files. The independent local reviewer recomputes timings,
calibration, observer conservation and numerical reports after transfer. The
review omits source/binary/checkpoint/VTK payloads; their cluster-computed hashes
and reports are checked, without rerunning those payloads locally.
