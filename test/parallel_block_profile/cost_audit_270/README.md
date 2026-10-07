Stage270 prepares the next focused audit of unchanged Stage266 versus main.
It reduces observer frequency after Stage269's failed main overhead calibration,
then splits the large tendency, recomposition and migration regions into their
actual computation, transport, packing, installation and preparation costs.
No numerical or ownership optimization is included yet.

The exact controls remain Stage266 ff740744 and read-only legacy main 6a245db7
with the established checkpoint I/O hardening. Separate diagnostic copies carry
60 fixed record slots, including three retired slots that must remain empty.
There are no Domain boundary, per-patch velocity, node or physical-layer clocks.
New timings wrap whole routines, existing horizontal-level batches covering all
physical layers, and whole packing/unpacking/kernel batches. No synchronization
or timed-step output is added. Records flush once after the complete trajectory.

The major new splits are:

- Tendency: shared compatibility producer, native mass/velocity plan preparation,
  scalar capture preparation/publication/MPI, velocity remainder and RK assembly.
- Recomposition: direct/restricted flux and divergence kernels, ghost packing/
  local copies, posting/waiting and received installation, temperature edge
  transactions/waiting, final physical kernel, metric preparation and retention.
- Migration: manifest/checks, packing, size and payload exchange, unpacking,
  complete-store installation, ghost/writeback plans and workspace preparation.
  Payload MPI, full byte-checksum scans and proof reductions have separate scopes.

Counters retain inverse/cache work and add actual plan rebuild/reuse events,
logical migration bytes and remote scalar-ghost/temperature requests and doubles.
They do not measure physical memory bandwidth or all solver messages. Complete
records, conserved self times, retired slots, stage transactions and plan
accounting are independently verified. Stable float sums and the tested reviewer
rounding rule avoid cross-Python derived-report mismatches without loosening
raw hashes, integer counts or numerical acceptance gates.

One exclusive nk/bb13 allocation makes five serial builds, sharing the same 23
physics objects: Stage266, main, two optimized observers and checked block
observer. It needs only the existing J5 fixture at
`$HOME/stage187-fixture-20260911-145307`, plus your usual compiler/MPI modules.
No old results or kit directory is required on bbserv. Partition nk contains
only bb13, so no bb04 exclusion is needed.

There are fifteen solver runs. Seven correctness runs use 16 ranks for both
controls, both optimized observers, checked block observer and block oracle;
a 167-rank block oracle covers zero-Domain ranks. Checkpoint, exact same-family
trajectory/62 VTK members, restart/remap, storage, placement and independent
oracles remain. The oracle retains 48 exact inverse proofs and 84 installations.

Eight 83-rank calibration runs use the fixed order:
`main, main-observer, main-observer, main, Stage266, block-observer, block-observer, Stage266`.
Each binary gets two independent timing runs. Ordinary steps exclude steps 1,
8 and 9. Both calibrations must pass before cross-version phase attribution is
used. The unchanged limits are max(3 ms, 1%) mean shift, max(5 ms, 2%) per pair,
and 2% control/observer spread. Failures remain visible as
`observer_costs_unqualified`. No automatic repeat or promotion occurs.

Self times add within a local step; inclusive regions overlap. Signed differences
and negative offsets remain visible. The longest-local-step composition selects
one rank per step, not independent region maxima, and does not prove the global
critical path. CPU includes MPI spinning. Measurements locate hypotheses;
an unprofiled single-change comparison must establish recoverable gain. The
three-hour Slurm limit is a ceiling, not a runtime estimate.

On murier:

```bash
scp /Users/kevlahan/.codex/worktrees/26d8/wavetrisk_hydrostatic/.stage270-validation/rk270-reduced-cost-kit-v1.tar.gz kevlahan@bbserv.math.mcmaster.ca:~/
```

On bbserv:

```bash
tar -xzf "$HOME/rk270-reduced-cost-kit-v1.tar.gz" -C "$HOME"
bash "$HOME/rk270-reduced-cost-kit-v1/launch.sh"
```

Check status:

```bash
bash "$HOME/rk270-reduced-cost-kit-v1/status.sh"
```

After completion or failure, collect the review to murier Desktop:

```bash
bash "$HOME/rk270-reduced-cost-kit-v1/collect.sh"
```

The launcher saves the submitted ID and fresh results path in
`$HOME/rk270-latest.json`. It rejects another submission until the prior review
is collected. Status handles purged IDs; collection preserves early failures
and all raw observer files. The independent reviewer recomputes timings,
calibration, conserved records and science reports. Review archives omit
source/binary/checkpoint/VTK payloads; their cluster-computed hashes and reports
are verified without locally rerunning those payloads.
