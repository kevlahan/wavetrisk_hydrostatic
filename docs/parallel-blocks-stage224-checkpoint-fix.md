# Stage224: immutable checkpoint metadata position

Fix the concrete restart failure from Stage223 job 54725. This is a correctness
repair and local regression, not a performance campaign or promotion of the
above-Domain-count prototype.

## Correction

`load_adapt_mpi` saves `global_prefix_position` before reconstruction. It is
zero on ranks with no local Domain and otherwise holds the original position
of the first Domain's metadata. `distribute_checkpoint_global_state` accepts
that immutable scalar instead of the reconstruction cursor array. The source
reads from it; receivers with position zero restore the broadcast metadata.

The fix adds one scalar and preserves the checkpoint format and existing
metadata-only broadcast. It does not add Domain payload replication. The
production startup task-count guard and Git `main` remain unchanged.

`test/parallel_block_profile/stage224-checkpoint-prefix.patch` contains only
the correction relative to the frozen Stage223 reference checkpoint source.
The current candidate generator's restart-event anchor is updated for the new
call. Previously packaged kits, source manifests and binaries remain frozen;
a cluster rerun will require a new source identity and rebuilt candidate.

## Local regression

`test/parallel_block_profile/test_checkpoint_prefix_224.py` creates a disposable
copy of the preserved Stage223 prototype. Only this test tree changes Domain
ownership and exports restart metadata. The real block distributor still runs.
Four MPI ranks own all 40 J4 Domains in two deliberately uneven layouts:

| Layout | Domains per rank (0, 1, 2, 3) | Metadata source | Empty receiver |
| --- | --- | --- | --- |
| single-source0 | 1, 20, 19, 0 | 0 | 3 |
| single-source1 | 0, 1, 20, 19 | 1 | 0 |

Both use checked builds with runtime bounds checks, initialization checks,
floating-point traps and contraction disabled. The original checkpoint source
is a negative control: both layouts must fail at its erroneous prefix read.

The corrected cases must complete the full adaptive J4 fixture, reload input
checkpoint 2 and newly written checkpoint 3, and complete both vertical remaps.
On every reload the test independently reads the original record prefix from
the source scratch stream, serializes each rank's restored global state, and
requires byte-for-byte equality on all ranks. These bytes cover cumulative
step, time, integer clocks/output index, thresholds and any case dump payload
(the climate case callback is empty).

Per-rank lifecycle evidence must show exactly one Domain on the metadata
source, zero on the designated receiver, useful block work on every rank,
agreement in global clocks, and a complete dynamics/step sequence. The final
checkpoint and adaptive trajectory are compared against the prior validated
local fixture; the established cross-layout numerical policy is unchanged.

The test driver can reuse its own checked build after comparing every Fortran
source with a freshly generated test tree. This was needed because the sandbox
initially blocked MPI's local sockets; the subsequent runs use the verified
same binaries outside that restriction. No solver result is taken from the
failed socket launches.

Both complete trajectories passed. Results are recorded in
`parallel-blocks-stage224-local-results.json`:

- Original code: both layouts reproduce EOF at the erroneous prefix read.
- Corrected code: both layouts complete 12 steps, reloads 2 and 3, and two remaps.
- All four ranks match the independent 1,008-byte metadata prefix exactly on
  both reloads in both layouts, including each zero-Domain receiver.
- Both final checkpoints are exactly equal to the previously validated local
  reference; adaptive trajectories agree.
- Every rank completes block dynamics; metadata sources remain single-Domain
  owners and the designated empty ranks remain without Domains.
- All 13 Stage223 helper tests and four checkpoint-reader tests pass. The
  prototype generator accepts both frozen and corrected source forms.

Reproduce the corrected full regression (fresh output directory required):

```bash
python3 test/parallel_block_profile/test_checkpoint_prefix_224.py \
  --prototype /private/tmp/stage223-local/final-tree \
  --checkpoint "$PWD/src/checkpoint.f90" \
  --fixture /private/tmp/stage223-local/empty-root-checked \
  --reference /private/tmp/stage223-local/empty-root-checked \
  --out /private/tmp/stage224-repeat
```

For the negative control, use the preserved prototype's
`src/checkpoint.f90` as `--checkpoint`, a different fresh `--out`, and
`--negative-control`. MPI needs permission to open local communication sockets.
 The local regression targets the
source/receiver cardinality that escaped earlier tests. It is not evidence
that the full 41/40 and 167/160 cluster capability matrix has passed.
