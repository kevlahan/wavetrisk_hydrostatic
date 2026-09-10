# Restart mismatch investigation — 179d preparation

## Scope and status

2026-09-10 initial investigation, not a performance stage. This document records
the **pre-correction baseline**. At that point no production arithmetic changes
or commits had been made. Implementation and subsequent validation are recorded
in `parallel-blocks-restart-corrections-179d.md`; do not interpret the baseline
findings below as the current validation status.

Authoritative `main` remains `b82647b36eb825aed0c0aae4e92adcf4ab46d608`.
Block HEAD remains `2eee9646312b3b4ed25c2b6f137d519bb8906ca9`.
Separate source/build copies were instrumented; neither `main` nor the original
checkpoint-3 fixture was modified. No tolerance, grid-level cap, vertical layer
count, or initialization change was used to hide the discrepancy.

## Completed experiment

Both isolated optimized RK4 builds completed the standard checkpoint-3-to-4
restart experiment sequentially with four MPI ranks, J7, 30 atmospheric layers,
and the original adaptation tolerance. End time was set to 0.3105 day to include
checkpoint 4 and its immediate reload. These heavily instrumented local runs
are not timing measurements.

The new `test/parallel_block_profile/restart_probe.f90` writes read-only solution,
wavelet, and mask snapshots before/after checkpoint operations. It does not
refresh boundaries, synchronize MPI, or materialize block values itself.
The comparator uses global Domain IDs and canonical patch-tree keys, independent
of local patch numbering. It includes all vertical layers, but does **not** dump
boundary aliases or special pole storage. Thus conclusions below are scoped to
captured patch interiors. The diagnostic scaled threshold is not an oracle
acceptance tolerance.

| Comparison | Matched patches | Finding |
| --- | ---: | --- |
| Initial checkpoint-3 restart-ready | 2,894 | All captured values and masks numerically exact |
| Before checkpoint-4 serialization | 3,144 | Same topology/masks; fields already differ |
| Block before/after final consumer synchronization | 3,144 | No captured value changes |
| After serialization | 3,144 | Coarse scalar solution differences increase |
| After checkpoint-4 reload | 2,387 | Same topology/masks; scalar reconstruction differences persist |

No missing patch keys or nonfinite values occurred in these comparisons.
At checkpoint 4, both runs also had the same time, cumulative step 259, and
remap counter 4. The final writeback is not the first source of the captured
discrepancy. Initial restart reconstruction is not the first observed source
either; the fields become different during the intervening evolution.

## Concrete scalar coverage defect

Immediately before serialization, the block run has zero coefficients where
legacy has **230,400 nonzero level-5 wavelet coefficients for mass, and another
230,400 for temperature**. Maximum absolute differences are respectively
2.0407329691852283 and 717.3876394510444. The count equals
640 root patches × 12 wavelet sites × 30 atmospheric layers.

In contrast, level-6/7 scalar solution differences before serialization are at
floating-point scale (maximum mass about 3.8e-12 and temperature about 1.5e-9).
Large coefficient differences are concentrated at the block-root level.

`produce_block_scalar_wavelets` in `src/parallel_block_mpi.f90` iterates retained
parents and their children. This does not produce the root patch's coefficients
when its parent is outside the extracted block subtree. The final transform
requests the `level_start-1` parent range, but lowering that range does not add
the absent parent. The structural-zero import therefore leaves those root
coefficients unproduced. Legacy's corresponding traversal includes that parent.

Checkpoint `Restrict_scalar` uses wavelet coefficients in its lifting operation.
The observed increase in coarse solution differences during serialization is
consistent with consuming the missing coefficients. Reload then reconstructs
different fine values. This is an identified coverage defect, **not yet a tested
production correction or a complete explanation of every discrepancy**.

There are also velocity differences before serialization (maximum captured
active solution difference about 0.01655). These need independent localization;
do not declare the whole mismatch fixed after correcting scalar root coverage.

## Next correction and acceptance gates

1. Produce scalar coefficients over eligible **target patches**, including roots
   whose parent is outside the block. Preserve provisional/final transform
   ranges, masks, stencil arithmetic, and existing boundary provenance. Do not
   introduce Domain production shadow computation as the fix.
2. Add coverage checks derived independently from the requested transform range:
   expected target sites per level, explicitly including block roots. Existing
   visited-parent counts cannot detect parents omitted from the traversal.
3. Add a root-patch regression with a nontrivial scalar field and external parent;
   compare native coefficients against the legacy stencil before compression.
4. Run one-step/phase-localized legacy/block comparisons before repeating the
   checkpoint-4 experiment. Check numerical state at authoritative/materialized
   boundaries; stale intermediate Domain values alone are not failures.
5. Localize any remaining velocity/soil discrepancy, then require RK3 and RK4
   oracle checks and production restart equivalence. Only then resume performance
   work and the matched 83-rank benchmark. Do not relax tolerances to accept
   missing coefficients or change the original fixture.

## Reproducibility and files

Helpers:

- `prepare_restart_probe.py`: creates and instruments a new external source/build copy.
- `run_restart_probe.py`: sequential isolated four-rank restart runs, all validation/profile
  switches off and `WAVETRISK_RESTART_PROBE=1`; records only relevant environment settings.
- `compare_restart_probe.py`: compares binary snapshots by Domain/tree key.
- `test_restart_probe.py`: synthetic identity, field-layout/sensitivity, nonfinite,
  missing-key, and missing-file tests. Full profiling helper suite: 18 tests pass.

Completed builds:

- `/private/tmp/wavetrisk-restart-probe-legacy-02`
- `/private/tmp/wavetrisk-restart-probe-block-02`

Completed run data and full binary snapshots:
`/private/tmp/wavetrisk-restart-probe-runs-01/{legacy,block}`.
These temporary directories are diagnostic scratch, not permanent evidence storage.
Compact comparison JSON is retained in `docs/restart-mismatch-179d/`.

Example comparison (run from repository):

```sh
python3 test/parallel_block_profile/compare_restart_probe.py \
  --left '/private/tmp/wavetrisk-restart-probe-runs-01/legacy/probe-before-dump-step-259-rank-*.bin' \
  --right '/private/tmp/wavetrisk-restart-probe-runs-01/block/probe-before-dump-step-259-rank-*.bin'
```

Do not use whole compressed checkpoint hashes as a numerical equivalence test:
metadata/load estimates and representation need separate treatment. Compare
keyed field values, topology, and relevant run state instead.
