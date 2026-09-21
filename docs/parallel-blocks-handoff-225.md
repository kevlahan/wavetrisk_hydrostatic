# Stage225 handoff: empty-Domain oracle correction validated locally

## Cluster completion update: job 54729

The prepared campaign has now run and passed all eleven gates on bb13.
All 3,052 archive hashes were verified; raw lifecycle, placement and adaptive
trajectory evidence was rechecked against the report. All 17 numerical
comparisons are exact. Both 41/40 and 167/160 pass with checked, optimized and
oracle builds; the empty-root/nonzero checkpoint-source case also passes.
Ranks 160–166 own zero Domains but 1–6 blocks and complete all twelve steps.
See `parallel-blocks-stage225-review.md` and `.json` for the completed review.
The preparation-era instructions below remain as provenance; another repeat
of the unchanged matrix is not needed. Next is a bounded production integration
of the validated capability, retaining restrictions on unvalidated paths and
checking the integrated boundary cases. No production guard change, commit or
merge was made during this review; main remains unchanged.

Implemented the allocation-consistent transform-stage extent check in
`src/parallel_block_mpi.f90`. The sole production change relative to Stage224
is one expression (`max(1, terminal_displacement)`) plus explanatory comments.
The patch is `test/parallel_block_profile/stage225-oracle-extent.patch`.
The checkpoint prefix fix from Stage224 remains in place.

Both local four-rank J4 layouts pass full dynamics/adaptation/remap oracles:
source rank 0 owns one Domain and rank 3 owns none; source rank 1 owns one
Domain and rank 0 owns none. Each completes 12 steps, two checkpoint reloads,
two remaps, exact restored metadata, an exact final checkpoint comparison and
a matching adaptive trajectory. All ranks complete useful block work. No
additional solver correction was needed. See
`parallel-blocks-stage225-oracle-fix.md` and
`parallel-blocks-stage225-local-results.json` for reproducible commands,
source/build/log hashes and the ten-case predicate regression.

Final local runs use DEBUG=check with -O1 and all runtime, initialization,
FP-trap and contraction-off settings retained. Five optimization-dependent
warnings were reviewed and retained as nonfatal; details are in the report.
The initial O0 run was deliberately stopped after successful oracle progress
because of cost, and is not counted as a complete validation. Production
Makefile and cluster compiler policy are unchanged.

The Stage225 cluster kit is now prepared (not submitted):
`/Users/kevlahan/Desktop/wavetrisk-stage225-cluster-kit.tar.gz`.
The source is the frozen Stage224 reference plus only the oracle extent fix;
all 136 packaged Fortran files match the completed local regressions after
reproducing their test instrumentation. Archive re-extraction and candidate
patch reversal/reapplication reproduce both source inventories. The Stage224
kit remains frozen. Source and archive identities are recorded in
`parallel-blocks-stage225-kit.json`.

All 17 packaged launcher/collection/capability/checkpoint tests passed. The
isolated complete workflow passed real extraction, state backup and duplicate
guards, same-node launch selection, failure archive creation/verification and
repeat transfer reuse with simulated Slurm/scp. All packaged Python helpers
parse and all Stage225 shell/batch scripts pass bash syntax checks.

Next: transfer and run the fresh Stage225 build and eleven bounded correctness
gates following `parallel-blocks-cluster-225.md`. This repeats quick J4 anchors,
then oracle/empty-root and J5 160/167 gates. Six cluster builds retain established
checked O0 / optimized policies without local O1 warning overrides. Do not
skip or disable oracle gates or infer performance gains from these tests.
The repaired nk placement uses 2x128 physical cores, the same build/run node,
and no explicit step CPUs-per-task/exact restriction. `collect_225.sh` combines
archive creation, verification and scp FROM bbserv TO murier, and safely reuses
the archive after an interrupted transfer. It does not query purged job IDs.

The production startup task-count guard remains unchanged. Git main remains
b82647b36eb825aed0c0aae4e92adcf4ab46d608. Existing uncommitted Stage222/224 work
is preserved; no commit or merge was made in this turn.
