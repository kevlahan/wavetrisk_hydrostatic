# Stage225: empty-Domain oracle extent correction

Stage224 job 54727 completed its non-oracle 40/41-rank gates exactly, but its
oracle-enabled run aborted during `check_block_writeback_plan` on rank 40.
The allocator creates `max(1, logical_count)` vector-transform storage, including
one dummy element on a zero-Domain rank. The checker incorrectly required
allocated capacity to equal the unpadded logical count.

The correction compares capacity with `max(1, terminal_displacement)`, matching
the allocator and adjacent scalar/vector buffer checks. Initial displacement
and matching transform/wavelet/coverage capacities are still checked. No
numerical comparison, collective participation or validation switch is removed.
The production change is one expression plus explanatory comments in
`src/parallel_block_mpi.f90`; the patch against the frozen Stage224 reference
is `test/parallel_block_profile/stage225-oracle-extent.patch`.

## Validation scope

The exact current predicate is tested in checked Fortran for valid empty and
nonempty buffers, a missing empty-rank dummy, a wrong nonempty capacity and a
wrong coverage-array size. The corrected predicate accepts valid cases and
rejects the malformed cases; the reconstructed original predicate reproduces
the false rejection of the valid empty rank. See
`parallel-blocks-stage225-oracle-extent-probe.json`.

The existing full-solver regression driver now supports `--oracles` and an
explicit `--mpi-source` overlay. Test-only ownership and metadata exports stay
in disposable trees. Both four-rank J4 layouts use the real block distributor:

| Layout | Domains on ranks 0,1,2,3 | Source / empty receiver |
| --- | --- | --- |
| single-source0 | 1,20,19,0 | 0 / 3 |
| single-source1 | 0,1,20,19 | 1 / 0 |

The final local build uses `DEBUG=check` with an explicit `-O1` override,
retaining all runtime/bounds checks, initialization checks, floating-point traps
and disabled contraction. All 55 solver compilation commands were checked for
these flags. An initial O0 run passed its first full dynamics step but was
intentionally stopped because each step took about 140 seconds; it is not
counted as a completed regression.

GCC 15's optimization-dependent `maybe-uninitialized` and `function-elimination`
warnings remain visible but are nonfatal for this local O1 build. The five
warnings were reviewed: one concerns an allocatable descriptor allocated in
the first pass of a fixed two-pass loop, and four concern short-circuiting of
already initialized validation-option queries. Other warnings remain errors;
no runtime check or oracle is disabled. No production Makefile is changed. All three validation switches (dynamics, adaptation,
remap) are explicitly enabled. Required completion is 12 steps, initial and
new-checkpoint reloads, two remaps, useful block work on every rank, exact
restored metadata on every rank and comparison with the validated checkpoint
and adaptive trajectory. Full-oracle runs are correctness tests, not timing
measurements.

Reproduction (use a fresh output directory):

```bash
python3 test/parallel_block_profile/test_checkpoint_prefix_224.py \
  --prototype /private/tmp/stage223-local/final-tree \
  --checkpoint "$PWD/src/checkpoint.f90" \
  --mpi-source "$PWD/src/parallel_block_mpi.f90" \
  --fixture /private/tmp/stage223-local/empty-root-checked \
  --reference /private/tmp/stage223-local/empty-root-checked \
  --out /private/tmp/stage225-oracles-repeat --oracles --check-optimization O1
```

MPI needs permission to open local sockets. Both full runs passed; results and source/build identities are recorded in
`parallel-blocks-stage225-local-results.json`.

- Both layouts completed 12 steps, checkpoint reloads 2 and 3, and two remaps.
- Dynamics, adaptation and remap oracles were enabled throughout; the logs
  confirm both remap oracle completions in each run.
- All four ranks matched the independent 1,008-byte checkpoint metadata prefix
  exactly on both reloads in both layouts.
- Both newly written checkpoints are exactly equal to the validated reference,
  and the complete adaptive trajectories match.
- Every rank completed useful block work, including the designated zero-Domain
  receiver, while each metadata source retained exactly one Domain.
- No further solver/oracle correction was needed after the extent fix.
 The old Stage224
cluster kit/binaries remain frozen. The production task-count limit and Git
`main` remain unchanged; the cluster oracle, empty-root and J5 160/167 gates
still require a rebuilt candidate after local validation.
