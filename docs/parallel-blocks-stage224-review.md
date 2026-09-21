# Stage224 review: job 54727

Verified all 567 inventoried files in `stage224-review-54727.tar.gz`. The user
initially attached the outgoing kit; the correct review archive was already on
the Desktop and was subsequently attached explicitly.

## Established results

All four completed solver cases passed: J4 40-rank reference and optimized
candidate, followed by 41-rank checked and optimized candidates. Each completed
12 steps, checkpoint reloads 2 and 3, and two remaps. All four final checkpoints
have the identical SHA256
`270a3d6aec3ce3bc58c19fcad8f1c932adeeeb69d0e5674fa89dd234c93f25ab`.
All six numerical comparison gates passed exactly.

In both 41-rank runs, rank 40 owned zero Domains throughout the recorded
lifecycle while owning useful blocks and completing every dynamics/step event.
This is successful native-layout evidence for 41 ranks over 40 Domains,
including checked execution and checkpoint continuation. It is not yet the
complete intended capability matrix or a performance result.

## Failure and source cause

The oracle-enabled 41-rank case ran 5.93 seconds and completed initial restart
and installation before rank 40 explicitly aborted with:

`parallel_block_mpi_mod: complete vector transform stage extent mismatch`

The log records `MPI_Abort(..., 1)` on rank 40, followed by Slurm killing all
ranks. Exit 137 is a consequence; the primary failure is not an OOM termination,
a timeout, or a reported numerical mismatch. Failure occurs inside
`check_block_writeback_plan`, before the later numerical oracle exchanges.

`build_block_writeback_plan` correctly allocates transform, wavelet and
coverage arrays with `max(1,total_vector_transform)`. On a zero-Domain rank,
`total_vector_transform` and the terminal displacement are zero, while each
buffer has one dummy element. The oracle compares the raw terminal
displacement to allocated buffer size and therefore compares zero with one.
Neighbouring payload/Domain-stage checks already account for this minimum-one
allocation convention.

All four candidate source files were reproduced from the frozen kit and match
the cluster build hashes. `probe_oracle_extent_224.py` extracts the exact
failing predicate from that source and runs checked Fortran cases. It confirms
that the original check rejects valid empty-rank storage and accepts an empty
zero-length allocation inconsistent with the allocator. A generated candidate
predicate comparing `max(1,terminal displacement)` with capacity accepts valid
empty and nonempty cases and rejects the tested malformed capacities/coverage
lengths. This is a predicate-level experiment, not full MPI oracle validation;
no production solver source has been changed during this review.

## Next step

Make the plan's empty-rank extent check consistent with its allocation policy,
retaining consistency checks rather than disabling oracles. Run local full
oracles in both one-Domain-source/zero-Domain-receiver layouts, including source
rank 1 with empty rank zero, and retain complete adaptive/restart/numerical
gates. Earlier Stage224 local runs deliberately had oracles disabled, so they
did not exercise this checker.

After local validation, resume the cluster capability matrix. The 41-rank
oracle run, 41-rank empty-root case and every J5 160/167 solver gate remain
incomplete or unrun. Keep the production startup restriction until the planned
coverage has passed. Do not increase memory/time limits to address this failure.
