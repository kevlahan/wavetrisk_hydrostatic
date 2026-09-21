# Stage223 review: job 54725

The archive passes verification of all 277 inventoried files. The launcher
repair succeeded: 40/41/160/167-rank placement probes all passed, each with a
unique pinned CPU per rank. This is placement evidence, not solver completion
at 160 or 167 ranks.

Both J4 40-rank baseline and candidate complete 12 steps, checkpoint reloads
2 and 3, and two remaps. Their final checkpoint hashes are identical; all
numerical gates pass exactly. These short correctness runs support no speedup
claim.

The J4 41-rank checked candidate fails after 6.1035 seconds during its initial
checkpoint reload. The primary error is **Fortran end-of-file** on rank zero,
`src/checkpoint.f90:1710`, in `distribute_checkpoint_global_state`. Rank zero
exits 2; Slurm then terminates the other 40 ranks under kill-on-bad-exit. The
reported 143 is consequential, not evidence of a 15-minute timeout or deadlock.
No 41-rank dynamics, later empty-root cases, or J5 solver gates completed.

## Source cause and independent reproduction

`load_adapt_mpi` records initial Domain positions in the local scratch stream,
then calls `load_domains_stream`. That loader deliberately updates each
`domain_pos(d)` as it consumes the record. The subsequent global metadata
broadcast reads its prefix using the now-advanced `domain_pos(1)`. It should
read from the original first-record start.

At 41 ranks for 40 Domains, the source rank owns one Domain. The advanced
position is past its only record, so the prefix read reaches EOF. The earlier
four-rank empty-root test assigned many Domains to each nonempty rank. Reading
past the first record could then land on the next concatenated record's prefix,
which contains the same global metadata, masking the cursor error.

`probe_checkpoint_prefix_223.py` extracts the actual failing subroutine from
the preserved candidate tree. Its source SHA256 matches the cluster build's
checkpoint source. It exercises actual checked Fortran stream I/O with stubbed
MPI broadcasts:

| Prefix seek | One Domain record | Two Domain records |
| --- | --- | --- |
| Original advanced cursor | EOF, exit 2 | Exact prefix passes accidentally |
| Original record start | Exact prefix passes | Exact prefix passes |

This isolates the I/O defect and explains the previous false assurance. It is
not a full solver or MPI validation. The experiment changes only generated
test source; production source and `main` remain untouched.

## Next bounded step

Preserve the first record's original position independently of the mutable
reconstruction cursors (or retain an immutable copy of the prefix). Add a local
MPI correctness gate combining a source with exactly one Domain and a receiver
with none. Check global clocks/model state on the receiver and both initial
and newly written checkpoint reloads. Retain the source-not-rank-zero case too.
Then rebuild the changed candidate and repeat the Stage223 correctness gates,
starting at J4 40/41 before J5 160/167. Do not increase the timeout or lift the
production task-count guard on these results.
