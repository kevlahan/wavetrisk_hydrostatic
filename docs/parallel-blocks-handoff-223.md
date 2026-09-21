# Stage223 handoff

Stage223 is packaged, verified and ready for bbserv; it has not been submitted.
Use `parallel-blocks-cluster-223.md` for the complete transfer/build/run/retrieval
commands. The archive is `~/Desktop/wavetrisk-stage223-cluster-kit.tar.gz` on
murier. Its checksum and preparation evidence are in
`parallel-blocks-stage223-kit.json`.

The user confirms nk has two 128-core sockets (256 physical cores). Reserve
one exclusive node with all 256 cores; discover and verify physical CPU IDs.
Mappings alternate sockets: 20+20, 21+20, 80+80 and 84+83 ranks. Build and run
on the same node. Do not assume CPU IDs are contiguous or confuse SMT siblings
with physical cores.

The local checked empty-root case passes the complete 12-step trajectory with
an exact Stage215 checkpoint match. Rank zero owns blocks without Domains;
rank one supplies checkpoint metadata. The kit freezes the validated Stage222
sources as baseline and applies only a bounded prototype patch in a separate
tree. The production startup guard and Git `main` branch remain unchanged.

The eleven-run sequence tests J4 40/41 first, including checked, optimized,
oracle and empty-root cases, then J5 160/167. Stop at the first failure, retaining
the results and per-rank lifecycle files. This is a capability gate, not a
performance campaign. Do not apply a speedup threshold or promote the prototype
automatically. Review zero-Domain ranks' completed block work, global clocks,
restart metadata, adaptive trajectory, topology, thresholds and numerical
comparisons before replacing the production restriction.

The packaged workflow was tested with fake Slurm and scp commands: real source
extraction, preservation of old state, duplicate submission guards, same-node
selection, failure collection and safe archive reuse all pass. Retrieval must
run FROM bbserv TO murier. `collect_223.sh` combines archive creation/reuse,
verification and scp, and queries the user's queue rather than a purged job ID.

### Launcher repair, 2026-09-19

User reports job 54724 on bb13 failed at the first 40-rank placement probe,
before any solver runs. See `parallel-blocks-stage223-launcher-repair.json` and
the repair section of `parallel-blocks-cluster-223.md`. A small repair archive
reuses all six builds; no numerical/source change. Original helper stderr was
hidden; the explicit `--cpus-per-task=1` step restriction has been removed and
probe diagnostics are preserved. Await the repaired cluster run before drawing
any conclusion about ranks without Domains.

### Reviewed job 54725, 2026-09-21

This supersedes the pending-run status above. See
`parallel-blocks-stage223-review.md` and `.json`. All placement probes pass;
40-rank baseline/candidate complete with exact checkpoint equality. At 41 ranks
rank zero fails after 6.1 seconds with EOF reading the global checkpoint prefix;
143 is subsequent peer termination. Root cause is reusing `domain_pos(1)` after
`load_domains_stream` advances it. The source-extracted Fortran probe confirms
single-record EOF and how multi-record ownership masked it locally. No source
fix or rebuilt kit has yet been applied. Next: immutable prefix location and
local one-Domain-source/zero-Domain-receiver validation before cluster rerun.


### Correction and local regression, 2026-09-21

The checkpoint cursor correction is now implemented and locally validated;
this supersedes the previous "no source fix" status. See
`parallel-blocks-stage224-checkpoint-fix.md` and
`parallel-blocks-stage224-local-results.json`. Both one-Domain-source layouts
(source rank 0 / empty rank 3; source rank 1 / empty rank 0) reproduce EOF with
the old source and pass full checked MPI trajectories with the correction.
All restart metadata bytes and both final checkpoint comparisons are exact.
The frozen cluster kit has not been changed or rebuilt; the complete 41/40
and 167/160 capability matrix remains pending. Do not resubmit the old binary
or lift the production task-count guard on the local regression alone.

### Next cluster campaign prepared

Stage224 is now packaged with the corrected checkpoint metadata position and
the repaired whole-node launcher. See `parallel-blocks-handoff-224.md` and
`parallel-blocks-cluster-224.md`. It uses fresh builds/state/results, preserving
all Stage223 evidence. Cluster submission is pending; do not reuse Stage223
binaries for the corrected test.
