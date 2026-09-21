# Stage224 handoff

Fresh cluster campaign prepared; not submitted. Desktop artifact:
`/Users/kevlahan/Desktop/wavetrisk-stage224-cluster-kit.tar.gz`.
See `parallel-blocks-stage224-kit.json` for checksum/build provenance and
`parallel-blocks-cluster-224.md` for complete transfer/build/run/retrieval.

The reference is the frozen Stage223 reference plus only the corrected
checkpoint metadata position. Candidate retains the opt-in Domain/task limit
bypass and lifecycle evidence. The locally proven source is packaged; test-only
single-source ownership maps and metadata exports are not in the cluster solver.
Reproducing that instrumentation matches all 136 Fortran files in the local
validation source inventory. The cluster candidate uses real native Domain and
block ownership, plus the existing empty-root stress case.

Eleven gates: J4 40 reference/candidate, 41 checked/optimized/oracle/empty-root;
then J5 160 reference/candidate and 167 checked/optimized/oracle. Stop on failure.
The exact local fix validation is in `parallel-blocks-stage224-local-results.json`.
The campaign is correctness/capability only; checkpoint saving is required.
Do not infer timing gains or remove the production task-count guard automatically.

nk: two 128-core sockets, one exclusive node. Preserve fixed whole-node Slurm
step policy (no explicit step CPUs-per-task or exact), manual physical-core
pinning, visible probe diagnostics, same build/run host and six fresh binaries.
The Stage223 runtime opt-in selector/event names are intentionally retained;
all orchestration state, results and source/build identities use Stage224.

Package re-extraction, candidate patch reverse/forward application, helper
hashes, imports, eight launcher/collector tests and the full fake-Slurm/scp
workflow pass. `main` remains b82647b36eb825aed0c0aae4e92adcf4ab46d608.
Use `collect_224.sh` after either success or failure: one command creates,
verifies and transfers the archive FROM bbserv TO murier without sacct or
queries for purged job IDs. Keep previous Stage223 artifacts.

## Job 54727 reviewed

This supersedes the prepared/not-submitted status above. See
`parallel-blocks-stage224-review.md` and `.json`. All 40/41 non-oracle cases
pass with exact final checkpoints, including the zero-Domain rank's completed
block dynamics and both restarts. The oracle run explicitly aborts in the
writeback-plan extent validator: logical zero is compared with a one-element
dummy allocation. A checked source-predicate reproducer confirms the issue;
no production fix has yet been applied. Next: correct the validator and run
both local one-source/empty-receiver layouts with full oracles enabled, then
continue the cluster correctness matrix. No evidence yet for empty-root or
J5 160/167 solver completion.

### Oracle correction locally validated

Stage225 now corrects the empty-Domain extent check. Both one-Domain-source /
zero-Domain-receiver layouts pass complete checked local runs with all three
oracles enabled, including source rank 1 and empty rank zero. Metadata and
new checkpoint comparisons are exact. See `parallel-blocks-handoff-225.md`.
The Stage224 cluster kit remains frozen; a fresh candidate rebuild is needed
before continuing the outstanding cluster gates.
