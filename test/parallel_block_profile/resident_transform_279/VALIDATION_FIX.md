# Stage279 v3 — transform census validation correction

Job55506 ran for10:28 on bb13. All five builds succeeded; five16-rank validation
runs were marked passed. The83-rank solver also completed, but the evidence
checker raised an exception before marking that run passed. No bare performance
measurements or160/256 runs began. Preserve the original failed result unchanged.

The checker incorrectly required each canonical owner to have nonzero scalar,
inner-vector AND outer-vector descriptor counts. Placement assigns each scalar
or inner operation by its output and outer dependency components as a unit.
Canonical ownership therefore does not imply possession of every kernel class.
Ranks10,41,42,43 have legitimate outer-only plans during24 rank-transactions.
For example, rank42 transaction25 owns352 canonical nodes with0 scalar,0 inner
and573 outer records. Calling this no transform work was incorrect.

V3 requires positive total descriptors for a canonical owner, still checks every
per-block census, and requires all three classes globally at every transaction.
The runner and independent reviewer additionally require identical global class
counts across all audited validation runs of the fixed trajectory. All36 vectors
already agree exactly across production16, checked16, oracle16 and production83.
Descriptors are compiled records, not measured FLOPs or proof of active arithmetic
by themselves; the production branch, science/oracles and retirement checks remain.
Missing records, incorrect traffic, contexts, epochs, reference flags and absent
zero-Domain participation remain failures. No check on scientific output is relaxed.

Read-only reanalysis passes all six completed cluster runs and the three local
v2 audits. Checkpoint/binary/VTK hashes were independently rechecked in place on
bbserv; all62 VTK members per run and exact block/main-relative science reports
pass. The actual83-rank records are included in the regression tests.40 helper
tests pass, including outer-only ownership, absent work, missing global classes,
cross-count drift, and shared current revision identity. A stale v1 identity in
the final reviewer was also corrected; runner/preflight/reviewer share one value.

Solver source archive, source manifests, compiler flags and scientific fixture
are byte-identical to v2. Its GCC13 checked/optimized builds and three fresh local
full replays remain valid. No solver rerun or remote mutation was required for
this correction. V3 uses the unchanged49-run schedule; no performance result,
scaling improvement, memory improvement or promotion follows from this fix.
The v1/v2 kits and both failed-job records are retained. Main is untouched.

V3 is prepared for user transfer/submission. Collect job55506 with the v2 kit
before launching v3. The new run uses a fresh directory and independently built
binaries; the failed validation job is not silently relabeled as a completed
performance campaign.
