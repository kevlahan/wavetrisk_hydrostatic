# Stage325 corrected inverse verification

The later job 55576 and independent review completed this gate for all saved
256-rank, four-step fixed level-5/tol0 captures, including all1024physical
inverse bridges. See
[native qualification](parallel-blocks-stage325-native-v5-qualification.md).
The design and earlier bounded attempt below remain historical evidence.

The next gate reuses the immutable job55569 captures. It fixes the verifier's
comparison of local inverse storage slots and recomputes the optimized pair's
physical program, value and provider evidence. It changes no Fortran source,
binary, numerical work, solver inputs or timing instrumentation.

The previous full replay55573 passed768streams and254providerfamilies perarm but
stopped at the cross-arm comparison. Two hash-bound rank101captures show six
storage slots exchanged while representing the same308physical nodes. Across four
steps, the corrected feature comparison resolves all eight reported feature
mismatches: all1596commonfeatures agree. This local result does not qualify the
complete256rankworld.

## Physical values require a proof beyond renumbering

Whole-array value hashes cannot be permuted into physical-node hashes. The audit
initially exposed24such missing direct comparisons at rank101. The existing
captures nevertheless include the required fixed-physical publications. The
corrected checker uses those publications under a closed source/schedule contract:

* Every moved slot represents an uncatalogued physical node and belongs to the
  unique consumer-alias snapshot. Moved rows cannot enter scaffold or restore
  arrays; those physical axes must remain identical.
* Only the configured level5/5 optimized final transaction is accepted. Exact
  scope order, options, flags, adjacency, zero-work counters and absence of scalar,
  inner or outer arithmetic are checked. Unsupported writes or schedules fail.
* The immutable consumer snapshot is published to fixed Domain positions. Its
  exact physical field hashes establish family1 values at phases1 and2.
* Family2 is unchanged after phase1 and is published to fixed Domain positions,
  establishing all three phases. Family1's ordinary publication establishes
  phase3 because subsequent operations do not write uncatalogued workspace rows.
* All original raw-array hashes, shapes, counts and source identities remain
  mandatory. A changed publication, source, provider, value or execution order
  cannot be accepted by merely relabeling slots.

The exact source contract pins inverse, MPI, Domain capture and shared sources;
full native metadata checks also bind all258original compiled input files. Every
spawned worker verifies actual source bytes independently. This proof applies to
this simplified configuration; it does not qualify adaptive or multilevel runs.

## Transfers and overwrite order

Inverse interior transfers can gather into workspace and scatter to blocks;
boundary/outer transfers scatter only. Certificates use physical destinations in
the actual direction. Unique block destination keys are required in the supported
non-edgewise configuration. The source emits each catalogue/boundary/row once,
level plans filter without duplication, and outer plans combine disjoint patch
and boundary families. One workspace node may supply several different block
locations without imposing an artificial read order.

Remote alias sources snapshot before local copies. Local alias chains retain
complete execution order. For remote repeated writes, positive and negative edge
codes refer to one destination; their signs remain in ordered payloads. The
adversarial review caught and corrected a potential bug that would have grouped
those as separate destinations.

## Bounded execution

Recovery-v4 uses one16CPUtask onbb13/nk. Sixteen private cache writers reconstruct
both optimized arms, with rank-strided shards. Each SQLite file and its four
ordered tables are certified. Rows are rechecked after committed merging before
source shards are deleted. Two complete provider worlds reconcile concurrently.
The checked-L256 fullstream/providerclosure from55573 is reused through its
sealed certificate and all historical input hashes; it is not evidence for the
new optimized cross-arm comparison.

The coordinator compares complete feature maps, all1024inverse bridges and all
provider certificate families. It retains both proof dictionaries and a complete
mismatch census before deciding the gate. Work and preservation are independent
reports. Aggregate scratch is sampled against4GiB and conservative aggregate
virtual memory against16GiB; individual16GiB process limits also apply. All
verifier-owned process groups are terminated on failure or alarm.

There are zero builds, solver runs and performance measurements. Native source,
compiler/runtime, binary, checkpoint, science and placement qualification remain
bound to original job55569. Collection is limited to64MiB and excludes raw captures,
checkpoints, work TSV and full peer vectors.

A300secondallocation and270secondalarm preserve the original2700secondtotal:
2392seconds were already consumed and the maximum cumulative target time is2692.
Observed parsing/provider costs make a parallel completion plausible, not assured.
Timeout is an incomplete result; no retry or budget extension follows.

## Validation before submission

All183package tests pass locally and after fresh sealed extraction. The Linux
aggregate-memory test is skipped locally and runs in bbserv's runtime preflight.
The source/identity tests include16independent adversarial bridge cases. Final
spawned-worker probing uses two retained captures in a private local binding
fixture, verifies all four actual source files in each worker, and passes the
actual compressed-read and full cache paths. Its scope is local checker validation.

The sealed archive is `rk325-verification-recovery-v4.tar.gz`, SHA256
`de10508111f8d9adba4317cb6ea887401a8c9a03e035179547a1f8dc056254bb`;
manifestSHA256 is
`0c5b5b0d5a1046d51539474a4232da83d556e66bdfd9e5c4aa005f6d7c87d030`.
The reproducer is `test/parallel_block_profile/locality_325/recovery_v4`.

A complete native pass and independent review are required for the previously
authorized conditional commit. This gate supplies no performance benefit,
accepted-source promotion or cause of the original negative scaling by itself.

## Native preflight compatibility correction

The initial v4 package failed four constructed merge tests on bbserv because its
SQLite connection did not enable file URI interpretation explicitly. Submission
was stopped; no Slurm job ran. The corrected v4p1 package reopens the merge database
with uri=True before attaching read-only shards. All183local tests and the final
actual spawned-worker fixture pass. The original v4 archive and failed native
preflight log remain preserved. The allocation budget is unchanged.

## Stage325 corrected inverse continuation — job55574,2026-10-26

The user authorized the corrected campaign. Source review and adversarial tests
establish a physical-value bridge using existing fixed-Domain publication records,
exact raw array hashes and the closed5/5optimized inverse schedule. This avoids
weakening the value requirement or requiring a new solver capture. The saved
rank101pair now passes all1596commonfeatures and four physical bridges.

Recovery-v4p1 is RUNNING from14:59:47 onbb13/nk withone16CPUtask,300secondcap,
270secondalarm and zero compiler/solver/performance runs. All183tests pass both
in fresh extraction and onbbserv. The initialv4unitpreflightcaught aSQLiteURI
compatibility issue and stopped before submission; that package/log is preserved.
V4p1 explicitly enablesURIhandling. Earlier allocation elapsed2392seconds plus
this300secondcap remains2692seconds within the original2700secondbudget.

Monitor55574to terminal, collect the bounded review and independently inspect all
optimized feature/provider/physical bridge, work and preservation evidence. The
checked native closure is reused from55573; both optimized worlds are recomputed.
Commit only if the completegatepasses. No retry, budget extension, solver or
performance campaign is automatic. See docs/parallel-blocks-stage325-corrected-inverse-gate.md.

## Stage325 corrected continuation terminal review — job55574

Job55574endedFAILED after264seconds (4:24), with a handled termination signal
while awaiting final provider outputs. All512optimized streams,16source-bound
shards and16exact merges completed. Both provider loops recorded254/254 without
a new semantic assertion failure. Final proof exports, cross-arm feature/provider
comparison,1024physical-value bridge comparisons and work replay did not finish.
This is interrupted verification, not a demonstrated new numerical mismatch.

The976,591-byte review has75verified entries and the sealed62filepayload matches.
Independent post-terminal checks confirm258sources,79originalkitfiles,3889metadata,
binaries/checkpoints and all1536capture/work hashes unchanged. All183native unit
tests passed before submission. Resource peaks are sampled and below declared
caps. Slurm logged cancellation at255seconds, before the270secondalarm; the exact
signal number was not recorded. Do not claim the270secondalarm/300secondlimit elapsed.

Hold full qualification, conditional commit, performance and promotion. The
original total now uses2656seconds (44:16), leaving44seconds. No retry or budget
extension is authorized automatically. Terminal review is complete; stop monitoring.
See docs/parallel-blocks-stage325-corrected-inverse-gate.md and
 test/parallel_block_profile/locality_325/validation/independent-review55574.json.

The subsequent offline software repair is documented separately in
`parallel-blocks-stage325-resumable-verification.md`. It preserves closed derived
phase artifacts and removes full feature/bridge JSON materialization. The
original v4p1recipe and deployed identities remain unchanged; no new native
verification or solver run is implied by the offline test pass.

The complete stream report records per-arm row counts65024route,256Domain,
66048context (including1024bridge contexts) and508840feature rows. All16worker
source attestations and all512actual compressed-read hashes were independently
checked against the historical anchors. These premises do not replace the final
pair comparison. Both provider progress reports record254/254, but the parent had
not received/exported final payloads before interruption. The sealed cleanup
removed temporary caches/payloads; original captures remain available. The bounded
review therefore cannot replay the incomplete final comparisons by itself.

This continuation's five-minute budget was too tight for the complete pipeline
including finalization and independent work. A future proposal must preserve
completed per-phase artifacts and quantify finalization cost before requesting
another budget. The present evidence establishes no performance benefit or cause
of the solver's160-to256negative scaling. No additional campaign was launched.
