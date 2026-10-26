# Stage325: locality candidate qualification and changed-work review

**Current qualification: job 55576 and independent terminal review passed the
full native gate for the saved level-5, tol=0, 256-rank, four-step case.**
Performance remains unqualified. The historical failures below are preserved;
see [the completed gate](parallel-blocks-stage325-native-v5-qualification.md)
for the later semantic, physical inverse, work and preservation closure.

Job 55569 completed both exact GCC 13.2.0/MVAPICH 4.0 builds and all three
256-rank, four-step runs in the simplified max_level=5, tol=0 configuration. Per-run
science, placement and work gates passed. The recorded optimized C0 and L256
initial and final checkpoints are bit-identical. The job then failed in its
semantic verifier. The original failed ledger, sealed kit, sources, binaries,
checkpoints and captures remain preserved.

Two verifier assumptions required source-derived corrections. RHS phases 4, 6
and 11 use two spatial components; the corrected gate checks the closed
14-phase source contract, channel capacities, coordinates and logical bindings.
Transform placement also deliberately retains owner=-1 on unowned dependency
rows. The corrected gate requires their entire sentinel tuple, requires owned
rows to name the actual receiving rank, and permutes only real rank labels.
The adapter route additionally checks that all sources are in its owned
numerical prefix, while its Domain tail remains destination-only. These are
verification changes; neither defect calls for another solver run.

Recovery 55571 stopped after 13 seconds at the second defect. The subsequent
nine-capture local preflight covers fixed and moved groups, storage-order
changes, owned scaffold rows and ranks without Domains. All 15 focused
owner/adapter tests pass independently; all 5,450 available cross-arm feature
comparisons in those captures match. This selected preflight does not replace
the full 768-stream and world-provider replay.

## What the recorded work establishes

The comparison below uses the optimized C0/L256 pair. Checked-L256 has extra
oracle calls and clones instead of optimized migration storage moves, so its
full work totals are a different validation workload.

An independent audit recomputed every returned phase summary and projected
hardware endpoint summary from all 768 compact rank files, checked their hashes,
and reconciled all 640 physical source-block extents with the exact maps.
All 16 numerical stages match under the same 256-group permutation. Fixed Domain
counts stay on physical ranks. Full native peer arrays were not included in
the compact return; their pairwise proof remains part of native replay.

| Work or capacity measure | C0 | L256 |
|---|---:|---:|
| Retained / remote blocks per reconstruction | 6 / 634 | 400 / 240 |
| Remote migration bytes per reconstruction | 408,995,792 | 155,278,016 |
| Peak rank migration send bytes | 2,580,040 | 1,964,756 |
| Peak rank migration receive bytes | 1,964,756 | 1,905,304 |
| Peak aggregate migration payload allocation | 4,544,796 | 1,964,757 |
| Imported publication rows per preparation | 59,165 | 24,240 |
| Peak transform workspace rows | 726 | 469 |
| Peak transform exchange scratch values | 69,452 | 40,832 |

There are five catalogue reconstructions and four transform/RHS preparations.
Migration serialization and matching unpack/checksum byte extents fall 62.03%.
Publication routes 9301/9303 each reduce snapshot, receive and install values
by 59.03%, with lower rank peaks and fewer remote requests. These are actual
recorded reductions in named operations and arrays, beyond a change in socket
labels. Array capacities are not peak process memory.

The reduction does not eliminate the surrounding preparation chain. All 640
source blocks are still constructed, with peak 4 per source rank. All 17 direct
migration collectives per rank/reconstruction remain. RHS geometry still
selects and packs 352 packets per preparation, with peak 4; total RHS and transform
geometry packing is unchanged. More geometry is sent to self, and those self
packets still require serialization/parsing.

## Retained costs and offsets

Boundary traffic illustrates the necessary replacement accounting. Across all recorded calls in state
family 1, remote scalar gather/install values fall 54,915,072→20,908,032 while
retained direct values rise 503,808→34,510,848. Their total remains 55,418,880;
the vector total remains 83,128,320. The same conservation holds for family 2
and for combined local/remote verification values. Wire work is removed, but
the source gathering and consumer installation remain.

Local writeback staging and reverse installation grow. Mixed input/adapter
routes gain self copies while retaining much of their gathering/snapshot and
installation work. The candidate overrides 640 ownership entries and performs
640 catalogue plus 256 permutation guard iterations per rank/reconstruction,
in addition to the original catalogue scan. Optimized block storage moves and
initial selector consensus remain unchanged.

The six numerical rank peaks are identical: nodes 348, primitive events 1,980,
scalar targets 1,950, velocity targets 1,440, exported values 85,230 and messages
4,530 per RK stage. The map changes locality and co-location of fixed and
moving work; it does not reduce any of these numerical rank peaks.

Actual placement also exposes competing effects. Transform adapter 9501/9503
cross-socket bytes increase 3.97%/3.94%; numerical transform 9001/9003
cross-socket bytes increase 15.14%/15.22%. The LLC-group velocity-work peak
rises 9,600→10,560 (+10%), although most socket numerical balances improve.
These are traffic and work exposures, not measured cache contention or elapsed
regressions. No sum of unrelated phase maxima is treated as a critical path.

## Decision and remaining question

While full replay is pending, hold performance and repair promotion. If the
complete native semantic/provider and preservation gates pass, the evidence
justifies designing a separate bounded performance evaluation of this exact
qualified locality intervention. The large measured migration/publication
reductions are a reason to test its net benefit, not a speedup forecast.

That evaluation should use the qualified optimized binary and runtime/placement
identities, the same simplified configuration, and disabled qualification
recorders. Its fixed run budget, order balance, useful-benefit threshold and
uncertainty calculation must be declared before execution. A paired C0/L256
comparison at 256 ranks evaluates this candidate's net benefit. Recovering good
160→256 scaling additionally requires an appropriately qualified, matched 160
anchor; old point estimates cannot supply that conclusion.

Even a net gain would identify a beneficial placement intervention rather than
separately attribute elapsed cost to migration, publication or hardware locality.
No gain would not by itself exonerate communication, because competing costs
can offset removed work. The main elapsed-time cause, useful 256-worker scaling
recovery and adaptive-case benefit remain unresolved until their own gates
are satisfied.

## Evidence retained

The original 6,585,848-byte review archive has SHA256
`a9df88afe7f9fce4d8421bf60e6facfee8103fb9242c33331fc38276761481e6`;
all 3,975 returned file hashes were verified. Its 768 compact native captures
total 2,554,376,142 compressed bytes and remain native. Recovery 55571's review
archive has SHA256
`8e79441572322c0df6c170a55efc07ff10f555e5c0cdf01789f154e3d8b6d065`.
The historical capture/work hashes were first recorded before 55571 verification;
they must not be described as hashes recorded before job 55569.

Reviewer-owned accounting evidence is in
`.stage325-development/accounting/native55569-accounting-audit.json`,
`native55569-review.md`, `recovery-wrapper-review.md`,
`transform-owner-review55571.md` and `transform-owner-probe55571.json`.
Final replay identity, stream/family/feature counts, resource usage, preservation
checks and independent acceptance must be appended after the result is reviewed.

## Terminal review of verification-only job 55572

Job 55572 ran on bb13/nk from 2026-10-26T13:02:59 to 13:10:55 and ended
FAILED after 476 seconds. It completed all 256 checked-L256 stream checks and
148 of 254 provider families for that arm. It failed on the next family,
the initialization ghost plan, with `KeyError('ghost_exchange_plan%local_order297')`.
The optimized arms, complete provider comparison and independent work replay
were not reached. No build, solver or performance run was launched.

This is a verifier option/allocation contract defect. All three recorded run
environments have `WAVETRISK_SHARED_GHOST_SOURCE=0`. The exact bound Fortran source
only constructs grouped local-order tables when this option is enabled. When
disabled, the audit emits the unallocated shape `[-1]` and no payload, and the
scalar consumer traverses original request-index order filtered to local owners.
The verifier unconditionally indexed the missing optional payload. All 50 ghost
plans in ten saved captures show the expected unallocated shapes without payload.
The source's original loop order is available from the captured request arrays;
ignoring local order or accepting a missing field indiscriminately would weaken
verification and is not a valid repair.

The earlier nine-stream preflight exercised per-rank parsing, coverage, context
construction and shared-feature comparison. It did not execute the complete
`route_certificates` provider branch. Existing report tests mock that function.
This was a concrete preflight gap: the valid disabled-option state was already
present in the saved captures, but the failing consumption path was not tested.
A further replay should not be proposed on the strength of those tests alone.

The 491,032-byte returned archive contains 38 hash-verified files; its SHA256 is
`78070fb6a7817d6e6565953e06b16e3240b7634868504add15771dc163ec058f`.
All 1,536 capture/work hashes matched before replay. All 256 actually consumed
compressed streams matched again during parsing. A separate read-only check after
failure revalidated the original 79 kit files, 258 source inputs, 3,889 metadata
files, both binaries, checkpoints, and every one of the 1,536 capture/work hashes.
Original solver evidence and sealed packages remain unchanged.

Peak process RSS was 220,380 KiB and sampled SQLite storage 1,103,671,296 bytes;
this was not a memory or time-limit failure. The one-CPU job used Slurm `--mem=0`
with a 16 GiB process virtual-memory limit. Combined elapsed time for jobs
55569, 55571 and 55572 is 881 seconds (14:41), leaving 1,819 seconds (30:19)
inside the original 45-minute total. No new allocation or extension is authorized
by this terminal review.

Hold native qualification, promotion and the validation-conditional follow-up
commit. The next justified work is an offline audit and repair of the provider
checker, with end-to-end tests of every remaining provider family and supported
option/allocation state before any proposed native replay. It must check actual
ungrouped request order for the recorded configuration and preserve the distinct
grouped-order contract. No retry, solver rerun or performance campaign was
submitted. The completion monitor is stopped after this terminal review.

The independent review and post-failure preservation records are saved in the
parent checkout's `.stage325-development/independent-review55572.json` and
`postfailure-preservation55572.json`. The ten-capture inspection is retained as
`ghost-contract-55572.json`. The current verifier and sealed recovery-v2 payload
remain unchanged for reproducibility of this failure.


## Subsequent offline repair

The user-authorized provider audit and repair are now completed locally.
See [the provider preflight](parallel-blocks-stage325-provider-preflight.md).
The terminal failure and original artifacts above remain preserved; the new
checker has not yet completed the full native replay.

## Stage325 corrected native provider gate — job55573,2026-10-26

The user authorized the next gate after the completed offline provider audit.
At submission, recovery-v3 was RUNNING on bb13/nk from13:56:07
(scheduler timestamp, no offset supplied). It reuses the three original capture
sets: zero builds, solver runs or performance measurements. All112 packaged tests
pass both before sealing and after fresh archive extraction. The package verifies
both failed-replay histories and the original1,536 input hashes. Its30-minute
allocation and29-minute runner alarm bound cumulative target time to2,681seconds
within the original2,700-second budget;881seconds were already consumed.

Full native semantic/provider and work qualification remains pending. Monitoring
is active every3minutes, with bounded collection and independent review on terminal
state. Commit only if the full gate supports it. Never automatically retry, extend
the budget, modify original evidence or submit a performance run.
Reproducer: test/parallel_block_profile/locality_325/recovery_v3/build_recovery325.py.
ArchiveSHA256: `39750bf8db0ac8f71a56673724adc554efa744bb74f5e47bc38177ce1374fc21`.
ManifestSHA256: `386333bfa0bb686f4dbec9cdc8f940df50d8c020dc02c76ff5d2a8e7003d70ad`.

## Terminal result of recovery-v3, job55573

The scheduler reports FAILED,13:56:07–14:21:18,elapsed1,511seconds (25:11).
No build, solver or performance run was executed. The30-minute allocation limit
and29-minute alarm were not reached. Peak RSS was408,904KiB and sampled SQLite
usage3,303,763,968bytes, within the16GiB process and4GiB scratch caps. The failure
was a comparison assertion, not resource exhaustion.

The returned archive `rk325-replay-v3-review-55573.tar.gz` is592,618bytes,
SHA256 `bf893a02891ee6a8b9d13145b962a9de94764bdbd46a7b938ef354fed82d23f2`.
All46members and the original sealed recovery-v3 package manifest verify.
The exact verifier matches the preflight version. All768compressed capture reads
match the historical55571hashes, in addition to the1,536pre-verification bindings.
Independent post-failure metadata/source/binary/checkpoint and1,536input-hash
checks pass. Both original optimized checkpoints remain byte-identical, and the
original science, placement and per-run work gates remain valid. This review did
not independently decode science again or complete the interrupted work replay.

All256streams and254providerfamilies completed in each of checked-L256,qual-C0,
and qual-L256. All five storage frames were checked, including48changed storage
groups. Per-arm resolved consumers were6,732,148 /6,734,118 /6,732,148 respectively.
Those counts do not establish cross-arm provider identity: the comparison of
provider certificates follows the failed feature comparison and was never reached.
The report retains only the first three unequal feature keys, not a full mismatch
census. The scratch database was removed under the sealed cleanup policy. Do not
claim a complete global diagnosis from the first reported rank.

### What the first mismatch establishes

Two additional native rank101captures (10,278,999compressed bytes total) were
copied locally and checked against the unchanged historical hashes. Parsing and
inspection took about2seconds locally; no native verification retry was made.
For each of four steps the only differing fields within the inverse plan and node
map are `node_id` and `slot`. The six positions exchange slots297–299 with300–302,
for physical node IDs210–212 and258–260 in local Domain1 (global Domain102).
Each arm contains the same308unique physical nodes. The maps are reciprocal,
and translating slot IDs makes the complete node-to-slot maps identical. Physical
keys and all other fields in these sampled inverse scopes are identical.
Among1,596common feature keys at this rank there are eight mismatches: these two
scopes in each of four steps. This is a narrowly established storage-renumbering
case, not proof that all inverse consumers or all ranks are equivalent.

The source file matches the original candidate manifest (SHA256
`65b936677196c97fabf0f4ec8363132c3e48e9e4ff1f17ba4e789d7781d60db0`).
`node_slot` assigns the next slot on first request at lines259–270. The boundary
manifest is received in MPI peer order at lines382–388, then calls `node_slot`.
Changing owner labels can therefore change appended slot numbering even though
Domain ownership and physical nodes are unchanged. The current checker hashes
these raw arrays as fixed-Domain features, explaining why this representation
change fails that check. This source interpretation does not license discarding
the arrays: a correct comparison must map every live reader, writer, alias,
operation, flag and value to the same physical identity while preserving execution
and overwrite order. The failed preflight did not cover this fixed-rank mismatch;
its selected captures did not include rank101 in either arm.

### Decision

Qualification remains FAILED/incomplete. Keep the solver and checker unchanged,
hold the conditional commit and performance campaign, and stop the monitor after
this review. Total target elapsed is2,392seconds (39:52), leaving308seconds (5:08)
inside the original45-minute budget. Another full25-minute replay would not fit;
there is no automatic extension or retry.

The justified next proposal is an offline audit of the inverse storage namespace
and all consumers, followed by adversarial physical-identity tests against the
saved mismatch. It must preserve distinctions between representational changes
and true operation/value/order changes and expose all mismatches without stopping
at three keys. This review has not implemented a checker repair or authorized a
new campaign. No scaling benefit or primary slowdown cause is established.

Machine-readable evidence is in locality_325/validation/independent-review55573.json,
inverse-failure55573.json and postfailure-preservation55573.json. The parent
`.stage325-development` retains the bounded review, two captured inputs and the
inspection script; full original captures/checkpoints remain on bbserv.

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
