# Stage270: reduced-frequency producer and transport cost audit

Preparation follows the [verified Stage269 review](parallel-blocks-stage269-results.md).
The production source remains accepted Stage266 `ff740744`; legacy main remains
`6a245db7` with its established checkpoint I/O hardening. Only separate diagnostic
copies change. No numerical omission, ownership change or speedup is established
by preparing this audit.

## Why this test

Stage269 verified all fifteen runs and measured a bare-control main gap of
212.308 ms per ordinary 83-rank step. Its main observer added 21.147 ms and failed
calibration; candidate calibration passed at +6.929 ms, near its 7.512 ms limit.
Cross-version phase attribution therefore remains unqualified. Candidate-only
leads were substantial unsplit tendency work (~164 ms), recomposition (~109 ms),
migration (~64 ms), inverse waits (~53 ms) and scalar forward transforms (~52 ms).
Cache rebuilding (~2.1 ms) and surface pressure (~1.5 ms) are lower priorities
for recovering this gap.

Those large labels cover several different operations. Native recomposition
includes scalar restriction and its exchanges before the final physical kernel.
Migration includes packing, installation and plan/workspace preparation. A
replacement must address actual redundant work rather than optimize an assumed
meaning of a profiler label.

## Observer changes

Remove all clocks from Domain boundary routines and per-patch native velocity
source/gradient routines. Stage269 counted about 3,150 boundary scopes per step
in main and 2,419 in the candidate, plus candidate velocity scopes. Their clock
cost is a plausible reason for the failed main calibration; this campaign
measures whether reducing frequency resolves it. It does not subtract an assumed
overhead correction or relax the overhead thresholds.

Add bounded whole-routine and existing horizontal-level batch scopes, covering
all physical layers at once. There are no node/physical-layer clocks, added
synchronization or timed-step output. The fixed record schema has 60 slots;
three retired slots must remain completely empty. Per-rank records flush once
at the end of the trajectory. Self-time conservation, stage transactions,
retired slots and plan accounting are independently checked.

The new decomposition covers:

| Large region | Separately measured operations |
|---|---|
| Tendency preparation | Shared compatibility producer; mass/velocity plan preparation; scalar capture preparation/publication/MPI; velocity remainder preparation; RK assembly |
| Scalar restriction | Local flux/restriction/divergence batches; ghost packing and same-rank copies; MPI posting/completion; received installation; temperature edges and waits |
| Final physical recomposition | Metric preparation; whole physical kernel batch; accumulator retention; any Domain-to-block payload |
| Migration/reconstruction | Manifest/check; packing batch; size/payload exchange; unpacking batch; complete-store installation; ghost/writeback plan preparation; production workspace preparation |
| Payload exchange | Payload MPI; byte-checksum scans; proof reductions |

The checksum split follows a specific source-backed lead:
`exchange_block_migration_payloads` unconditionally scans outgoing and incoming
serialized bytes into integer checksums and performs three proof reductions.
These currently remain part of the production transaction. No checks are removed
here. Price their scans and collectives separately before proposing an omission
or integrating checks into an existing traversal. Small unused per-vector
`difference_moment` arithmetic remains a lower-priority hypothesis until the
major kernel/transport costs are isolated.

Retain prior inverse/cache counters and add actual ghost/writeback/metric plan
rebuild/reuse events, migration logical bytes and scalar-ghost/temperature remote
request/payload counts. These counters measure named work, not physical memory
bandwidth or every solver message. A ready-plan call is counted as reuse only
when it actually returns through the existing reuse branch. No lifecycle flag,
plan reuse policy, field value or exchange schedule changes. Rebuild counters
do not prove topology was unchanged: any proposed reuse must compare actual
geometry, ownership and payload epochs, including partially changed blocks.

## Campaign

The self-contained kit bundles four sealed source families: untouched Stage266,
untouched main, block observer and main observer. It requires only the established
J5 fixture on bbserv, not old campaign/result directories. One exclusive nk/bb13
allocation makes five serial builds with the same 23 physics objects, including
a checked block observer. Partition nk contains only bb13.

Seven correctness runs cover both controls and optimized observers at 16 ranks,
checked block observer, and block oracles at 16/167 ranks. Keep checkpoint,
trajectory, 62 VTK members, restart/remap, storage, physical placement and
zero-Domain-rank gates. Each oracle retains 48 exact inverse proofs and 84
installations. Eight 83-rank runs use the same fixed calibration order as Stage269:
`main, main-observer, main-observer, main, Stage266, block-observer, block-observer, Stage266`.
There are fifteen solver runs total; both controls and both optimized observers
have two independent timing samples. Ordinary steps exclude steps 1, 8 and 9.

The unchanged prespecified calibration limits are max(3 ms, 1%) mean shift,
max(5 ms, 2%) per paired shift and 2% spread. Both observers must pass before
cross-version attribution informs priority. Failures remain visible and science
checks continue to apply. No automatic repeat or promotion is offered.

## Interpretation and next decision

New child scopes change the meaning of parent **self** times. Use parent inclusive
time for the complete transaction and self time for each newly isolated component.
Do not compare Stage270 parent self times directly with Stage269's less-decomposed
parent self times. Do not sum overlapping inclusive scopes. Signed self-region
mean differences remain additive within the new schema, but they do not equal
the difference of whole-run ordinary medians. Longest-local-step composition
selects one rank per step and does not establish a synchronized global critical
path. CPU includes MPI spinning; transport waits can reflect upstream imbalance.

Review the calibrated current decomposition and actual rebuild/payload counts.
Select one substantial, demonstrably redundant operation with a bounded complete
replacement; preserve numerical order, accepted boundary epochs and independent
oracles. Qualify the single intervention with the short unprofiled main/Stage266/
candidate comparison before accepting it or stacking another Domain removal.
The Stage267 shared-reader constraints still apply. No speedup forecast is
justified by this diagnostic campaign alone.

The [kit instructions](../test/parallel_block_profile/cost_audit_270/README.md)
provide no-edit transfer, launch, status and collection commands. The launcher
records a unique result path and submitted job ID, rejects duplicates until
collection, and handles purged IDs and early failures. The reviewer uses stable
float sums plus the tested eight-ulp rule for derived report comparison; raw
hashes, integer counts and scientific acceptance remain exact. Reviews omit
source/binary/checkpoint/VTK payloads and verify cluster-computed hashes and
reports without rerunning the omitted payloads locally.

## Completed local validation

Three new builds passed serially with `make -j1`. Four real J4 replays ran one
at a time with four MPI ranks and one library thread: optimized main observer,
optimized block observer, checked block observer and independent block oracle.
All match their source-bound controls exactly, including checkpoint, step states
and all 62 VTK members. The oracle retains all 48 inverse-return proofs and 84
installations. Every rank has complete 12 × 60 records, conserved self times,
empty retired slots and consistent rebuild/reuse accounting.

Ordinary-step scope entries fell from Stage269's 3,173.222 to 23 in main and
3,101.889 to 470.444 in optimized/checked block runs, reductions of 99.275% and
84.834%. The block oracle records 585.444 entries. These are matched local J4
work counts, not bb13 overhead measurements or speedup estimates. The optimized
block records two ghost/writeback plan calls per step, one rebuild and one reuse,
and four edge-metric calls, one rebuild and three reuses. Reuse/epoch validity
remains the existing implementation; this audit changes no preparation policy.

Twelve workflow tests cover the complete fifteen-run/five-build schedule,
duplicate/failed-job handling, collection and independent review, calibration
failure and signed offsets, raw conservation and retired-slot rejection,
source schema and removed hooks, multiline batch anchoring and report rounding.
A real Fortran test at O0/O2 verifies nested accounting, disabled output and
invalid selector/scope/stack/step states. Fresh delivery extraction repeats
these tests and checks exact controls, helper imports, shell syntax, retired
hook absence and absence of macOS extended-attribute tar headers.

The [sealed preparation record](parallel-blocks-stage270-preparation.json)
contains the archive identity, local evidence and limitations. Current bb13
calibration and recoverable performance gains remain pending.
