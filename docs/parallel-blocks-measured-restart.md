# Restart Domain migration from measured Stage230

The user authorized a new development sequence on 2026-10-07: start from the
last fast qualified Stage230 source, remove Domain pathways progressively,
and resolve any repeatable performance penalty within the stage that introduces
it. Performance optimization cannot be deferred to later migration stages.

## Primary strategic goal: efficient execution at four or more ranks per Domain

On 2026-10-07 the user explicitly reaffirmed that the block data structure is
intended to support efficient execution on **four or more MPI ranks per geometric
Domain**, `P >= 4*D`. The limit to break is the number of geometric Domains, not
the number of blocks. Correct execution on ranks without Domains establishes
capability, not performance success. Legacy parity below the Domain count is a
necessary regression check, but cannot by itself qualify the refactoring's goal.

Prioritize measured work still confined to Domain owners, duplicated Domain and
block computation, representation transfers, global synchronization, and block
granularity/load imbalance. Each substantial implementation change must retain
correctness and numerical/oracle checks, compare against its accepted baseline
and legacy main where supported, and include performance beyond the Domain
count. Resolve a regression before stacking another migration change.

The proposed J5 83/160/256-rank comparison is an intermediate diagnostic only:
with 160 geometric Domains, the fourfold target is **640 ranks or more**. That
requires multiple nodes on this cluster, a placement policy appropriate to those
nodes, and matched comparisons within the same hardware/allocation conditions.
Use partition bb or all with bb04 excluded for such campaigns; nk/bb13 alone has
256 usable cores. Report total ordinary-step time, speedup and relative parallel
efficiency alongside per-rank work/imbalance evidence. No minimum efficiency
percentage has been agreed yet, and no fourfold scaling success is claimed.

Current root-child splitting yields at most four blocks per Domain. At `P=4*D`
that permits only one block per rank when every Domain is split, with little
freedom to balance unequal subtree work; beyond that count the present block
capacity is insufficient. Finer subdivision is therefore a strategic design
requirement to investigate alongside removing Domain-bound work, with each
change measured separately. Stage270 provides cost evidence at 83 ranks and
correctness above the Domain count; it does not qualify the fourfold goal.

## Restart baseline

The branch `codex/domain-migration-measured` starts at `4ea64012`, which removed
production Domain tendency storage and measured 0.754893 s/ordinary J5 step at
83 ranks on bb13. The existing `codex/performance-investigation` implementation
and Stage231–262 evidence remain the numerical and design reference. The later
owner-graph implementation is not imported wholesale into this branch.

Stage230 still uses Domain solution, primitive, diffusion, physics and lifecycle
consumers. Its predecessor's matched main comparison measured a 215.386 ms gap;
the Stage230 campaign itself did not run main. Stage264 therefore measures exact
Stage230 directly against the established legacy main snapshot before changes.

## Stage264: baseline campaign

`test/parallel_block_profile/reset_264` prepares a self-contained source and
script kit. It requires only the established J5 fixture on bbserv; it does not
source previous job state or import other campaign kits. Both source trees are
sealed, Stage230 numerical inputs are verified against the committed source,
and the legacy snapshot is the same main source used in Stage249–251, with only
the established checkpoint I/O hardening. The same 23 external physics objects
are used by optimized Stage230, optimized main and checked Stage230.

One exclusive nk/bb13 allocation builds serially, validates J5 at 16 and 167
ranks, and runs two independent 83-rank main/Stage230/Stage230/main batches.
Placement, fixture, arithmetic controls and ordinary-step selection match the
established campaigns. No bb04 exclusion is needed on nk. Local correctness
uses four MPI ranks, serial builds, and one thread per library. Local timings
cannot establish performance.

### Reviewed baseline: job 55374

The transferred `rk264-review-55374.tar.gz` passed independent review of 2,239
hashed files and all fourteen runs. Optimized, checked and oracle J5 validation
passed at the scheduled 16/167-rank configurations, including ranks without
Domains. Exact block checkpoint comparisons, main-relative numerical policy,
same-count VTK hashes, two reloads/remaps and production Domain-trend absence
passed. Binary/checkpoint/VTK payloads are omitted from the review archive;
their cluster-computed hashes and comparison reports were verified, not rerun.

The mean of four independent 83-rank ordinary-step run medians was 0.542412091 s
for legacy main and 0.762952294 s for Stage230. The gap is 0.220540203 s, or
40.659% above main. All four paired gaps are positive (217.753–223.546 ms).
Stage230 run medians span 0.759168296–0.768837219 s; main spans
0.540949456–0.545291611 s. This is consistent with the historical fast baseline
and confirms a substantial existing legacy gap. The historical 0.754893 s was
measured in another allocation and is not a matched candidate comparison.

This qualifies the reset baseline, not legacy performance parity. The first
change below remains the target; timings alone do not establish how much of
the 220.54 ms it can remove. Subsequent candidates must compare to this exact
source in the same allocation and retain the main anchor.

Independent review initially exposed a Python tuple/list serialization mismatch
in its storage-census comparison. Comparing JSON forms fixes the false failure
without changing solver code or relaxing storage checks. Review now passes;
a resealed archive with an altered census is still rejected. Machine-readable
evidence is in `parallel-blocks-stage264-results.json`.

## First substantial change: remove the inverse-exit SOL return transport

The target is the repeated complete SOL block-to-Domain return transaction in
`activate_block_native_inverse_transform`, immediately after
`install_inverse203`. Stage230's inverse workspace lives on the geometry rank
which also owns the required Domain consumer. The return transaction retrieves
that result from final block owners even though the producer already holds it.

The proposed replacement publishes canonical atmospheric results directly from
the native inverse workspace to the existing consumer. Keep both current block
installations during this first stage. Their later coalescing is a separate
change with separate timing evidence. Removing the return route removes four
full scalar/vector return transactions in the validated RK4 schedule; the net
gain includes the new direct publication and all downstream effects.

### Required equivalence proof

1. Publish canonical interior values from `restore203`, together with the latest
   required boundary view. `value` contains alias-overwritten interior values
   and cannot be copied indiscriminately. `publish_native_inverse` currently
   publishes boundary nodes only; it is not an interior replacement as written.
2. Preserve the provisional all-level boundary refresh and the accepted-step
   physics consumer aliases. These represent different value epochs.
3. Preflight generation, addresses, field extents, coverage and finite values
   before authoritative mutation. Reuse compact maps without a second complete
   Domain numerical image or per-value owner-program lookup.
4. Account for untouched scaffold values and non-atmospheric state. The generic
   writer also writes global `sol` outside the supplied atmospheric slice,
   including soil. Prove equivalence or provide a bounded direct replacement
   before removing that behavior.
5. Keep zero-Domain ranks in every still-required block transfer. A rank without
   geometry can own final blocks. Verify block state and boundary generations
   after every provisional and accepted inverse.
6. Retain the old return reconstruction in an independent test-only comparison
   path. Assert exact canonical/scaffold/alias/soil equality at the transaction
   boundary, then check all twelve states, restart, remap, checkpoint and VTK
   outputs against Stage230. Count the removed production return calls.

This removes an unnecessary Domain transport pathway while retaining consumers
that still need Domain-addressed inputs. Complete consumer migration proceeds
only after measured stages establish a faster or neutral replacement.

### Performance acceptance

The first candidate must pass four-rank optimized/checked/oracle local gates,
then matched bb13 J5 checks including ranks without Domains. Screen cheaply,
confirm the candidate against exact Stage230 with independent paired runs, and
include legacy main in the milestone campaign. Compare complete ordinary-step
time, including preparation, communication and publication.

A repeatable regression blocks the next stage. A mixed-sign, near-noise result
is inconclusive and receives at most one predeclared extra batch; it is not
accepted automatically under a 3% allowance. Acceptance requires evidence of a
gain or neutrality within the measured run variability, with no accumulated
increase in the main-relative gap. No automatic promotion occurs.

On 2026-10-07 the user accepted a shorter default for future campaigns: one
six-run `main, immediate baseline, candidate, candidate, immediate baseline,
main` batch initially. This gives two order-balanced candidate/baseline pairs
and two main anchors. Repeat that batch once only when the change is small,
mixed-sign or near measured variability. Preserve correctness/oracle gates.
Do not edit already-delivered kits or reinterpret their historical schedules.

Stage265 job 55381 passed all scientific checks but averaged 0.764427491 s
versus its matched Stage230 baseline's 0.760609614 s, a 3.817877 ms penalty.
Four paired differences include both signs, so performance qualification is
inconclusive. The matched main anchor is 0.540708576 s. Keep Stage265 unqualified
and investigate its new publication cost before another Domain removal.

The prepared [Stage265 tiled-publication refinement](parallel-blocks-stage265-tiled-plan.md)
changes only inverse publication relative to Stage265. Optimized, checked and
oracle four-rank replays passed exactly, including all 48 inverse-return proofs
and 62 VTK members. Its shorter campaign compares the complete refinement with
accepted Stage230 and main.

Job 55382 verified all thirteen runs and supports accepting the refinement at
83 ranks: candidate 0.753104779 s versus Stage230 0.767293997 s and main
0.542783628 s. Both paired gains are positive (15.273 and 13.106 ms); no repeat
is recommended. The remaining main gap is 210.321 ms (+38.749%), so legacy
parity remains unachieved. Recommend candidate `fab72486` as the next immediate
baseline, with main retained. See `parallel-blocks-stage265-tiled-results.json`.

The [Stage266 audit and separate campaign](parallel-blocks-stage266-audit.md)
omit the first full block installation only on provisional Stage203 production
exits. Boundary refresh/final installation and the oracle's original two-call
sequence remain. Optimized, checked and independent oracle four-rank replays
passed exactly. The campaign compares this single omission with accepted
Stage265t and main using seven validation runs and six timing runs.

Job 55398 verified all thirteen runs. Both 16/167-rank oracles retained all
48 exact-return proofs and 84 installations. At 83 ranks, candidate
0.748448779 s beats the matched Stage265t baseline 0.758481670 s by 10.033 ms
(1.323%), with both pairs improving. Recommend accepting `ff740744` as the next
immediate baseline without a repeat. Main is 0.541753610 s; the remaining gap is
206.695 ms (+38.153%). See `parallel-blocks-stage266-results.json`.

## Following stages

The provisional duplicate inverse installation is now independently measured.
The [Stage267 consumer audit](parallel-blocks-stage267-audit.md) found that the
remaining provisional canonical workspace-to-Domain interior publication is
still required by shared primitives, diffusion and native velocity gradients.
Native boundary refresh alone does not establish independence from that image.
The provisional surface-pressure reader proposal remains available, preserving
exact legacy arithmetic, coverage and boundary epochs, but its cost is not
isolated. The [Stage268 priority audit](parallel-blocks-stage268-priorities.md)
supersedes choosing that reader as the immediate next optimization. Investigate
wavelet/inverse and adaptation/reconstruction overhead first, alongside eager
hydrostatic cache construction and repeated thermodynamic recomposition. A
promising bounded omission is lazy provisional cache construction after proving
the cache has no production numerical reader. These are proposals, not measured
gains; keep full publication for the consumers identified by Stage267.

The [Stage269 bounded cost campaign](parallel-blocks-stage269-audit.md) is now
prepared. Exact Stage266/main controls and separate diagnostic observers cover
all major transaction families in one campaign, with overhead calibration for
both families. Four real observer replays and workflow checks passed locally;
Job 55466 completed and passed all scientific gates. Its [review](parallel-blocks-stage269-results.md)
measures a bare-control gap of 212.308 ms. Main observer calibration failed;
Stage266 observer calibration passed and identifies large broad tendency,
recomposition, migration and inverse-wait costs. Cache rebuilds and surface
pressure are small. Reduce observer frequency and decompose the large producer/
transport regions before choosing one numerical intervention; retain the failed
cross-version attribution as unqualified.

The [Stage270 review](parallel-blocks-stage270-results.md), job 55467, verifies all
2,540 archive members and all fifteen runs. Total Slurm elapsed time was
**00:15:38**, including five serial builds; this is recorded in the
[runtime ledger](parallel-blocks-campaign-runtimes.json) for future estimates.
Bare controls retain a 209.311 ms gap (+38.737%). Candidate calibration passes
at +2.758 ms; main still fails at +8.715 ms, so cross-version phase attribution
remains unqualified. Candidate-only decomposition identifies 193.987 ms of
shared compatibility producer self time, 87.574 ms of scalar restriction
inclusive time, and 64.341 ms of migration inclusive time. Migration includes
20.540 ms of checksum scans/proof reductions; this is not a promised saving.

Prioritize the Domain-bound producer/capture/restriction chain for the scaling
goal and audit migration proof/preparation work as an independent bounded
optimization. Include an uninstrumented 83/160/256-rank J5 baseline and matched
candidate comparisons in the next performance campaign, with legacy main where
supported. This is an intermediate screen toward at least 640 ranks, not final
qualification. Do not use the 167-rank oracle as production scaling evidence or
repeat the unchanged broad audit solely to obtain a passing main calibration.

The [Stage271 candidate and scaling campaign](parallel-blocks-stage271-plan.md)
are prepared and locally validated. One payload proof reduction of three integer
sums replaces three scalar collectives, preserving checksum scans, equality
checks, byte totals and numerical work. This bounded experiment is separate
from the larger Domain-bound producer replacement. Four builds and twenty runs
compare uninstrumented Stage266/candidate at 83/160/256 ranks, with main anchors
at 83/160 and an oracle at 256. Three exact local replays, twelve actual-routine
MPI proof cases and nine workflow tests passed; the sealed delivery was checked
from a fresh extraction. Accounting capture is included. Stage266 remains the
accepted performance baseline until the cluster results qualify or reject this
candidate; no speedup or fourfold-scaling success is claimed by preparation.

The [Stage271 review](parallel-blocks-stage271-results.md), job 55468, now verifies
all twenty runs and 5,314 archive files. Slurm elapsed was **00:21:47**. The fused
candidate is not promoted: changes are +5.164 ms at 83 ranks, −2.357 ms at 160,
and −1.194 ms at 256 with mixed pair signs at 256. The exact Stage266 solver
source has been restored; the experiment remains in commit `ab9c8463` and its
sealed kit. Do not stack further changes on the unqualified candidate.

Stage266 itself slows from 549.270 to 584.620 ms when increasing 160 → 256 ranks,
a 6.436% penalty despite 60% more ranks. Source inspection confirms shared
operators and even native mass/velocity workspaces still follow `size(grid)`
and Domain ownership. The next investigation must distinguish these producer
costs and their waits from growing communication/replicated preparation using
targeted, separately calibrated 160/256-rank phase/work evidence. Whole-process
CPU/RSS records cannot provide that attribution. Prioritize the complete shared
producer/capture/restriction ownership replacement, preserve Stage267 epochs and
oracles, and keep the >=640-rank strategic target; do not resume serial small
optimizations or another main-profiler calibration campaign.

Then migrate the remaining shared primitive consumers and reassess adaptation
and the remaining complete producer/consumer chains. Choose each change from its
actual removable work and replacement cost. Preserve formulas and operation
order where exact checks require it; use compact block-owned stencil plans and
bounded exchanges rather than rebuilding the generic owner graph by default.

Completion requires no non-essential Domain numerical dependency in the
qualified production step, full-candidate speed comparable to or better
than legacy main where supported, and measured efficient execution at four or
more MPI ranks per geometric Domain. Retain necessary
geometry and explicit output/checkpoint adapters until their own replacements
are validated and measured.

## Stage272 targeted scaling diagnostic (prepared 2026-10-08)

The next campaign uses exact accepted Stage266 plus a separate observer overlay
at 160/256 ranks. Production source remains at `ff740744`; rejected Stage271
fusion is absent. It records actual Domain/block ownership, Domain producer
patch visits, native restriction work and selected remote payload counts, with
separate local producer and MPI-wait scopes. Compare Domain owners with ranks
without Domains using recorded membership. Calibrate at both counts before
qualifying cross-count phase attribution; 83-rank calibration is insufficient.

Four serial builds and fourteen runs retain exact science, checked/oracle and
placement gates. Timing order is 160 A-B, 256 A-B-B-A, 160 B-A. Main is only a
16-rank correctness anchor in this diagnostic. Runtime estimate is 16–22 minutes
after allocation; accounting is captured automatically. See the
[Stage272 plan](parallel-blocks-stage272-plan.md) and
[commands](../test/parallel_block_profile/scaling_272/README.md).

Use the resulting evidence to choose one complete producer/consumer ownership
replacement or address measured preparation/communication growth first. MPI
scope durations combine waiting and transfers; inventories are not executed
work; no automatic speedup or causal conclusion is claimed. The strategic goal
remains efficient execution on at least 640 ranks for J5's 160 Domains.

## Stage272 review (job 55471, 2026-10-08)

All fourteen runs pass correctness; runtime is 18:52. Untouched Stage266 averages
547.507 ms at 160 ranks and 594.543 ms at 256 (+8.591%). Both preset observer
calibrations fail: +9.516 ms mean overhead at 160; a 12.105 ms control spread at
256 exceeds its 11.891 ms bound. Phase timing attribution remains unqualified.

Exact counters nevertheless confirm two limitations. Domain producer visits
remain confined to 160 owners. At 256 ranks, those owners also carry 77.21% of
block patch/divergence work, leaving the extra 96 ranks 22.79% despite being
37.5% of ranks. The groups are balanced across sockets. The increasing-rank,
fixed-target next-fit assignment can front-load work; an equal-weight 640-block,
256-rank counterexample yields 192 ranks with three blocks and 64 with one,
despite a feasible two/three-block distribution with the same maximum.

Next proposed isolated intervention: balanced contiguous block assignment,
initially preserving catalogue order and the weight definition. Qualify actual
work distribution and uninstrumented 83/160/256 performance; retain the full
Domain producer/consumer replacement as the architectural target afterward.
Do not promote from unqualified phase times or claim all capture waiting is
recoverable. See [Stage272 results](parallel-blocks-stage272-results.md).

## Stage273 isolated partition candidate (2026-10-08)

The candidate replaces next-fit assignment with an exact minimax contiguous
partition of the existing ordered leaf weights, then chooses feasible cuts near
cumulative equal-work targets. A naive nearest-target rule was insufficient:
nonuniform-weight tests showed it could increase the heaviest rank's load.
The capacity constraint guarantees the smallest possible peak for this weight
model. It does not guarantee lower complete-step time, communication, or Domain
producer cost. Catalogue order, granularity, weight formula and numerical kernels
remain unchanged; only `src/parallel_block_mpi.f90` differs from Stage266.

Actual-routine Fortran tests cover optimized/checked assignment against an
independent dynamic-programming reference. Scientific validations record actual
catalogue weights and ownership; review reconstructs the old assignment on those
same weights as an explicitly labeled counterfactual. Timing runs disable these
records and added phase profiling. Four builds/twenty runs compare accepted
Stage266 and the candidate at 83/160/256 ranks, with supported-count legacy main
anchors, checked/oracle coverage and retained exact scientific gates.

Budget 20–25 minutes after allocation, based on Stage271's measured 21:47 for
the same build/run counts. Candidate is unpromoted pending bbserv evidence;
Stage266 remains the accepted fallback. A successful partition correction still
requires complete Domain producer/consumer migration and finer granularity for
efficient execution at >=640 ranks. See the
[Stage273 plan](parallel-blocks-stage273-plan.md) and
[commands](../test/parallel_block_profile/partition_273/README.md).

## New direct block architecture program: Stage274 (2026-10-08)

The user authorized beginning a new migration program informed by the failed
owner-program approach. Stage251 confirmed a 386 ms velocity-owner penalty;
Stage252 attributed most of the observed all-step increase to source execution.
Subsequent empty-collective and traversal optimizations recovered only part of
the regression. The replacement must change the execution structure rather than
repeat that architecture with more frequent benchmarking.

Stage266 already has direct velocity kernels, but its workspaces and upstream
inputs remain on Domain owners. Stage274 exports actual source/restriction/
gradient plans from a disposable Stage266 build and compiles an ownership-cut
budget. It compares the completed tendency exported today with the inputs that
would have to move if velocity moved alone, and flags same-level remote source
dependencies and unmapped boundary/scaffold inputs. This selects a complete
producer/consumer slice before its numerical migration. It adds no numerical
interpreter, new MPI collectives or production ownership changes.

Four builds/eight full runs retain science, checked/oracle and zero-Domain-rank
coverage; 83/160/256 configurations provide real plan samples. This is not a
performance campaign. The first three plan builds per Domain are structural
samples, not a synchronized or call-weighted trajectory. Conservative packed
input budgets and conditional upstream block placement do not prove a speedup.
Stage273 remains separate and unpromoted. The architectural goal remains direct
fused block producers, bounded halo/restriction communication, retirement of the
displaced Domain computation/transfers, and efficient execution at >=640 ranks.
See [Stage274 plan](parallel-blocks-stage274-plan.md) and
[commands](../test/parallel_block_profile/ownership_274/README.md).

Stage274 v1 stopped in postprocessing after three completed J5 solver runs:
`Duplicate Domain snapshot`. The reader incorrectly treated the exporter's
rank-local `(Domain, sample ordinal)` counters as globally unique. Version 2
preserves observing-rank identity, separates unique Domain counts from observation
counts, and retains all cross-rank observations without assuming synchronized
epochs. It does not change solver/exporter source. The prepared recovery reuses
four verified builds and the first three solver outputs, revalidates their science
and contracts, then runs only the five remaining cases in a fresh directory.
Original evidence and both attempts' accounting are retained. Twenty tests and
reprocessing of 360 real local plans pass; actual cluster recovery is pending.
See `parallel-blocks-stage274-recovery.json` and the v2 recovery commands.
