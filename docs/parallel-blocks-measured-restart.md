# Restart Domain migration from measured Stage230

The user authorized a new development sequence on 2026-10-07: start from the
last fast qualified Stage230 source, remove Domain pathways progressively,
and resolve any repeatable performance penalty within the stage that introduces
it. Performance optimization cannot be deferred to later migration stages.

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

## Following stages

After the first result, consider coalescing duplicate inverse installation only
if no intervening consumer needs the first view. Then reassess adaptation and
the remaining complete producer/consumer chains. Choose each change from its
actual removable work and replacement cost. Preserve formulas and operation
order where exact checks require it; use compact block-owned stencil plans and
bounded exchanges rather than rebuilding the generic owner graph by default.

Completion requires no non-essential Domain numerical dependency in the
qualified production step and full-candidate speed comparable to or better
than legacy main across the accepted rank configurations. Retain necessary
geometry and explicit output/checkpoint adapters until their own replacements
are validated and measured.
