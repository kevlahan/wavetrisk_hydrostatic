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
