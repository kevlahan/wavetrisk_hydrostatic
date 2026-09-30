# Stage231 handoff

Stage230 committed: `4ea64012`. Main untouched at `6a245db7`.
The user authorized the remaining Domain-state migration and explicitly wants
to return to optimization and explain poor scaling with more blocks/tasks.
See `parallel-blocks-stage231.md` and the source-anchored dependency JSON.

Current uncommitted solver changes are in parallel_block.f90,
parallel_block_mpi.f90 and time_integr.f90. Production retained tendency assembly
has no Domain solution argument; oracle/non-native modes require reference
state. Production no longer allocates or refreshes the unused block hydrostatic
cache; cache creation fails when disabled. Explicit validation retains it.
The live Domain primitive, RK and wavelet consumers remain. Do not describe
this as completed Domain-free dynamics or a demonstrated performance gain.

Raw build/run evidence is kept under `.stage230-validation/stage231/`.
The `candidate` tree and early production-4 were cache-only preliminary work;
its remaining tests were stopped after the combined candidate was ready.
The authoritative new build tree is `final-candidate`, with build-final.log and
build-final-checked.log; `final-*` runs and final-check-state.json track its
validation against the existing Stage230 production-4/41 references and fixed
local fixture. No new cluster campaign is prepared or submitted.

Load-balance question: leaf-cell weighting with ordered next-fit is approximate;
the initial 5% overload tolerance is relaxed by a factor of 1.2 as necessary.
It excludes communication/internal-level work and cannot distribute the
remaining Domain-owner phases to zero-Domain ranks. At one block per rank,
unequal block weights require finer splitting, not a different owner map.
Distinguish catalogue weight balance from measured complete-timestep balance.

Final local status: all four full cases passed (optimized production 4/41,
optimized full oracles 4, checked production 4), with exact same-rank checkpoint,
step-state and 62-member VTK comparisons. Both negative guards passed, and the
one-step mixed validation control passed. Review with review_local_231.py;
source/binary/log identities are in parallel-blocks-stage231-local-results.json.
The preliminary cache-only tests remain superseded. No performance claim.

User sequencing: remove remaining Domain-owned numerical work first; return to
block load balancing and further splitting afterwards. The next implementation
unit is the complete primitive/diffusion -> mass/velocity -> retained-tendency
path currently entered through block_tendency_compatibility_ml. Preserve ordered
level/column evaluation and canonical aliases, execute on final block owners,
and reject Domain numerical reads in that path. Then migrate RK/wavelet stage
state, remaining supported physics/soil and lifecycle consumers, and delete
their last Domain arrays/transfers. Each migrated boundary needs intermediate
oracle checks plus full ordinary-step, adaptation/remap and restart coverage;
success is removal of the old production dependency, not just matching output
while retaining a fallback. Geometry metadata and bounded I/O adapters must be
distinguished from persistent Domain numerical state.

Step 1 implementation is now started as Stage232 but **not delivered**.
See parallel-blocks-stage232-work-in-progress.md. New standalone arithmetic
kernels in src/parallel_block_primitive_kernels.f90 import only kind_mod and
are not yet in the production build. Test-only comparison adapters/builders
are under test/parallel_block_profile/primitive_232. The 10-group one-step
checked comparison and pressure-corruption control passed. A separate long
checked/all-oracle comparison is in .stage232-validation/checked-oracles-4.
Do not confuse arithmetic extraction with moving execution or supplying inputs
from final-owner resident state. The production owner/operand/version plan,
diffusion and restriction deliveries, direct retained installation and the
explicit temporary RK return remain to implement. Do not send this as a
bbserv campaign or restore full Domain-shaped snapshots to close the route.

Stage232 local validation is now complete for the arithmetic component: the
long checked/all-oracle run passed 12 steps, two remaps, exact same-rank
checkpoint/step-state and 62 VTK-member comparisons; the ten-group short check
and corrupted-Exner rejection also passed. The durable review is
parallel-blocks-stage232-local-results.json. Step 1 is explicitly incomplete.
No production owner migration, Domain-read rejection or scaling result is
established by these arithmetic tests. No new commit or cluster kit was made.


Stage233 connects interstage pressure production to resident inputs and final
block owners, with an explicit compatibility result return. See
parallel-blocks-stage233-owner-columns.md and its source-bound local result
summary. BLOCK_COLUMNS defaults off; check enables an independent bitwise
pressure oracle. Missing/stale operands and forbidden Domain solution input are
rejected. This is a real runtime integration, but **not the complete Step1
producer-to-tendency deliverable**: stage-one lifetimes, reference aliases,
regular/SW/pentagon consumers, diffusion and mass/velocity scheduling still need
migration. Do not run a performance campaign on this intermediate return route.

Stage234 continuation: read `parallel-blocks-handoff-234.md` for the current
owner primitive/gradient implementation, direct velocity retention and remaining
full-tendency work. Earlier stage descriptions above are historical.
