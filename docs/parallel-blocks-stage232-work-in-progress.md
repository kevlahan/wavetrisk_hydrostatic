# Stage232: primitive arithmetic port — step 1 remains incomplete

The requested deliverable is the complete production tendency path on final
block owners, with its former Domain numerical inputs unavailable. **This
working tree does not yet deliver that path.** Production still calls
`block_tendency_compatibility_ml`, and its primitives/diffusion and native
mass/velocity execution still use geometry owners.

## Implemented component

`src/parallel_block_primitive_kernels.f90` provides explicit-operand numerical
kernels for compressible pressure/geopotential/Exner, regular-node primitives,
southwest primitives, four pentagon corrections, surface pressure and climate
scalar/velocity diffusion callbacks. The module imports only `kind_mod`. It has
no Domain types, field pointers, callbacks, MPI, allocation or saved numerical
state. At the Stage232 checkpoint it was not included in the production build.
Stage233 now connects its pressure kernel to an opt-in interstage block-owner
producer; see `parallel-blocks-stage233-owner-columns.md`. The other primitive
kernels are still awaiting owner/operand integration. This is not a hydrostatic
diagnostic cache.

The regular primitive returns the two PV intermediates consumed by its
southwest operation. Those values cannot be replaced by arbitrary edge-averaged
PV. Its compact slot order also differs from the existing nine-node velocity
operator layout; use explicit address mappings rather than reshaping one into
the other. Pentagon operations receive the prior vorticity explicitly. The caller
must supply the correct preceding writes and preserve physical storage aliases;
geometrically coincident points are not sufficient to identify an operand.
The compact slots are C,N,E,S,W,NE,unused,SW (0..7). Existing native mass and
velocity numerical kernels are also already storage-independent; their route
compilers and runtime workspaces remain geometry-owned.

The pressure routines retain the Domain arithmetic order: accumulate mass,
then multiply by gravity and add top pressure; integrate with the average of
lower and upper pressure. They do not substitute the different thermodynamic
recomposition formula or the removed validation cache.

## Local numerical checks

`test/parallel_block_profile/primitive_232/build_primitive_check.py` creates an
isolated, explicitly test-only checked build from the Stage231 frozen local
tree. The original expressions continue producing the reference values; new
kernels compute comparisons from explicit operands. The adapter has Domain
reads because it is a reference comparison. It is not the production input
route, and must not be presented as a forbidden-Domain-read proof.

The one-step adaptive-fixture check passed all ten groups exactly: regular and
southwest primitives, all four pentagon corrections, pressure integration,
surface pressure, scalar diffusion and velocity diffusion. Counts include
initialization and all executed stage calculations, not independent samples.
The complete checked/all-oracle trajectory also passed: 12 steps, two remaps,
two checkpoint initializations, exact final checkpoint fields/step states and
62 identical VTK payload members against the same-rank Stage230 reference.
That longer build checked the seven regular/southwest/pentagon/pressure groups;
the ten-group short build additionally checked surface pressure and the climate
diffusion callbacks. Deliberately corrupting Exner in a separate test-only
executable was rejected at the pressure comparison, as expected.

`parallel-blocks-stage232-local-results.json` records all results and frozen
source/binary/input/log identities. Recheck without rerunning the solver using
`python3 test/parallel_block_profile/primitive_232/review_primitive_check.py`.

Reproduce the complete comparison build and short check with:

```bash
python3 test/parallel_block_profile/primitive_232/build_primitive_check.py \
  --reference .stage230-validation/stage231/final-candidate \
  --out /private/tmp/stage232-primitive-check
python3 test/parallel_block_profile/primitive_232/run_short_check.py \
  --tree /private/tmp/stage232-primitive-check \
  --fixture .stage230-validation/20260929/fixture \
  --out /private/tmp/stage232-primitive-one-step
python3 test/parallel_block_profile/primitive_232/reject_corruption.py \
  --tree /private/tmp/stage232-primitive-check \
  --fixture .stage230-validation/20260929/fixture \
  --out /private/tmp/stage232-primitive-corruption
```

Use fresh output directories. These are local correctness tests, not a bbserv
kit or performance campaign. No production migration or scaling claim follows.

## Required work to finish the requested deliverable

1. Compile the actual primitive/diffusion operations and required operand
   closure from topology and masks. Assign operations and their numerical
   state to final block owners. Each operand needs its storage/alias identity,
   field/component, vertical layer, RK stage and producer version.
2. Supply those operands from resident block state and native-inverse aliases.
   Extend needed alias coverage explicitly, not by recreating a full Domain
   solution snapshot. Reject missing operands. Reference fields/coefficients
   need explicit lifecycle handling when the grid or vertical state changes.
3. Execute diffusion prepasses and the ordered primitive/restriction/source
   phases on final owners, using the new kernels and existing native
   mass/velocity arithmetic. The local copies in the old communication routines
   can cover all levels even during a single-level exchange; preserve that
   dependency/version contract when compiling new deliveries.
4. Install the results directly into retained block tendencies. Remove their
   old source-owner capture and return routes once no longer consumed.
5. Until step 2 of the larger strategy migrates RK, explicitly return exactly
   the tendencies consumed by `RK_sub_step_compatibility`. Its mass and velocity
   readers cover the boundary/scaffold range as well as interiors. Count this
   return as real remaining work; do not claim the resident RK cycle is closed.
6. Run the complete production path with old numerical inputs unallocated or
   poisoned, retain independent intermediate comparisons, and test adaptive
   remapping/restarts and ranks with blocks but no Domains. The present
   arithmetic comparison does not satisfy this acceptance condition.

Do not launch a performance campaign or mark step 1 complete at the arithmetic
extraction milestone. Block balancing and further splitting remain deferred
until the Domain-owned numerical work has actually migrated. Main is unchanged.

The surface-pressure kernel accepts an ordered list of mass contributions.
Its direct fixture comparison covers live patch interiors. Boundary/scaffold
slots with repeated visits still require their explicit contribution/write
sequence in the new scheduler; this helper does not infer that sequence.
