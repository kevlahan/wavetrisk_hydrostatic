# Stage233: first live connection from resident operands to block-owned primitives

Stage234 subsequently connects regular/SW/pentagon producers and final gradients
to these owner results; see `parallel-blocks-stage234-owner-primitives.md`.
The results below describe the frozen Stage233 implementation.

## Status and scope

This implements the **interstage hydrostatic column producer** on final block
owners. It does **not** complete the requested full tendency migration. The
feature is opt-in (`BLOCK_COLUMNS=on`); the normal build remains `off`.
`BLOCK_COLUMNS=check` executes the independent old pressure arithmetic and
compares its outputs bit for bit with the routed block-owner results before
installing those results. Both enabled variants require `STAGE203=on`.

This is an implementation/correctness step, not a new performance campaign.
There is no bbserv kit, speedup claim or change to `main`.

## What is connected

`parallel_block_columns_mod` compiles the actual surface-pressure and
pressure-integration visits from the existing patch topology. Integer metadata
retains physical storage identity and visit counts. Boundary/scaffold aliases
are not coalesced by coordinates. Each alias column is assigned to the final
owner of the first catalogue block rooted in its parent Domain; this explicit
closure assignment is not a load-balancing claim.

Interior column operations are assigned to their catalogue block owner. At
execution they read mass, temperature and means directly from `Block_Data`,
using its declared vertical layout. The pressure/geopotential/Exner kernel from
Stage232 is called on that owner. The surface sum preserves layer order and
repeated visits, then multiplies by gravity and adds top pressure, matching the
original arithmetic. No full Domain solution snapshot is constructed.

For aliases, the native-inverse plan registers only the pressure consumer
closure and exposes the completed provisional-inverse epoch. The owner receives
these explicitly requested scalar operands. Missing slots, nonfinite operands,
stale epochs, wrong phase order and unconsumed/duplicate results cause failures.
The current implementation still imports **alias reference means** and the
case-specific **surface-geopotential reference** from geometry owners at the
start of each transaction. Those are remaining Domain reference dependencies,
not removed state.

A column carries its lower pressure and geopotential across vertical layers and
horizontal-level phases. Per-phase result deliveries contain the five values
actually consumed by the existing primitive path. Consumers advance through
ordered visits. Production `cal_surf_press()` receives **no Domain solution
argument**; passing one is rejected. Production upward integration returns the
owner-computed results before any mass/temperature pointer read.

The compatibility return is real remaining traffic: Domain pressure,
geopotential and Exner storage still exists because the next primitive and
velocity-gradient consumers have not moved. It must be removed when those
consumers execute on the owner. The current per-phase collectives and result
buffers are not presented as a performance improvement.

## Lifecycle boundary

This route handles RK stages after a completed **provisional** inverse. Stage
one is deliberately still handled by the existing producer: restart/bootstrap,
final-transform consumer aliases, remapping and external physics have different
state lifetimes. The new epoch guard rejects use of those versions as if they
were provisional operands. Plan changes follow the inverse topology generation;
changing generation during a transaction is an error.

The remaining full-tendency work is concrete:

1. Supply the initial/final-transform and post-physics column input lifetimes,
   and migrate alias reference means/surface reference ownership.
2. Connect regular, southwest and pentagon primitive consumers directly to the
   owner column results, retaining the regular-to-southwest PV context and the
   earlier vorticity writes consumed by pentagon corrections.
3. Move diffusion prepasses, ordered restrictions, mass/velocity source and
   gradient execution to final owners, with exact intermediate comparisons.
4. Install tendencies directly into retained block storage. Keep only the
   explicit temporary return needed by the still-unmigrated RK consumer.
5. Run the complete cycle with all former numerical inputs unavailable before
   deleting the old allocation/capture/return pathways. Load balancing and
   further block splitting follow that dependency removal.

## Verified local results

The final frozen build passed seven checks:

- Checked independent pressure comparisons at 4 ranks: 12 steps, two remaps,
  restart coverage, exact checkpoint fields and step states, and 62 identical
  VTK payload members against the same-rank Stage230 reference.
- Optimized production at 41 ranks / 40 Domains: the same full numerical gates
  passed. Rank 40 owned no Domains and executed 531,270 pressure-layer operations;
  it handled both compact interior and assigned alias columns after adaptation.
- The checked feature-off control completed its step with the same reported
  physical state and no owner-column transactions.
- All four injected faults were rejected: missing operand, stale inverse epoch,
  supplied Domain solution argument, and corrupted delivered Exner value.

These operation counts demonstrate execution coverage, not independent samples
or expected performance savings. The complete tendency migration remains open.

## Reproduction and evidence

Use `test/parallel_block_profile/columns_233/build_columns.py` with a fresh
output path and the frozen Stage231 source tree. It builds an optimized owner
path, a checked independent-comparison path, and a checked feature-off control.
The build manifest records changed sources, complete Fortran source hashes,
commands and executable identities. Reusing an output tree requires all changed
sources still to match; use fresh output paths for solver runs.

`run_columns.py` runs the evolved fixture. Full runs additionally use the
established Stage230 checkpoint, step-state and VTK comparisons, with a reference
at the same rank count. For more than 40 ranks it requires a plan showing actual
column ownership on a rank with no geometric Domains. `reject_invalid.py` builds
isolated replacement objects to test missing operands, stale epochs, forbidden
Domain solution inputs and corrupted transported results. It never overwrites
the candidate executable or source.

The source-bound result summary is `parallel-blocks-stage233-local-results.json`.
Raw local evidence is under `.stage233-validation/`; local timings are not
performance or scaling evidence. Earlier failed build/layout checks and earlier
executables remain there for inspection and are not counted as qualification.
