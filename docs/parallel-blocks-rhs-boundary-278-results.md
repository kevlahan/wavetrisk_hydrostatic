# Stage278 boundary arithmetic development result

Reviewed 2026-10-09 in the isolated accepted266-based implementation checkout.
This is a numerical component result, not an accepted production slice.

## Concrete changes

`Rhs_Primitive_Stencil` now carries optional southwest mass, two additional
velocity addresses, three explicit edge destinations, two triangle destinations
and twelve geometry values. The direct producer executes `comp_SW` immediately
at the same node before the next producer. It reuses that node's corner PV and
six dual velocities from the fused direct kernel; it does not attempt to recover
corner PV from averaged stored QE. Its physics residual is a required per-event
input for enabled SW actions. Missing input and malformed destinations are
rejected before output mutation; disabled SW actions do not read their addresses.

`Rhs_Pentagon_Stencil` and `rhs_pentagons` implement all four established corner
cases5/6/7/8 from accepted `post_step1`. They use compiled local mass/velocity,
PV and vorticity slots and geometry arrays. Reads and writes retain legacy
statement order, including circulation reconstructed from previous vorticity,
combined-area sums, paired triangle assignments and corrected edge PV. The
caller must supply completed producer/final-writer halos before this phase.

No geometry-owner cache, ownership interpreter or ordinary-stage reference
packet was added. Production still has no caller for this numerical draft.

## Verification and evidence

Both optimized and checked serial builds pass512 nonuniform connected trials.
Each build checks direct and SW flux/PV/vorticity against original accepted
`step1`/`comp_SW` expressions extracted into the independent oracle. Tests include
SW overwriting direct outputs, final shared-triangle writes, disabled SW sentinels,
and nonzero oriented physics residuals.

Each build additionally checks all four pentagon cases on512 nonuniform inputs:
2,048 individual corrections and512 ordered four-correction sweeps. Full PV and
vorticity fields are compared bit-for-bit against original accepted `post_step1`
expressions, including untouched outputs and later reads of earlier corrections.
Both modes pass the six expected malformed-input rejection cases. Existing
pressure/restriction/diffusion/source/gradient/RK component checks still pass.

The project's strict checked Makefile module target passes (`-Werror`, bounds,
initialization and floating-point traps). This is module compilation, not a full
climate binary build. The established unused physics-include override for this
standalone module target remains explicit in earlier build notes.

Source and binary hashes, extracted oracle, build logs, successful logs and
expected failure logs are retained in
`.stage278-validation/boundary-numerical-units` and
`test/parallel_block_profile/connected_rhs_278/local-results.json`.
Earlier unit and SW-only evidence remains intact. Source identity is accepted266
plus the documented inactive draft overlay; Stage273 and the Stage277 cache are
excluded.

## Dual-goal evaluation

This closes the missing SW and `post_step1` pentagon **arithmetic** in the draft.
It does not qualify actual boundary traversal, resident ghost addresses,
final-writer delivery, `post_vort`, complete diffusion sweeps, mask/topology
invalidation, RK solver integration or retirement of any Domain work/transfers.
There is no full-step science, zero-Domain production or performance evidence.
No speedup or efficiency estimate follows from these numerical tests.

The next connected implementation dependency is the real resident plan/route
compiler and its producer/halo epochs, followed by complete diffusion and the
actual multistage RHS/RK cutover. Test that whole independently selectable
transaction immediately against accepted/main controls before promotion.
Stage266 stays accepted, and no Stage278 cluster kit is ready.
