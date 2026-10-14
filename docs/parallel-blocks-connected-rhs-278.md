# Connected block-owned RHS contract — Stage278 development

Current completion status: the connected production selector and final local
qualification are complete; the sealed bbserv kit is ready for user submission.
See [the production result and run commands](parallel-blocks-rhs-production-278-results.md)
for current evidence and remaining goal gaps. The design/reassessment text below
is the chronological pre-implementation record; its draft-readiness statements
are superseded by that result. Accepted source remains266; no cluster promotion.


Prepared 2026-10-09 from accepted Stage266 and Stage274–277 evidence. This is
an implementation contract, not a new profiling campaign or an accepted cutover.
The isolated implementation checkout is
`/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic`, based on
`ff740744608731a05d53771c32a3055309f3a586`. Stage273 and the Stage277 cache are
excluded. The original checkout, branch, uncommitted documents and validation
artifacts are preserved. The **single authoritative scorecard** remains
`/Users/kevlahan/.codex/worktrees/26d8/wavetrisk_hydrostatic/docs/parallel-blocks-scorecard.md`.

## Execution boundary and placement decision

Move the connected mass/velocity transaction, including its prerequisites and
resident RK outputs. Use a local slot space for owned nodes plus compact required
halo/alias slots; compile stencils and peer routes per topology/ownership epoch.
Run direct fused kernels and explicit level phases. Never schedule individual
arithmetic expressions by dependency depth. No ordinary-stage input packet may
carry all geometry, state and reference values as Stage276's validation packets
did. No new whole-primitive return to Domain numerical consumers is authorized.

An owned block slot is authoritative only for its compiled owned node/edge.
Ghosts and boundary/scaffold aliases have explicit producer identity and phase;
they cannot be assigned the geometric Domain owner merely because an address
is unmapped. Keep required residual scaffold computation until its reader and
replacement are accounted for. Stage275's sampled absence of same-level source
cross-owner dependencies does not prove this for every adaptive topology.

## Fields, readers, lifetime and replacement cost

The table covers the supported restarted compressible, non-mode-split RK4
climate path with Simple physics and without NCAR topography/SSO. External
mode-specific readers remain outside this cutover and must stay gated.

| Field / current producer and owner | Production readers and required order | Replacement owner / epochs | New cost and work retired at activation |
|---|---|---|---|
| State mass, temperature, velocity; means in `domain_stage` / `sol_mean`, with resident block state already installed | `cal_surf_press`, pressure integration, `step1`, scalar/velocity diffusion, velocity gradient, RK base; boundary/scaffold and inverse consumers | Final block owners; stage state and mean epochs, topology, masks, ownership and completed state halo | Reuse resident state/means. Compile local and ghost gathers; refresh only required state halos. Retire RHS-driven interior Domain state publication only after all affected readers move. Existing inverse/lifecycle publication is separately retained. |
| Surface pressure and lower pressure/geopotential; `cal_surf_press` on Domain owners | `integrate_pressure_up`, then Bernoulli and Exner; Simple physics has a separate `cal_surf_press(sol)` transaction | Column's block owner; current state and surface geometry epoch; bottom-to-top physical order | Local column scratch and surface geometry prepared per topology. Retire RHS `cal_surf_press` interior pass and RHS pressure integration. Do not remove Simple physics' separate transaction by assuming its cached pressure has the same epoch. |
| Pressure, layer Exner, upper/lower geopotential; `integrate_pressure_up` in `basic_operators` | Compressible Bernoulli uses midpoint geopotential. Final Exner is injected fine-to-coarse before velocity gradient. Pressure/Exner are also used by Simple physics, physical diagnostics and VTK in their own transactions | Block-owned RHS column scratch; current RK state, physical layer, grid-level phase; separate raw and restricted Exner lifetime | Preserve sum-then-gravity surface-pressure order and midpoint pressure arithmetic. New halo traffic for required raw/restricted values only. Retire displaced RHS interior arithmetic; keep justified diagnostic/I/O adapters without recurring bulk numerical return. |
| Mass flux, QE/PV, KE and vorticity; `ops::step1` and `post_step1` on Domain owners | Native mass restriction/divergence; `compute_native_velocity_source`/Qperp; Bernoulli; `cal_Laplacian_vector_rot` and climate velocity diffusion | Block owners; state halo, masks, topology/ownership; ordered direct and SW writes followed by pentagon correction before readers | Geometry/stencils installed once per topology/ownership epoch. Compact level halos for flux/PV; final vorticity writer routes. Retire corresponding `step1` interior assignments and primitive capture/publication. Stage277's 19 assignment oracle is validation only. |
| Bernoulli / Exner fine-to-coarse injection; `cpt_or_restr_scalar` / `scalar_cpt_restr` on Domain owners | `compute_native_velocity_gradient` consumes FINAL restricted fields after complete level sweep | Destination block owner; ordered child injection, parent mask, fine-level completion and restricted halo epoch | Route injection writes by identity in established child order. New B/Exner halo slots and messages; retire Domain injection for moved interior. Raw pressure Exner must remain distinct if column/physics readers still need it. |
| Scalar Laplacian and flux residual; `cal_Laplacian_scalars` and climate `physics_scalar_flux_case` on Domain owners | Mass and temperature flux formation, native temperature closure; order2 reads Laplacian, order1 reads perturbation state (not state+mean) | Block owners; state halo, per-layer diffusion sweep, edge restriction then divergence then node halo | Reuse flux-restriction/divergence kernels with dedicated diffusion storage. Preserve boundary labels970/971/972 semantics. Retire displaced scalar diffusion passes and residual transport. A retained adapter must declare scalar residual payload, frequency and remaining Domain work; it is not implemented by this draft. |
| Divu and its Laplacian; `cal_divu_ml`, `cal_Laplacian_divu` on Domain owners | Climate `physics_velo_source_case`; order2 selects Laplacian; order1 selects `grid%divu` | Block owners; current velocity halo, diffusion order, restriction and node halo phase | Dedicated divu/diffusion workspace and compact halos; preserve pentagon-neighborhood viscosity reduction. Retire corresponding Domain passes and physics-source gathers. |
| Final vorticity / rot(rot(vorticity)); `post_step1`, `cal_Laplacian_vector_rot`, `post_vort` on Domain owners | Climate velocity-source curl; order2 overwrites vorticity with rot(rot(vorticity)) | Block owners; last shared-triangle writer, pentagon corrections, then rotational diffusion epoch | Ordered shared writes, triangle/curl halos, dedicated intermediate vector. Preserve boundary973 semantics. Retire corresponding Domain numerical passes; direct producer completion alone is insufficient. |
| Velocity physics source; climate callback on Domain owners | Direct source combines Qperp and length-scaled physics; restricted sources follow ordered mixed-mask action tape | Destination block owner; divu/rotu completion and per-layer/grid-level source phase | Local physics arithmetic and local precomputed viscosity coefficients (including pentagon reduction), followed by source halos. Retire callback traversal/gathers and global temporary `C_visc` mutation for this path. |
| Mass RHS; `native_mass(size(grid))` with fine-node/edge exchange | Coarse flux restriction and block scalar publication; `RK_sub_step_compatibility` | Destination block owner; completed restricted edge halo, active node mask, layer/grid-level RHS epoch | Local mass tendency storage and fine-node halo. Retire Domain native mass numerical workspace, boundary exchanges and `publish_native_mass_rhs` interior payload together. Existing mass publication substitutes RHS into the physics channel and must stop when the resident consumer reads its own RHS. |
| Velocity source and RHS; `native_velocity(size(grid))` source/restriction/gradient | Coarse source restriction; `prepare_block_native_velocity_remainder` / tendency export; `RK_sub_step_compatibility` | Destination block owner; completed mass-flux/PV/physics phase then final B/Exner halo | Local ordered source/restriction and gradient kernels. Retire Domain velocity numerical workspaces, native tendency export and corresponding remainder capture. No Domain-fed full flux/PV/B/Exner packet. |
| Mass/velocity RK closure; `RK_sub_step_compatibility` plus resident block multistage accumulation | Next substage's shared RHS, block candidate stage, inverse and publication, eventual final solution | Resident block owner; base state and matching completed RHS transaction; output stage distinct from saved base | Direct resident `base+h*RHS` with existing stage accumulation semantics. Retire duplicated geometry-owner closure/materialization for moved interior. Preserve temperature boundary preparation/application and oracle's final-stage alias behavior until their interfaces move. |

## Exact phase schedule

1. Finish the resident stage state/mean halo. Bind topology/ownership/mask epochs;
   no cached plan survives a changed address map or changed numerical mask.
2. Sum physical column mass in increasing layer order, then multiply by gravity
   and add `p_top`. Integrate each column upward with midpoint pressure, retaining
   raw Exner and midpoint geopotential. Simple physics and diagnostics retain
   separate transactions; they cannot share the RHS cache without matching state.
3. For each physical layer, prepare scalar diffusion and divu diffusion. For
   each geometric level from fine to coarse, finish the previous fine-node RHS
   exchange before a coarser mass restriction reads it.
4. Run pressure/primitive work in the established patch/node order. Preserve
   direct and SW edge orientation, all shared-vorticity writes and pentagon
   corrections. Inject B/Exner in established child order; restrict mass flux.
5. Finish the restricted edge halo, compute mass divergence, and launch/finish
   fine-node exchange as required by the next level. Rotational diffusion uses
   final shared-vorticity values, including its order2 overwrite semantics.
6. Evaluate velocity physics and direct/source restriction in the existing
   action order using this exact level's restricted mass flux. Finish source
   routes before a coarse action reads a fine source. A detected same-level
   cross-owner dependency needs an explicit bounded phase split, not an interpreter.
7. After the complete level sweep, finish FINAL B/Exner halos and compute velocity
   gradient. Mark the layer's mass/velocity RHS ready only now.
8. Consume authoritative resident RHS in actual RK stage/accumulator code. Never
   send it to a geometry owner just to reconstruct an interior RK stage. Keep
   boundary/scaffold and inverse publications until their independent readers move.

## Costs and retirement checks

A source plan accounts for owned/halo slots, static stencil/metric bytes, rebuild
allocation, transient column storage, restricted fields, diffusion intermediates,
RHS and saved RK base/output. Every exchanged payload is identified by field,
producer, consumer, layer, grid-level phase and topology/ownership/state epoch.
Count preparation and all new exchanges in complete-step timings. Do not use
Stage276 reference packets or Stage277 duplicate arithmetic in bare timings.

For each field f and phase e, report remote bytes as
`8 * sum(unique required (producer slot, consumer owner, component) for f,e)`;
report messages, peer counts and reused/local slots separately. Never treat this
structural byte formula as measured performance. Stage274/275 already establish
that Domain-fed velocity operands cost over3x the current RHS export for sampled
budgets; that is why upstream producer ownership is part of this cutover.

The activation diff must pair removal of each interior legacy computation with
its replaced transfer: mass publication, velocity export, duplicate RK closure,
primitive inputs and RHS-driven state publication. A residual adapter is named
and budgeted by exact field/phase; zero-Domain ranks must execute authoritative
production work. Removing `Domain` names or calling a local kernel is insufficient.

## Initial implementation status and planned cutover

The development draft is `src/parallel_block_rhs.f90`. It adds compressible column
pressure, fused mapped direct/SW primitives and all four mapped `post_step1`
pentagon corrections, ordered B/Exner injection, climate
scalar/velocity diffusion, existing mass/source/gradient operator entry points
and a two-field RK closure. It uses no geometry-owner cache, per-value owner
interpreter, reference stream or MPI collective. The default build compiles the
module, but **no production caller selects it**. No old work has been retired.

Numerical development tests run512 nonuniform connected cases in optimized and
checked builds, with an independently extracted accepted `step1` primitive
and SW/pentagon oracles, pressure/diffusion formula checks, final-writer checks and malformed
extent/stencil rejection. Source/gradient operators are reused accepted kernels,
so their connectivity checks are not new independent numerical oracles. The
fixture, adaptive topology, reload/remap and complete timestep are not qualified
by these tests. Evidence is stored in `connected_rhs_278/local-results.json`.

Still required before an independently activatable production slice exists:

- Compile/remap owned, halo, boundary and scaffold addresses into actual resident
  block workspaces. Connect the now-tested direct/SW and `post_step1` pentagon
  arithmetic to actual producer traversal, aliases and final-writer delivery.
  North/east positive-boundary mapping and rotational `post_vort` remain required;
  the new arithmetic does not establish those routes or their epochs.
- Connect full diffusion sweeps, scalar injection, level halo routes and field
  epoch readiness to these kernels. The pointwise diffusion functions do not
  constitute the Laplacian/rotational hierarchy implementation.
- Connect actual `prepare_block_multistage_tendency` and RK stage consumers;
  retire the matched Domain numerical paths and transports in that selector.
- Qualify that entire selectable slice immediately: optimized/checked/oracle,
  exact accepted-family science, main-relative tolerances, adaptation/reload/remap,
  zero-Domain production, and bare accepted/main83/160 plus accepted256 timings.
  Resolve repeatable regression before stacking the next slice.

The [boundary development result](parallel-blocks-rhs-boundary-278-results.md)
records the additional arithmetic evidence and its remaining limitations.

No Stage278 campaign is sealed or launchable yet. Local numerical tests cannot
justify cluster timing or source promotion. Once the connected selector and
retirement assertions exist, reuse the proven sealed user-run bbserv workflow,
with immutable fixtures and original failure evidence preserved. The eventual
J5>=640 multi-node target and finer subdivision remain separate unqualified work.

## Phase exchange implementation status

The [phase exchange development result](parallel-blocks-rhs-exchange-278-results.md)
adds a topology-time distributed final-writer compiler, direct bounded peer
exchanges and schedule-to-writer metadata export. Optimized/checked synthetic MPI
and large-message tests pass. Real alias/order distribution, transformed field
halos, full numerical phases and authoritative RK consumers remain unconnected.
This development is not an independently activatable production slice and does
not change the campaign entry decision or either dual-goal score.

## Diffusion and pressure integration update

See [the integration result](parallel-blocks-rhs-diffusion-278-results.md). Full
diffusion sweeps, actual ordered physical halo binding, and per-patch pressure
scheduling now pass local complete-step tests, including authoritative diffusion
substitution in the validation overlay. This replaces the earlier pointwise-only
diffusion status. Final block-owner placement, mask-epoch integration and the
complete production RHS/RK cutover with paired retirement remain required.
