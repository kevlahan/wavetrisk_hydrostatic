# Stage267: provisional inverse publication consumer audit

Status: source audit and next-change proposal only. No solver change, new timing
measurement or performance kit is included. The immediate baseline is Stage266
`ff740744608731a05d53771c32a3055309f3a586`, supported at 83 ranks by job 55398.
Its ordinary-step mean of run medians was 0.748448779 s; same-allocation main
was 0.541753610 s, leaving 206.695 ms (+38.153%). Keep that baseline and main
`6a245db7f1632b0205c2bebe95e5f0cd535c2582` in the next comparison.

## Decision

Do not skip the remaining provisional canonical interior publication. The
immediate native boundary refresh does not re-import Domain interiors, but the
next production tendency still consumes provisional Domain mass, temperature
and velocity. Removing the publication would leave these consumers stale.
No complete field-wide omission is justified by this audit. This does not rule
out future smaller node-range omissions after their readers are proven absent.

The next proposed change is one consumer migration: evaluate provisional
compressible surface pressure from the completed inverse's local canonical
mass view. Retain the complete publication for the other readers. This removes
one numerical dependency before attempting to remove its shared adapter.
Its runtime benefit is a hypothesis to test, not a measured saving.

## Producer and next-stage sequence

`activate_block_native_inverse_transform` in `src/parallel_block_mpi.f90:14682`
publishes canonical scaling, completes the all-level scalar/vector boundary
epoch, publishes those boundaries and performs the retained final installation.
It marks Domain boundaries current, calls `end_inverse203`, and records a
provisional completion. Stage266 defers only the earlier installation in
production; its independent oracle retains both.

`native_provisional_WT` in `src/time_integr.f90:224` then calls
`refresh_parallel_block_candidate_boundary_state` with `native_inverse=.true.`.
The production closure consumes the completion marker, refreshes block ghosts,
and invalidates/rebuilds block hydrostatic state without importing Domain
interiors. Oracle comparisons still read the Domain image.

The next RK4 tendency, however, passes `q1` to
`prepare_block_multistage_tendency` (`src/time_integr.f90:140`), which calls
`block_tendency_compatibility_ml` in production. RK3 has the same dependency.
Replacing Domain tendency outputs did not replace its shared prognostic inputs.

| Subsequent consumer | Actual production dependency | Publication decision |
|---|---|---|
| Native candidate boundary refresh (`parallel_block_mpi.f90:11385`) | Resident blocks, ghost exchange and derived hydrostatic cache; Domain comparisons are oracle-only in the native closure | This closure alone does not require interior publication |
| Surface pressure (`diagnostics.f90:313`, `column_mass:344`) | Provisional Domain mass plus fixed mean mass; incompressible fallback also reads temperature | Keep until a separate native reader is qualified |
| Shared primitives (`multi_level.f90:1032`, `ops.f90:35`) | `basic_operators` binds provisional mass, temperature and velocity; `integrate_pressure_up` and `step1` consume them for pressure, Exner, Bernoulli, fluxes and potential vorticity | Keep all three fields |
| Scalar diffusion (`multi_level.f90:1592`) | Provisional Domain scalars when `Laplace_sclr==2`; callback stencils also use scalar inputs | Keep scalar interiors and boundary views |
| Velocity divergence (`multi_level.f90:2068`) | Provisional Domain velocity when `Laplace_divu/=0`; derived divergence feeds viscosity | Keep velocity interiors and boundary views |
| Native velocity gradient (`multi_level.f90:1011`) | Explicitly forms `rho=mean_mass+q_mass` and `rho_theta=mean_temp+q_temp` from full Domain ranges | Keep mass and temperature; native tendency destination does not eliminate these inputs |
| Velocity source callback (`multi_level.f90:964`) | Binds Domain mass/velocity for the callback interface; the current climate callback (`test_case_module.f90:202`) actually reads derived divergence/vorticity and geometry | Do not count pointer binding alone as a direct mass/velocity read; its derived inputs still depend on preceding primitives/diffusion |
| Independent complete tendency and inverse comparisons | Domain image, original installation sequence and exact-return checks | Retain unchanged |

The diffusion conditions matter to the established climate configuration:
`src/shared.f90:380` sets scalar and divergence Laplacian order to two, and the
climate driver does not enable its commented divergence-order override.
The shared primitive and explicit gradient reads independently block omission
even in configurations without diffusion.

Several apparent dependencies need a different classification. The scalar
divergence capture's interior Domain oracle is conditional on validation
(`parallel_block_mpi.f90:21286`); production captures native products and needed
boundary/restriction inputs. `prepare_block_native_velocity_remainder`
(`parallel_block_mpi.f90:17715`) transports the prepared native velocity tendency;
its Domain argument does not by itself prove another production prognostic read.
Conversely, an absent Domain tendency output argument does not prove native
prognostic inputs in `basic_operators`.

`RK_sub_step_compatibility` (`time_integr.f90:517`) reads the timestep-base `q`,
not the preceding provisional `q1`, for its low-storage RK formula and scaffold.
Do not attribute those base reads to provisional publication. Later native
temperature publication in `retain_block_native_multistage_candidate`
(`parallel_block_mpi.f90:29938`) supplies the *new* provisional transform; it
does not repair an earlier tendency that consumed stale `q1`.

## Coverage and lifetime constraints for the proposed pressure reader

The existing canonical publisher (`parallel_block_inverse.f90:1409`) resolves
catalogued interiors through `restore203` or the solution workspace. It also
restores uncatalogued patch scaffold from timestep-base global `sol`, checks
finite sources and destination extents, and preserves the complete return
contract. Boundary publication (`parallel_block_inverse.f90:1365`) reads a
different, refreshed workspace epoch. These sources cannot be interchanged.

`end_inverse203` (`parallel_block_inverse.f90:1731`) clears `resident203`; retained
allocation is not a valid consumer lifetime contract. The next
`prepare_native_inverse` (`parallel_block_inverse.f90:202`) reseeds data and can
replace storage when the generation changes. A proposed read-only view needs
an explicit completed-provisional token covering generation, producing RK stage
and final boundary epoch. Invalidate it before overwrite, accepted-step physics,
grid/adaptation changes or reset. An unexpected missing/stale token in the
eligible path must fail rather than silently read provisional Domain interiors.
First-stage, accepted/output, fallback and incompressible paths remain explicit
legacy paths.

Compile integer source addresses on the geometry consumer's rank, using the
already local canonical inverse workspace. Cover canonical interiors, the
refreshed boundary slots, timestep-base scaffold, and every other location
visited by the original pressure iterator. Do not assume `publication_patch`
covers the pressure stencil: it covers live catalogued interiors only. Fixed
mean fields and necessary geometry/derived pressure adapters may remain.
Use bounded reusable scratch if needed; do not introduce a complete numerical
Domain mirror, per-value generic owner lookups or a new final-owner round trip.

Preserve the legacy calculation exactly:

- Reset surface geopotential with the existing geometry operation and initialize
  the complete surface-pressure range to zero.
- Preserve each node's visits in increasing vertical level order, including any
  repeated patch/alias visits. `cal_surf_press` iterates patches from 3, passing
  `p-1` and `(st,en)=(0,1)` to `apply_onescale_to_patch`
  (`domain_ops.f90:352`). Its plus-side boundary coverage and inner-boundary
  suppression cannot be replaced by an interior-only `PATCH_SIZE**2` traversal.
  Do not introduce a deleted-patch skip absent from this caller.
- Form `mean_mass + canonical_mass`, add it to the running column sum, then
  multiply the complete sum by gravity and add `p_top` once. Preserve untouched
  nodes as well as populated columns; copy the result to `press_lower`.

The resident block hydrostatic `surface_pressure` is not an exact substitute.
`compute_local_block_hydrostatic_patch` (`parallel_block.f90:5220`) starts with
`p_top` and adds `gravity*rho` for each layer. Legacy multiplies once after
summing; the rounding sequence differs. The cache also covers block interiors
on final owners rather than the complete legacy geometry/boundary traversal.
Do not change shared hydrostatic arithmetic or add a pressure export as part of
this one migration. Likewise, retain `integrate_pressure_up`, diffusion, shared
`step1`, native gradients, the full `q1` publication and accepted I/O adapters.

## Separate implementation and qualification plan

Implement only the provisional pressure reader and its explicit lifetime and
coverage contract. Carry RK context to the shared tendency caller so the new
path consumes the preceding provisional stage, including the final tendency
stage. First-stage pressure retains its existing source. The new path must also
execute under the independent oracle and compare with a separate legacy
calculation; selecting legacy pressure whenever the oracle is enabled would
leave the new consumer untested.

Before cluster timing, require optimized, checked and independent-oracle J4
replays against exact Stage266. Verify pressure, `press_lower`, surface
geopotential and column-source bits at every eligible consumer entry. The
twelve-step RK4 fixture has 36 such entries; RK3 has two per step. Verify
coverage/visit multiplicities including plus-side boundaries, generation and
stage invalidation, and a test-only isolated pressure call with poisoned Domain
interiors to prove the new reader does not use them. Restore the compatibility
image before other readers run. Retain all existing scientific, boundary,
checkpoint/remap and VTK checks, the oracle's 48 return proofs and 84 actual
installations for the established RK4 fixture. Build serially, use at most four
local MPI ranks and one library thread, and do not overlap compilation with MPI.

Use the established short campaign on one exclusive `nk`/bb13 allocation,
without a bb04 exclusion. Build main, Stage266, optimized candidate and checked
candidate serially with the same 23 physics objects. Preserve seven correctness
runs: main/baseline/candidate/checked/oracle at 16 ranks and candidate/oracle at
167 ranks, including ranks with zero geometric Domains. Run six unprofiled
83-rank timing cases: `main, Stage266, candidate, candidate, Stage266, main`.
Compare ordinary-step run medians excluding steps 1, 8 and 9, both
order-balanced baseline/candidate pairs and the same-allocation main gap.

Optional cost diagnostics must be separate from promotion timing. Determine
whether compiled traversal/locality savings exceed view lookup and scratch
costs; removing this reader does not remove the retained publication's copies.
Repeat only for small/noisy/mixed evidence. A repeatable regression blocks
promotion and further removal. After qualification, audit the next shared
primitive consumer before attempting another publication omission.

Source hashes and indexed evidence are recorded in
`parallel-blocks-stage267-audit.json`. This audit supplies no new performance
claim and does not qualify a pressure implementation that has not been written.
