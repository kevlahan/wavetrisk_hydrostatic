# Stage278 connected diffusion and pressure integration

The implementation remains based on accepted Stage266. This is local numerical
integration evidence, not a qualified block-owned RHS/RK candidate or a bbserv
campaign kit. The accepted solver and both final goals remain unchanged.

## Changes

`parallel_block_rhs_diffusion` executes complete scalar Laplacian, velocity
divergence, divergence Laplacian and rotational Laplacian hierarchies. It includes
fine-to-coarse flux restrictions, unmasked diffusion divergence, direct/SW
writes, curl, and all four rotational pentagon corrections. The kernels accept multiple independent columns with explicit initial scratch;
the real climate replay preserves one physical layer at a time and batches the
two scalar variables. Batching all physical layers is not qualified here.
The legacy diffusion phase schedule is preserved. In particular,
velocity divergence has an edge exchange but no node exchange.

`parallel_block_rhs_compile` binds these kernels to actual adaptive geometry and
ordered physical storage aliases. Diffusion divergence is compiled independently
of masked mass divergence: inactive mass stencils intentionally lack complete
edge addresses. Source Domain namespaces remain distinct when merged. The
compiler allocates its own offset array, including the final sentinel.

`parallel_block_rhs_boundary` compiles the actual signed pack/unpack topology into
physical source identities. It preserves three separate effects: remote sends
sample values before local copies; same-rank copies propagate in their original
order across **all** levels; remote receives overwrite afterward. The compact
exchange then performs the resulting parallel copies without calling the Domain
boundary routine. This avoids incorrectly snapshotting sequential alias chains.
The tested binding retains the original geometric namespaces and placements;
redistribution to final block owners remains required.

Pressure compilation records each patch's integration and producer intervals.
The new surface and per-patch pressure kernels preserve repeated visits and
intermediate lower/upper pressure states. A once-per-node column pass is not a
sufficient general ordering contract. The tested J4 surface traversal has one
visit per covered node; separate synthetic tests exercise repeated and untouched
aliases.

The existing mass restriction kernel now accepts an optional work-accounting
flag. Accepted callers retain their behavior; diffusion calls do not increment
mass-retirement counters.

## Qualification and limits

The reproducible driver starts from the preserved accepted local J4 source,
applies a validation-only overlay, and builds serially with no more than four
MPI ranks. The overlay has two modes:

- Compare new diffusion fluxes/results against every corresponding accepted
  routine, with independent scratch initialized **before** legacy computation.
  Pressure checks compare the actual intermediate values after each patch.
- Skip those legacy diffusion routines within the candidate RHS transaction and
  install compact results for the real physics consumers. The complete accepted
  tendency remains independent in the oracle run. This tests consumer effects;
  its Domain installation adapter is validation-only and earns no ownership or
  retirement credit.

Final optimized shadow, checked installed, and optimized installed/oracle runs
all pass: each covers48 RHS transactions,30 layers,12 topology generations,
12 exact timesteps, two reloads/two remaps, exact checkpoint fields and62 exact
VTK members. Each run checks 60,396,000 pressure-node visits. Per rank there are
1,440 scalar-Laplacian sweeps,1,440 divu sweeps,1,440 Laplacian-divu sweeps and
4,320 rotational level sweeps. The512-case optimized/checked unit suite also
passes, now including repeated pressure visits and untouched aliases.

In shadow comparisons, finite scratch values are checked bit-for-bit; matching
NaNs in untouched scratch are permitted. Active scientific outputs contain no
nonfinite differences. The installed cases skip legacy diffusion only in the
candidate RHS; the independent tendency oracle still uses the accepted routines.

All source identities, commands, binary hashes, exact scientific comparisons,
and per-rank operation counts are recorded in the associated evidence inventory.
Only the final corrected runs are qualification evidence. Earlier runs and
compiler failures are preserved separately. In particular, an initial harness
allocation failure and a later underallocated offset sentinel invalidate those
initial development runs, even though the latter produced matching science.
A stale Fortran module dependency error was resolved with a fresh checked build.

## Scorecard decision

These changes resolve numerical and physical-halo integration questions. They
do **not** complete final-owner metric/mask/plan placement, resident RHS production,
actual RK/inverse consumers, or paired removal of Domain computation and exports.
No ordinary production selector activates the new chain. No speedup, P>D scaling,
legacy-gap reduction, or >=4D qualification is claimed. Validation timings include
duplicate checks and are not a production performance comparison.

The campaign entry gate remains closed until the complete production change and
its replacement costs are qualified. Stage266 remains accepted. No cluster job
has been submitted, and no further user guidance or approval is needed for the
remaining implementation. Stopping at the earlier component milestones was an
execution shortfall, not an unresolved user decision.

Evidence: [source-bound inventory](../.stage278-validation/diffusion-validation.json).
