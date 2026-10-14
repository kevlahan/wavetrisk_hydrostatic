# Next repair: resident transform closure and finer numerical ownership

2026-10-14. **Original implementation specification; the provisional slice is
now implemented and locally science-qualified.** See the
[implementation and qualification](parallel-blocks-stage279-resident-transform.md).
Cluster performance remains unqualified. Base: job55504 v7, commit `de558ed03798f721991f38eeecc00f00eafc8370`.
Stage266 remains accepted. No new cluster campaign is justified by this document
alone. The existing numerical and performance gates remain in force.

## Selected next production change

Close the connected **resident RK → forward wavelets → compression → inverse →
next RHS** path on numerical block owners. Include the boundary/scaffold versions
that its consumers actually read. The first independently activatable slice is
the three provisional RK4 transitions per timestep. Keep the final transition's
adaptation/remap/physics contract explicit, then qualify its closure separately.
Do not treat the provisional slice as complete whole-timestep independence.

This selection addresses a source-confirmed dependency on Domain owners in a
material measured region. It is not a proposal to move inverse arithmetic and
return all results to unchanged Domain consumers. Vector wavelet production,
compression, inverse operations and the next-RHS input contract must close as
one slice before a production performance claim.

Finer, cost-aware numerical ownership is the separate scaling companion. Evaluate
it independently of the transform cutover. The current block count and attached
scaffold work impose a concrete limitation even if the placement algorithm improves.

## Quantitative case and its limits

[Derived budget](parallel-blocks-resident-transform-repair-case.json) records exact
source hashes, region IDs and the job55504 results hash. The profile groups below
sum disjoint self regions from the same slowest whole-window rank on each ordinary
step. They are broad screening envelopes, **not measured excess over main or
removable costs**. In particular, boundary/publication regions include callers
outside the selected transform slice. Bare gaps use mean run medians; profile
budgets use means, with known instrumentation perturbation.

| Quantity, ms/ordinary step | 83 tasks | 160 tasks | 256 tasks |
|---|---:|---:|---:|
| V7 bare time | 632.99 | 505.52 | 544.53 |
| Matched main bare time | 540.90 | 365.13 | — |
| Remaining cut for main parity | 92.09 (14.55%) | 140.40 (27.77%) | — |
| Inverse including preparation/staging | 99.26 | 69.26 | 76.68 |
| Domain boundary/publication | 59.94 | 41.37 | 31.08 |
| Native wavelet/compression work | 48.91 | 31.62 | 44.80 |
| Broad transform/boundary/publication envelope | 208.11 | 142.25 | 152.56 |
| RHS preparation | 62.74 | 58.56 | 72.52 |
| RHS exchange | 129.16 | 130.66 | 145.05 |

The broad transform envelope is material at both83 and160. Removing replicated
representations and ownership transit has a plausible absolute-performance case,
not only a high-rank scaling case. Moving its numerical producers off the fixed
Domain ranks also addresses an architectural barrier that RHS-only tuning leaves.
However, there is no calibrated timing prediction for the replacement yet.

At160 the entire broad envelope is about142ms against a140ms remaining gap.
Necessary transforms, compression and dependencies cannot all disappear, and the
provisional slice is smaller still. **Do not promise full legacy parity from this
repair alone.** The performance plan must retain RHS execution/preparation as a
separate remaining target, selected after the new attribution shows which work
and dependencies survive. This is more defensible than assuming all MPI wait is
removable transport cost or choosing an arbitrary percentage saving.

The256 scaling companion has independent evidence: mean per-step maximum/mean
primitive work is1.989, and unchanged block/scaffold bundles permit at most5.49%
reduction of mean peak primitive work from reassignment alone. Splitting heavy
units is necessary for substantial flexibility. Neither number predicts full-step
speedup. At640 ranks, the current640-block inventory averages one indivisible
unit/rank and cannot flexibly balance unequal work.

## Exact current pathways to replace

| Current source pathway | Current cost/dependency | Replacement and retirement evidence |
|---|---|---|
| `validate_complete_block_vector_wavelets` in `parallel_block_mpi.f90` | Provisional block values reconstructed into Domain-indexed vector staging; numerical production follows that staging | Compile compact vector-transform kernels with geometry metadata delivered to numerical owners; produce authoritative block wavelets directly. Count zero Domain-owner vector-transform operations for the enabled slice. |
| `prepare_block_native_wavelet_compression` / `activate_block_native_wavelet_compression` | Vector coefficients imported into blocks; compressed family written back to Domain storage | Keep both coefficient families resident from producer through compression and inverse. Retire the selected vector import and compressed-family interior publication after all selected readers move. |
| `prepare_native_inverse` / `compile203` in `parallel_block_inverse.f90` | Plans and execution arrays built around local `grid`; initial scalar/vector solution and wavelet gathers assemble the geometry-owned workspace | Bind compact inverse operations and required operands to block owners; geometry owners supply metadata, not a numerical workspace. Reuse existing direct all-column numerical kernels/order. |
| `scalar_inverse203`, `inner_inverse203`, `native_inverse_outer` | Scalar, inner-vector and outer-vector work uses geometry-owned arrays/aliases | Execute each owned output once on its numerical owner, with phase-specific direct dependency routes; retain pentagon, shared-edge and coarse/fine ordering. |
| `install_inverse203` | Canonical interiors transferred from the inverse workspace to final block owners; boundary views installed separately | Commit already resident canonical results locally. Exchange only actual external reader dependencies. Count removal of the corresponding interior installation traffic. |
| `publish_canonical_inverse203`, `publish_native_inverse` and next-stage adapters | Canonical Domain mirror plus distinct boundary/consumer aliases | Next RHS and provisional transform read resident canonical and explicit alias versions; publish only remaining named reader footprints. Retire full provisional canonical mirrors only after their reader census is empty. |

Do not double-credit prior repairs. Stage203 already suppresses repeated in-phase
interior gather/scatter in production; Stage266 already defers the first
provisional installation. Stage265 canonical Domain publication is a **local
copy**, not an extant full-SOL network return. The current audit reports zero full
SOL returns per inverse. The new retirement ledger must distinguish initial
inverse gathers, final block installation, coefficient transfers, local copies
and remaining boundary routes. It must not claim to remove costs already absent.

Stage229 measured45.7ms benefit at83 from the compact inverse workspace and batching
on its then-current baseline. That is historical evidence that locality and direct
column execution matter, not another saving available today. Preserve those
properties while changing ownership; do not reactivate the slow per-value or
per-dependency-depth executor from the earlier rejected designs.

## Producer, consumer and epoch contract

| Field/version | Producer and new owner | Production consumers | Validity and lifetime |
|---|---|---|---|
| RK canonical candidate | Existing resident RK, block owner | Forward transform and required stencil readers | RK stage/state epoch; timestep-start rollback image retained until existing acceptance point |
| Forward scalar/vector coefficients | Complete compact transform kernels, block owner | Compression, inverse and later adaptation where applicable | Topology/mask/ownership epoch plus RK stage; preserve inactive structural zeros and parent-level range |
| Compressed coefficients | Existing native compression, block owner | Inverse | Same stage, compression complete before readers; no Domain round trip |
| Inverse canonical solution | Direct scalar/outer/inner kernels, target owner | Next RHS, next provisional transform, final lifecycle consumers | Coarse-to-fine phase and canonical state epoch; an alias overwrite must not alter canonical values |
| Inverse dependency aliases | Bound producer values, required reader owners | Subsequent inverse phases | Signed orientation and exact phase version; sends snapshot before ordered local copies, remote installs follow |
| Consumer-visible aliases | Explicit captured boundary/scaffold versions | RHS boundary adapters and retained final physics readers | Preserve existing post-prolongation/pre-inner or later version as appropriate; never substitute the freshest kernel alias indiscriminately |
| Fixed scaffold inputs | Explicit current source contract, then declared numerical owner | Inverse and boundary consumers | Preserve the distinction between timestep `sol` and provisional scaling; no implicit import of the wrong epoch |

The selected slice must include its actual next-RHS reader closure. Adaptation,
remap, split physics and output are named consumers, not permission to delete
publication early. Final-stage materialization remains until those readers move.
Checkpoint savings are secondary: the user prioritizes ordinary steps, with saves
about everyO(10^3) and remaps about every five steps.

## Granularity and assignment companion

Use the existing catalog order and identities as a control. Capture per-block
operation vectors from compiled plans, including internal patches, scalar/vector
inverse operations, RHS producers, scaffold assignments and boundary/import
footprints. Add retained work by rank/phase to the placement model. Do not infer
block costs from rank RSS, wait time or leaf weights alone.

Compare three distinct decisions offline: partition the existing units with
current weights; use revised work weights; split heavy units and then repartition.
Use the Stage273 minimax reference for its applicable contiguous-partition problem,
with its limitations explicit. It was locally qualified and unpromoted in the
reviewed evidence, not demonstrated as a current-cluster speedup.

Subdivision must partition *all* numerical work. Splitting a parent subtree into
children also leaves its internal/coarse patch and scaffold targets to assign;
do not silently drop them, duplicate them, or attach the entire residual to the
first child. Require a unique owner for every target and the necessary interlevel
reader routes. Prefer selective heavy-subtree refinement with a bounded minimum
useful-work size over an arbitrary uniform block multiplier. The exact cuts and
size limits need the real per-block work/footprint record; current55504 does not
contain that record. This bounded validation capture resolves a named design
choice and does not warrant a standalone broad performance campaign.

Test ownership closure and subdivision separately. Improved load balance cannot
compensate for added representation/plan costs at83/160, and a faster160 run can
make the160→256 relative-efficiency problem worse. Both absolute times and the
candidate's own scaling curve must pass their respective gates.

## Storage, communication and implementation boundaries

- Keep compact contiguous column workspaces and prebound direct stencils. Compile
  from immutable geometry and current masks/ownership; do not access remote `grid`
  or perform per-value owner lookup during numerical execution.
- Separate topology/ownership, mask and state/alias validity. Rebuild only when the
  relevant identities change. No reuse based on matching counts, and no unproven
  cross-epoch cache: the previous exact cache probe had zero hits.
- Give routes one backend/epoch communicator context with bounded tags; drain all
  requests before release/replacement. No communicator per block, field or column.
- Bound persistent storage by locally owned numerical units plus required imports,
  phase routes and reader snapshots. Retire the displaced geometry-owned numerical
  workspace for the selected slice; retain geometry only where still needed.
- Record construction peaks and old/new overlap as well as steady state. Current
  narrow numerical-buffer growth does not explain the43.79→61.30GiB summed-RSS
  increase from160→256. Do not claim this repair eliminates all unexplained memory.
- Preserve compact boundary-only version snapshots; do not retain full alternative
  solution images or unbounded epoch histories to simplify correctness.

Explicit shared-memory windows are not selected. `MPI_Comm_split_type` alone
neither shares buffers nor removes numerical ownership dependencies. A future
shared immutable-geometry proposal needs measured duplication and NUMA/first-touch
benefits, not a supposition that the current waits are shared-memory transport.

## Qualification and decision to run another campaign

1. Establish an exact local reference for every migrated transform phase and its
   consumer-visible versions, including repeated/signed aliases, pentagons,
   untouched/frozen entries and parent/coarse boundaries. Use the actual compact
   kernels and preserve arithmetic order. Serial builds, at most four local ranks.
2. Complete optimized, checked and independent-oracle full replays with the enabled
   production slice. Require zero displaced Domain numerical work and explicit
   retired traffic/copy counts. Above-Domain execution needs cluster qualification;
   local four-rank speed is not evidence for83/160/256 performance.
3. Account for preparation, imported operands, routes, temporary peaks and steady
   storage. A new ownership path that merely substitutes another full-size mirror
   or unbounded metadata compiler is not campaign-ready.
4. Compare locally qualified changes against exact v7,266 and main on matched
   bb13 placement. Keep ownership-only and granularity-only effects separable;
   measure a combination only after individual correctness/non-regression gates.
   Report repeated bare ordinary complete-step results and remap/non-remap groups,
   memory and a bounded coarse profile. Do not repeat an inconclusive small gain
   unless the uncertainty changes the implementation decision.
5. Require83/160 absolute improvement without repeatable regression, and improvement
   of256 against v7. Scaling progress also requires an improved160→256 speedup
   using the **same new candidate** at both counts. A256 time below the new160 time
   crosses the current negative-speedup boundary; it does not by itself establish
   acceptable final efficiency. No final efficiency threshold has been agreed.
6. Main parity is earned only by the matched full-step comparison. Required cuts
   are presently92.09/140.40ms. The next implementation does not get called a
   full-gap repair simply for reducing a named region. If the remaining gap lies
   in RHS preparation/execution, select its next change from the residual evidence;
   do not assume batching or subtract all waiting as removable network latency.

The source evidence selects this architectural repair, but existing profiles do
not establish a numerical speedup forecast. Implementation, paired retirement,
local qualification and measured replacement costs are required before preparing
its next runnable kit. This specification itself is not a campaign-ready claim.

## J6/J8 and the original target

Keep the current J5/D=160 workload and83/160/256 comparisons as the deciding
performance track, then qualify matched160/320/640 on bb hardware when ownership
and granularity permit. Do not splice bb13 and bb measurements into one curve.

The user's possible future J6 coarsest/J8 finest workload is complementary. With
four times the Domains, D becomes640:160 and256 tasks are both below the Domain
count, and640 tasks only reaches it. The original P>=4D objective would require
at least2,560 tasks there. That larger-work test may reveal useful strong-scaling
headroom, but cannot replace demonstrating scalable block execution above the
Domain count on the current workload. No larger workload campaign is scheduled.
