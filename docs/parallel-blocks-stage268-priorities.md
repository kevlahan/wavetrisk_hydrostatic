# Stage268: prioritize recoverable performance costs

Status: source and historical-evidence audit; proposals below are not implemented
or timed. The unvalidated surface-pressure prototype has been saved outside the
production source in `.stage268-validation/unvalidated-pressure-prototype.patch`.
Solver/build files remain exactly Stage266 `ff740744`. Do not promote or package
that prototype as a qualified candidate.

The next optimization should follow the largest *avoidable added work*, rather
than the next easy Domain reader. Removing Domain dependencies remains the end
goal, but a numerical migration can add transport, copies or duplicate derived
calculations. Measure the complete replacement cost before selecting it.

## Evidence for changing priorities

The [Stage266 result](parallel-blocks-stage266-results.json), job 55398, establishes
the current 83-rank comparison: candidate 0.748448779 s per ordinary step versus
same-allocation main 0.541753610 s, a 206.695 ms gap. It contains no isolated
surface-pressure measurement.

The earlier [Stage229 audit](parallel-blocks-stage229-review.json), job 54737,
located the inherited candidate/main excess principally in wavelet transactions
and adaptation/reconstruction:

| Historical region | Signed candidate minus legacy time per step |
|---|---:|
| Dynamics wavelets, including transaction preparation | +136.817 ms |
| Adaptation and reconstruction | +77.664 ms |
| Remaining dynamics | +19.765 ms |
| Tendency context | -27.882 ms |

These are investigation leads, not current savings estimates. The legacy
observer failed its overhead calibration; its whole-step shift was +11.223 ms.
The observers used different historical sources and a different legacy anchor.
Region means are not ordinary-step run medians, and MPI waiting is included.
Do not subtract Stage265t/266 gains from these numbers, sum positive regions as
a recoverable budget, or interpret the tendency offset as proof that every
individual tendency operation is efficient. Still, the evidence contradicts
prioritizing a shared pressure routine solely because its migration is obvious.

The later Stage251--262 owner-velocity regressions belong to the abandoned
post-Stage230 design. Their large owner-lookup timings must not be attributed to
the reset Stage266 source. They illustrate why a measured source transition is
more useful than carrying a hotspot label between implementations.

## Ranked investigation targets

This ranking is a priority for gathering evidence, not a measured ranking of
current net savings. Compare all major transaction families in one bounded cost
census before choosing one intervention.

1. **Wavelet/inverse transaction overhead.** Inspect `WT_after_step`
   (`src/adapt.f90:380`), provisional and accepted entry points, native inverse
   preparation, seed/scaffold staging, boundary epochs, canonical publication,
   final installation and ghost refresh. Compare with main's complete wavelet
   transaction, including its necessary reconstruction. Stage265t and Stage266
   already recovered work in this family; that supports investigating the
   remaining repeated work, without establishing a remaining budget. The
   Stage267 audit still prohibits blindly skipping complete publication.

2. **Adaptation and complete block reconstruction.** `src/main.f90` runs Domain
   physics/adaptation and inverse reconstruction, then rebuilds the catalogue,
   source blocks and migrated compact state, and calls
   `retain_post_grid_change_block_reconstruction` (`parallel_block_mpi.f90:3330`).
   Price topology construction, source extraction, packing/local installation,
   transport and boundary/ghost re-establishment separately. Audit unchanged
   topology and unchanged-block opportunities. Retain actual post-physics data
   refresh and genuine remapping; skipping all reconstruction would leave stale
   payloads. Check active selectors before counting serialization proofs as a
   production cost: the fast-install option can already suppress them.

3. **Eager hydrostatic products without numerical consumers.** Every provisional
   refresh explicitly invalidates and rebuilds the full block hydrostatic cache
   (`parallel_block_mpi.f90:11498`), including surface pressure, Exner and air
   temperature. `compute_local_block_hydrostatic_patch`
   (`parallel_block.f90:5220`) evaluates a power for each physical layer. The
   cache's private numerical arrays are exposed through getter/consumer APIs;
   current external numerical traversals are the hydrostatic oracle and the
   reconstruction validation suite. Production tendency checks cache readiness
   in `evaluate_candidate_block_exner_shadow` (`parallel_block_mpi.f90:28627`),
   but its kernel computes separate `native_columns`. Readiness requirements
   alone do not demonstrate consumption of the first cache's values. Audit
   lifecycle checks and all configurations before omitting work.

4. **Repeated thermodynamic/velocity recomposition.** Shared Domain primitives
   form pressure/Exner, `compute_native_velocity_gradient` forms the native
   velocity tendency, and `compute_block_exner_difference_kernel`
   (`parallel_block_mpi.f90:34224`) reconstructs native columns and gradient
   terms. The retained raw trend then participates in an add-then-subtract
   pressure-term round trip. This is a promising duplicate-producer chain, but
   exact rounding is part of the accepted result: replacing `(trend+D)-D` with
   `trend` or substituting another pressure cache is not an exact simplification.
   Price columns, addressing, arithmetic and transport; a later replacement must
   preserve the required numerical policy explicitly.

5. **Provisional surface-pressure reader.** Keep the Stage267 migration design
   available, but lower its priority until its cost is isolated. The existing
   pressure calculation is also present in main, and changing just this reader
   leaves full Domain publication and other shared readers intact. Its Domain
   dependency is removable independently; its performance benefit is unproven.

## First bounded change to investigate

The strongest new source-backed omission candidate is **lazy provisional block
hydrostatic cache construction**. Preserve ghost refresh and cache invalidation;
build the cache when an actual numerical consumer requests it, including every
oracle. Replace readiness-only production demands with the relevant input/epoch
contract instead of computing unused pressure, Exner and temperature to satisfy
the ready flag. Initially restrict the intervention to provisional RK exits.
Leave accepted/adaptation/remap cache handling for a later independently tested
change.

This is not yet an approved omission: prove all private-array/API consumers,
positivity/coverage contracts and stage epochs first. Count actual cache refreshes
and columns on every rank to show they disappear from the production path and
remain available to the oracle. Verify the independent native-column path does
not read these arrays. The final kernel currently uses different pressure
rounding, so do not make it reuse the eager cache as part of this change.

If the census shows this cache is a small cost, select a larger recoverable
wavelet or reconstruction operation instead. Do not launch a pressure migration
campaign merely because a prototype exists. Prefer interventions with a credible
route to tens of milliseconds; a smaller migration can still qualify if it is
needed to unlock a larger consumer/adapter removal and adds no repeatable penalty.

## Measurement and selection protocol

Keep exact Stage266 and read-only legacy main as controls. Add bounded
transaction scopes and work counters to diagnostic copies, with no arithmetic
changes, layer/node clocks, extra synchronization or timed-step output. Record
per-rank/per-step work, self/inclusive wall time and process CPU; flush after the
run. Distinguish production from oracle costs and unchanged from changed topology.
Report actual calls, columns, copied bytes, messages and plan rebuilds.

Cover wavelets/inverse, adaptation/reconstruction, hydrostatic cache construction
and tendency thermodynamic recomposition in the same diagnostic run. Calibrate
the observer against bare controls, using the existing ordinary-step window
(exclude steps 1, 8 and 9). Do not treat process CPU, average rank time or sums of
independent region maxima as the elapsed critical path. Preserve signed offsets
and attribute waiting to its upstream work instead of counting it twice.
Historical failed observer calibration must be visible in any reused report.

Rank an intervention by demonstrated redundant work, its critical-path location,
and the replacement's transport/storage cost. Profiling supplies candidates;
an unprofiled single-change comparison establishes the recoverable gain.
For the selected implementation, retain exact optimized/checked/oracle local
replays with serial builds, at most four MPI ranks and one library thread.
Then use the established short nk/bb13 campaign: seven correctness runs at
16/167 ranks and six 83-rank timing runs ordered
`main, Stage266, candidate, candidate, Stage266, main`. Preserve source/fixture
seals, shared 23 physics objects, placement, checkpoint/VTK/remap and zero-Domain
rank gates. Repeat only for small/noisy/mixed evidence; do not stack removals.

## Possible speedup

There is no defensible numerical forecast for pressure migration or lazy cache
construction before their current isolated costs are measured. A region's total
cost is only an upper bound on an omission, and local four-rank timing cannot
predict bb13's critical path. The following are conditional arithmetic scenarios,
not expected gains, using the measured Stage266 baseline:

| Net saving per step | Candidate step time | Elapsed-time reduction | Remaining gap to the measured main |
|---:|---:|---:|---:|
| 5 ms | 0.743 s | 0.67% | 201.7 ms |
| 25 ms | 0.723 s | 3.34% | 181.7 ms |
| 50 ms | 0.698 s | 6.68% | 156.7 ms |
| 100 ms | 0.648 s | 13.36% | 106.7 ms |

Matching the measured main requires about a 27.62% reduction in candidate elapsed
time. This is why the next audit should target major duplicate transactions;
there is no evidence that the surface-pressure reader alone can close that gap.
The evidence seals, current source anchors and scenario calculations are in
`parallel-blocks-stage268-priorities.json`. No new performance kit is sealed by
this prioritization audit.
