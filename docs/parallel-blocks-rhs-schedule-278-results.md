# Resident producer schedule — Stage278 development

Prepared 2026-10-09. Accepted solver remains Stage266
`ff740744608731a05d53771c32a3055309f3a586`. No candidate is promoted and no
bbserv job or launchable kit exists. This is another implementation prerequisite;
it does **not** finish the requested production integration or campaign preparation.

## Change and reason

`src/parallel_block_rhs_schedule.f90` compiles the actual shared-producer visit
sequence over resident blocks: interior direct visits, immediately adjacent SW
writes, positive north/east boundary visits and the northeast corner. A caller
supplies patch order and pentagon identities explicitly. The compiler preserves
that order rather than inferring numerical order from block DFS storage. It
provides flat owned/boundary/ghost offsets and separate local final-writer tables
for mass flux, PV and vorticity, including all four mapped pentagon corrections.
These offsets support direct derived-field addressing; the resident state bindings
and numerical kernels still need integration. No new full-state packing is
qualified or assumed free.

`parallel_block_rhs_layout.f90` now exposes topology-time flat-slot mapping,
rejects missing boundary storage explicitly and includes physical-boundary side
classification in its exact generation snapshot. A side-classification change
cannot reuse an old layout generation. Numerical masks remain a separate,
unimplemented contract for mass/velocity and diffusion plans; shared step1
traversal itself is not selected by the node/edge numerical masks.

An initial real-layout replay failed at the first positive-east visit. The
accepted `comp_offs3` uses the SOUTH strip height to construct the SOUTHEAST
offset. For a 4x3 SOUTH strip adjoining a 3x4 SOUTHEAST strip, a generic
`idx(PATCH_SIZE,-1)` picks a different row. The new schedule now reproduces the
accepted address, with a mixed-dimension regression case. This preserves science;
it does not change or claim to repair the legacy operator. Failed source,
binaries and raw logs remain in `.stage278-validation/schedule-failed`. A strict
build also rejected a guarded `corner(s-1)` access; the comparison now uses a
separate loop starting at two. Pre-repair evidence is retained under
`schedule-replay/before-warning-repair`, and all three replays were repeated
against the corrected source.

## Evidence

Two serial full climate builds, optimized and strict checked, plus three
sequential four-rank J4 replays (optimized, independent oracle, checked) pass.
Every replay retains 12 exact step states, two reloads, two remaps, exact final
checkpoint fields/topology and 62 exact VTK files against accepted266.

For **each replay**, all 48 RHS calls and 12 observed topology generations check:

| Schedule preparation observation | Count |
|---|---:|
| Direct producer events, including positive boundaries | 2,029,280 |
| Immediate SW events | 171,532 |
| Positive-boundary direct events | 187,808 |
| Supplied patch orders tested | forward and reverse |
| Address / visit / local-writer differences | 0 |

The event counts above count one order; both orders are checked. Schedule checks
are topology checks once per RHS, **not** all-layer primitive arithmetic checks.
The existing sampled Stage276 input/arithmetic checks also pass unchanged.
The independent schedule oracle extracts the accepted step1 traversal and
`comp_offs3` offset transformations directly from the accepted Git source.

Optimized and checked standalone units pass all 16 physical-boundary combinations,
mixed strip dimensions, owned/internal/ghost addressing, supplied patch ordering,
all four pentagon writer mappings, and eight expected rejection cases per build:
stale generation, duplicate/missing patch, wrong level, wrong pentagon order,
invalid corner, changed boundary side without generation and missing storage.

Evidence: `.stage278-validation/schedule-validation.json`, `schedule-replay/`
and `schedule-units/`. Original paths in raw replay reports record the temporary
execution location; their files are now retained under `schedule-replay/`.

## Replacement cost and dual-goal score

The schedule allocates44 bytes per direct event (seven offsets, patch identity,
coordinates and SW flag), excluding compiler descriptors, pentagon entries,
existing input-layout snapshots, geometry, state/derived fields and routes.
Observed maximum summed schedule payload per local rank is
471,284–536,844 bytes on this J4 fixture.
This is a cost to include in the connected candidate, **not a saving or a J5
memory forecast**. Dense final-writer tables are transient compilation scratch;
this validation deliberately recompiles plans to check them on every RHS.
Production epoch caching and its cost remain to be integrated and measured.

| Scorecard gate | Result |
|---|---|
| Actual resident boundary traversal prerequisite | Locally checked, including the previously missed mixed-dimension alias |
| Authoritative block production / real RK readers | Not connected |
| Retired Domain computation or transfers | None |
| Complete connected-candidate science | Not qualified; accepted production remains active in this overlay |
| Legacy gap / P>D performance | No new evidence or defensible speedup forecast |
| Useful J5>=640 execution | Unqualified |
| Seal / launch | Not ready |

## Remaining work before the planned bbserv campaign

1. Distribute accepted level order, physical pentagon identities, numerical masks
   and geometry with topology/ownership epochs; replace local writer provenance
   with actual cross-block final-writer delivery. Root/scaffold producers and
   halo state readiness must be covered explicitly.
2. Complete pressure/Bernoulli/Exner restriction and scalar/divu/rotational
   diffusion routes, including `post_vort`, using the actual resident workspace.
   Connect their outputs to the existing native mass/source/gradient kernels and
   the temperature pipeline with correct per-level and per-layer epochs.
3. Cut over `prepare_block_multistage_tendency` and actual resident RK consumers;
   retire the displaced Domain producers, mass RHS publication, velocity export
   and duplicate interior RK closure in the same independently selectable slice.
   Keep inverse/boundary and separate Simple split physics adapters explicit.
4. Qualify that active slice locally with optimized/checked/independent-oracle
   science, retirement assertions, actual route/storage counts and replacement
   costs. Then adapt and seal the proven launch/status/collection kit and verify
   a fresh extraction against the immutable sources and fixture identities.

The existing16-run bb13 specification remains unchanged: four serial builds,
eight validations and eight bare main/266/candidate comparisons at83/160/256.
Both scorecard goals remain mandatory: independent ownership with removed work,
and complete-step performance toward legacy parity and useful>=4D execution.
The49.196ms mixed pool at83 remains only a sensitivity budget; this schedule has
removed none of it. The later bb inter-node check and separately qualified finer
partition/160–320–640 campaign remain distinct from this numerical cutover.
