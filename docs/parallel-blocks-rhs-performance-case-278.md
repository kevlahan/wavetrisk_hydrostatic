# Connected RHS: quantitative performance case and revised implementation

Current completion status: the connected production selector and final local
qualification are complete; the sealed bbserv kit is ready for user submission.
See [the production result and run commands](parallel-blocks-rhs-production-278-results.md)
for current evidence and remaining goal gaps. The design/reassessment text below
is the chronological pre-implementation record; its draft-readiness statements
are superseded by that result. Accepted source remains266; no cluster promotion.


Prepared 2026-10-09. Accepted source remains Stage266 `ff740744`.
**Decision: revise the implementation before a cluster launch.** Existing evidence
supports targeting redundant capture and tendency transport while moving the
connected RHS to block owners. It does not justify adopting the inactive draft's
whole-array layout or forecasting closure of the legacy gap. This document adds
a quantitative gate to the [implementation/campaign plan](parallel-blocks-next-connected-changes.md).
No production change, accepted speedup or completed ownership migration is claimed.

## Reproduced evidence, with compatible denominators

The accompanying `performance_case.py` verifies **998 raw observer files** against
the archived Stage270/272 manifests, reproduces relevant Stage270 component costs,
and checks that Stage272 integer work inventories repeat exactly between runs
and conserve the published global inventory. Outputs and input hashes are saved
in `.stage278-validation/performance-case.json`. It uses Stage272 counts only:
its failed timing calibrations are preserved. No solver was rerun.

[Stage270](/Users/kevlahan/.codex/worktrees/26d8/wavetrisk_hydrostatic/docs/parallel-blocks-stage270-results.md)
has compatible candidate-only timing calibration at83 ranks. The two-run bare
accepted mean is749.653ms, main540.342ms, gap209.311ms. Cross-version phase
attribution failed calibration. The following are mean-rank costs averaged over
ordinary steps and two observer runs, not synchronized critical-path savings.

| Region | ms/step | Treatment in the performance case |
|---|---:|---|
| Shared compatibility producer, self |193.987| Necessary arithmetic plus existing waits/overheads; replaced, not deleted. No assumed compute speedup at83. |
| Scalar capture preparation, inclusive |20.130| Mixed plan/storage/coverage/producer-buffer work. Only the eliminated fraction earns credit. |
| Scalar capture publication, inclusive |23.776| Mixed mass-RHS substitution, transport, install and boundary finalization. Needed boundary work is retained or charged to its replacement. |
| Velocity remainder preparation, inclusive |5.289| Completed native tendency export and installation; retiring it requires direct resident tendency production. |
| Scalar restriction transaction, inclusive |87.574| Existing useful block-owned work; retain its numerical kernels and required halos. No blanket removal credit. |
| Resident RK state assembly, inclusive |1.178| Already resident; reuse. This is not a measurement of removable Domain RK closure. |
| Canonical inverse publication |1.092| Retained reader contract; no removal credit. |
| Native scalar forward, inclusive |55.486| Outside this replacement; no removal credit. |
| Migration, inclusive |64.341| Outside this numerical slice; no removal credit. |

Capture preparation, capture publication and velocity remainder are sequential
sibling calls at the tendency entry. Their combined **49.196ms** is a candidate
cost pool, not49.196ms of waste. In particular, capture preparation allocates/reuses
scalar storage and restriction/edge plans; publication also finalizes boundary
records. The2.244ms capture MPI cost is already inside publication. The Domain-to-
block payload timer overlaps several parents, including velocity remainder; it
must not be added. Domain RK compatibility is not isolated by these coarse timers
and receives zero additional credit. The193.987ms producer self cost excludes
nested measured children; unmeasured producer costs cannot be separated further.

Even the unrealistic case of deleting the entire49.196ms pool at zero replacement
cost corresponds to only6.56% of this allocation's accepted time, or23.50% of its
legacy gap. This is an accounting ceiling for that pool under unchanged other
work, not a bound on all possible architectural benefits. It demonstrates why
this slice alone cannot be sold as legacy parity from current evidence.

## Conditional break-even calculation

At83 define:

- R: fraction of that49.196ms pool actually eliminated, from0 to1.
- d: fractional change in the193.987ms producer self cost (positive is slower).
- N: incremental preparation, halo, adapter and memory-related execution costs
  not already included in d or the retained part of the old pool, in ms/step.

The screening equation is `saving = 49.196 R - 193.987 d - N`.
The new path must satisfy `N < 49.196 R - 193.987 d` for positive savings in this
accounting model. Allocate each cost once. Existing retained scalar restriction,
geometry metadata and adapter work must not be credited as removed and then
omitted from replacement costs. If they change, their net change belongs in N.

| Assumed pool removal | Producer change | Added N | Conditional saving | Fraction of same-allocation accepted time |
|---:|---:|---:|---:|---:|
|25%|unchanged|10ms|2.30ms|0.31%|
|50%|unchanged|10ms|14.60ms|1.95%|
|50%|5% slower|10ms|4.90ms|0.65%|
|50%|10% slower|10ms|−4.80ms|−0.64%|
|50%|5% slower|20ms|−5.10ms|−0.68%|
|75%|unchanged|10ms|26.90ms|3.59%|

These are explicit sensitivity assumptions, **not measured replacement costs,
probabilities, confidence intervals or a predicted gain range**. No evidence
currently establishes50% removal or10ms added cost. For that50% removal with a5%
producer penalty, the break-even N is14.90ms/step (3.72ms per RK RHS on average).
That is a derived design budget to test, not an agreed performance threshold.
The4.90ms outcome is comparable to Stage270's4.75ms spread between bare run
medians; that spread is not a confidence interval, but warns against a strong
claim from one small difference. Bare whole-step measurements remain decisive.

Do not scale the83-rank components by83/160 or83/256: valid phase costs at those
counts do not exist. In the latest [Stage277 controls](parallel-blocks-stage277-results.md),
accepted160 and256 times differ by39.388ms. If the candidate saves delta160 and
delta256 respectively, actually becoming faster at256 than160 requires
`delta256 - delta160 > 39.388ms` for that baseline observation. Merely saving time
at256 is insufficient if160 improves equally. Re-evaluate with matched controls.

## Quantitative communication boundary

Recompute the [Stage275](/Users/kevlahan/.codex/worktrees/26d8/wavetrisk_hydrostatic/docs/parallel-blocks-stage275-results.md)
field budgets for two alternatives. Start with its hypothetical upstream-owned
budget; remove rho/rho_theta transfers assuming values and halos are resident at
the real consumption epoch. Then either compute RHS physics locally or add back
the measured protocol budget for its Domain-produced residuals.

| Ranks | Upstream + resident density + local RHS physics / old velocity export | Same but retain Domain RHS physics | Remaining partial byte headroom with local physics |
|---:|---:|---:|---:|
|83|0.542–0.554|1.105–1.151|608,208–636,784 bytes|
|160|0.568–0.592|1.091–1.159|417,152–443,952 bytes|
|256|0.648–0.672|1.170–1.239|338,768–364,720 bytes|

Each range covers the existing three rank-local structural sample groups, per
layer. They are neither synchronized global epochs nor whole-step totals.
Do not multiply them by30 layers and48 calls to invent trajectory traffic.
They exclude new producer-input/derived-field halos, scalar-pipeline traffic,
message count/latency, packing and computation. Resident-density subtraction is
conditional on actual ghost/alias readiness, not just owner-interior agreement.

Thus a retained RHS physics adapter cannot be assumed cheap: this partial input
budget already exceeds the old velocity export, even before producer halos.
At256, the local-physics alternative leaves only about0.34–0.36MB of same-sample
byte headroom. Any extra required halos consume it. Equal bytes also need not mean
equal time; peer fanout and phase count must be compiled and counted explicitly.
No network bandwidth or latency is inferred from MPI waits.

This concerns the RHS callback's diffusion/physics residuals, not the separate
post-RK Simple split-physics transaction. The latter remains an explicitly
budgeted adapter until its own consumers migrate. It earns no savings here.

## What existing ownership can and cannot buy

New reduction of the raw [Stage272](/Users/kevlahan/.codex/worktrees/26d8/wavetrisk_hydrostatic/docs/parallel-blocks-stage272-results.md)
integer census gives the following maxima over ranks. Each count sums patch
inventories over the four RHS stages of a step. It is a work proxy, not a measured
critical path, and it excludes differences in patch/halo/physics cost.

| Ordinary step | Max Domain patch-stage count | Max block count at160 | Max block count at256 |
|---:|---:|---:|---:|
|2|204|88|88|
|3|204|88|88|
|4|212|88|88|
|5|212|88|88|
|6|212|88|88|
|7|236|96|88|
|10|104|64|44|
|11|128|72|68|
|12|224|92|88|

Both runs reproduce these counts exactly. Under the explicit equal-cost-per-patch
hypothesis, moving Domain work to the current block distribution lowers this
maximum proxy to0.407–0.615 of its former value at160 and0.373–0.531 at256.
That supports a load-distribution hypothesis at and above D, without predicting
seconds saved. Individual stage maxima can differ from these stage-summed maxima.

But the block maximum is identical at160 and256 for five of nine ordinary steps.
The mean of the listed maxima falls only4.71% (ratio0.95288). It does not support
assuming a37.5% reduction in producer duration from160 to256. Uniform patch cost,
sequential stage waits and equal per-rank throughput are all unproved. Keep
partition/subdivision separate; the unchanged catalogue is the control for this
numerical slice. Efficient640-rank execution still requires finer granularity.

## Revised implementation contract

The logical production boundary stays connected RHS through real resident RK.
The revised layout and retirement requirements are narrower and more explicit:

1. **Write resident results directly.** Keep accepted block temperature
   restriction/divergence kernels and resident RK registers. Replace the
   Domain capture/export/install chain with producer writes into the actual
   resident fields those consumers use. Do not retain full capture streams
   beside a second RHS workspace. Keep only required boundary/scaffold records.
   This is the source-level claim that could earn part of the49.196ms pool.
2. **Put RHS diffusion/physics with its producer.** Port its required ordered
   hierarchy and halos with direct kernels; do not fill physics_flux and
   sw_physics from recurring Domain residual exports. Retain the separate split
   Simple physics adapter explicitly. The traffic table supplies the reason.
3. **Use patch-contiguous production loops and shared geometry.** Reuse existing
   compact scalar storage and epoch-valid metrics where compatible. Compile
   boundary/exception addresses once per valid topology/layout/mask epoch;
   ordinary interior loops use direct indexed patch arrays. Avoid duplicating
   the draft's49-double per-event geometry and traversing generic slot resolvers
   or shape/address validation at every numerical operation. Keep checked bounds,
   plan-install validation and runtime epoch checks. No assumption of a measured
   cache bottleneck is required to avoid adding these costs.
4. **Bound temporary storage by phase.** Stream pressure integration and primitive
   work in the accepted vertical/level order. Preserve the raw versus restricted
   Exner distinction and final-writer order. Keep required persistent Exner/state
   fields; do not allocate another complete node-by-layer copy of every output
   just because the numerical helper accepts that shape. Supply a storage ledger
   with owned/ghost capacity and geometry replication before activation.
5. **Compile and count the real replacement route.** Enumerate producer-input,
   derived-field, restriction and boundary dependencies from the actual plan.
   Record bytes, messages, peers and synchronization phases by state/mask/topology
   epoch, including mask-only rebuilds. Charge producer-buffer zeroing, packing,
   plan rebuilds and retained publication once. Keep Domain RK duplication at
   zero savings credit until its actual removed reader paths are demonstrated.

These are implementation requirements, not a new standalone optimization or
permission to skip any dependency. If local residency requires more routes than
the partial headroom allows, revise that field's storage/consumer boundary before
activation. No fine-grained interpreter, bulk primitive return, bundled Stage273
partition change or separate Stage277 cache experiment is justified by this case.

## Decision gate before the proposed campaign

The immediate implementation task is the real patch/route compiler and direct
resident producer-to-consumer layout above, developed as part of the complete
selector. Its review must populate the storage/traffic ledger and identify actual
removed call paths. A narrowly scoped local timing comparison may reject a bad
loop/layout using identical arithmetic and the actual compiled plan. Local timing
cannot certify bbserv milliseconds or fill missing160/256 timing attribution.

Before sealing the full candidate, require a traceable estimate for R, the
producer throughput change and incremental costs, with unresolved values explicit.
A case with no positive margin under plausible costs should be revised rather
than sent as a presumed optimization. The named unknowns are (a) retained fraction
of capture setup/publication, (b) complete route fanout/halo cost, and (c) actual
fused-loop throughput and load distribution. Current evidence answers none with
a calibrated whole-step bound; a guaranteed speedup is not claimed.

The already proposed bounded bb13 campaign qualifies a completed candidate and
tests this hypothesis; it is not authorization for a new generic profiling run.
Keep exact science, zero-Domain authoritative execution, retired work and bare
full-step performance as joint gates. Reject repeatable regressions, report the
remaining main gap, and leave the640-rank goal open. **Launch readiness remains
false.** The revised implementation has a concrete cost target and falsifiable
failure conditions; it has not yet earned a positive performance forecast.
