# Stage269 result: useful candidate costs, failed main calibration

Job 55466 completed all fifteen runs. Independent review verifies all 2,540
archive files, matches the delivered kit seal, and reproduces the scientific,
placement, storage, checkpoint/VTK, raw observer and timing gates. Five serial
same-node builds shared 23 physics objects. Checked execution and both 16/167-rank
oracles passed; each oracle retained 48 exact returns and 84 installations.
Production source remains Stage266 `ff740744`; no optimization was tested.

The bare 83-rank controls give mean ordinary-step run medians of **0.538908767 s
for main** and **0.751216972 s for Stage266**. The gap is **212.308 ms (+39.396%)**.
This is a new allocation's comparison of unchanged source, not evidence of a
new implementation regression relative to the earlier Stage266 campaign.

## Calibration and usable evidence

| Observer | Mean added time | Relative shift | Prespecified decision |
|---|---:|---:|---|
| Main | 21.147 ms | 3.924% | Failed |
| Stage266 | 6.929 ms | 0.922% | Passed |

Main exceeds both mean and pair limits; its two shifts are 20.498 and 21.796 ms.
Stage266 passes all limits, though its mean shift is close to the 7.512 ms limit.
Cross-version phase differences remain **unqualified**. Do not subtract a flat
21 ms correction from main or use the signed phase ranking as a qualified
breakdown of the 212 ms run-median gap. Candidate-only broad costs are useful
investigation leads at this observer resolution; small costs are not precise
speedup budgets.

| Candidate scope | Mean self time per ordinary step | Interpretation |
|---|---:|---|
| Tendency context | 163.730 ms | Broad unsplit tendency work; children excluded |
| Native recomposition | 108.869 ms | Scalar restriction/transport plus final physical kernel/retention |
| Migration | 64.353 ms | Packing, exchange, installation and state/plan preparation |
| Inverse MPI wait | 52.782 ms | Waiting; upstream imbalance may contribute |
| Native scalar forward | 51.599 ms | Native scalar transform family |
| Inverse setup | 26.700 ms | Setup self cost; nested inverse work excluded |
| Hydrostatic cache | 2.116 ms | All measured actual cache rebuilds |
| Surface pressure | 1.549 ms | Shared pressure routine also exists in main |

Self times conserve each rank's local step; these are mean-rank costs, not
elapsed critical-path savings. Broad inclusive scopes overlap. Candidate
recomposition self times are 109.605 and 108.134 ms in the two runs, and migration
65.467 and 63.239 ms. Both remain substantial in the composition of each step's
longest local-duration rank. That composition selects one rank per step, not
independent maxima, and does not prove the global critical path.

The small cache/pressure costs lower the priority of lazy cache construction and
surface-pressure migration for recovering this gap. Domain migration remains
required eventually, but these measurements supply no basis for expecting those
omissions alone to recover tens of milliseconds.

## Source-backed next investigation

`evaluate_candidate_block_velocity_recomposition` first runs
`evaluate_candidate_block_scalar_restriction`, including direct scalar flux,
per-level ghost/temperature edge exchanges, coarse restriction and divergence
replay. It then runs the Exner/complete physical kernel and retains the tendency
in the accumulator. Therefore the 109 ms label does **not** isolate column or
Exner arithmetic. Price scalar restriction kernels and transport separately
from final recomposition, metric preparation and accumulator retention. Also
split the larger 164 ms tendency self region into its actual producer/transport
operations before assuming it is irreducible legacy work.

`migrate_blocks` includes manifest construction, size/payload exchange, packing,
unpacking, installation, ghost/writeback plan reconstruction and production
workspace preparation. It currently occurs once per ordinary step. The 64 ms
label does **not** isolate network transfer. Price these components, then audit
which geometry/plans can survive unchanged topology while still refreshing
post-physics prognostic values. Do not skip the whole migration/reconstruction
transaction or stale the numerical payload.

A smaller source lead is `difference_moment`: the final Exner kernel accumulates
three moments for every vector sample, but the source search finds no later
consumer of that member. That diagnostic arithmetic may be removable; its cost
is not isolated, so it should not displace the major producer/transport audit.
Keep exact add-then-subtract velocity ordering and independent numerical oracles.

For the next focused audit, retain whole-transaction scopes but remove clocks
from Domain boundary routines and per-patch velocity routines. Main currently
executes about 3,150 observed boundary calls per ordinary step; the candidate
about 2,419. Their clock/CPU-time overhead is a plausible cause of main's failed
calibration, not a measured attribution of that failure. Add a small set of
coarse scopes for restriction/local kernels/exchanges/final kernel/retention and
migration packing/transport/installation/plans. No node/layer clocks or new
synchronization. Recalibrate both bare controls; repeat neither this unchanged
high-frequency observer nor its failed comparison.

After that decomposition, select one demonstrably redundant operation and test
its complete replacement unprofiled against Stage266 and main. Preserve the
Stage267 reader/epoch requirements. No source removal is promoted by Stage269,
and no numerical speedup forecast is justified yet.

## Reviewer compatibility repair

The first local review rejected four derived longest-local-step sums differing
by 1.1e-16 to 1.8e-15 seconds. Raw records and hashes were identical; built-in
float summation differs between the cluster and local Python versions. The
repository reviewer now permits up to eight ulps only when comparing derived
float reports. Hashes, integer counts, keys/order, raw record conservation and
all numerical acceptance gates remain exact. Its entry point also resolves the
repository's shared helpers without requiring a custom PYTHONPATH. The original
delivered kit and archive remain immutable.

Ten workflow tests passed, including a new rounding test that accepts a one-ulp
report difference and rejects changed counts/hashes, NaNs, extra keys and material
time changes. The complete real archive then verified. The [machine-readable
review](parallel-blocks-stage269-results.json) preserves both calibrations,
unqualified cross-version attribution, candidate-only costs and limitations.
Review archives omit binary/checkpoint/VTK payloads; their cluster-computed
hashes and numerical reports are verified, without rerunning those payloads.
