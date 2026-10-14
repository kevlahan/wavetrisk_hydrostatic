# Stage278: assignment, granularity and strong-scaling assessment

2026-10-14. Development source: `de558ed03798f721991f38eeecc00f00eafc8370`
on `codex/stage278-structural-repair`. Accepted solver remains Stage266; legacy
main is unchanged. This is an offline assessment of verified job55504, not a
new solver change, campaign or promotion.

## Decision

The evidence supports **both a work-distribution problem and growing overhead**.
It does not support a standalone claim that better leaf-weight partitioning,
communication batching, or simply using more ranks will recover scaling.
At256, heavy single-block work bundles already approach the observed maximum
work. Better assignment is worth evaluating, but finer ownership units and the
boundary/scaffold work attached to them must be addressed for substantial gains.
Separately, retained Domain dependencies and preparation/exchange costs need
removal or reduction; distributing existing RHS arithmetic alone cannot plausibly
close the main gap. These are complementary requirements of the original goals.

## Evidence scope and reproducibility

This extends the earlier trajectory-level follow-up by selecting exactly the
nine ordinary steps (one-based2–7 and10–12). Remaps remain included. It reads the
existing validation rank counters and the separate coarse-profile replays at
83/160/256; no new cluster or local solver execution was performed.

[Offline analyzer](../test/parallel_block_profile/structural_278/assess_balance.py)
subtracts cumulative counters only within the same topology/mask epoch, groups
four RK calls per step, and verifies the result against the independently reviewed
trajectory totals. Every ordinary step has identical global primitive, scalar
and velocity counts at all three rank counts. The nine steps contain53,328,240
primitive events at each count. It checks complete rank/step profiles, self-time
partitioning, and agreement with the existing trusted attribution.

[Machine-readable result](../test/parallel_block_profile/structural_278/qualification/balance-assessment-55504.json)
contains each step's distributions, phase summaries and SHA256 hashes of its
input files. Reproduce from the preserved return:

```bash
python3 test/parallel_block_profile/structural_278/assess_balance.py .stage278-validation/review-55504/review/results /private/tmp/balance-assessment-55504.json
```

Counters and timings are from separate exact-science replays, joined by rank and
step. Their association is descriptive; it is not a synchronized causal trace.
Primitive events count diffusion-producer events times layers, not total FLOPs
or complete timestep work. Reported correlations are descriptive, with no claim
that ranks or successive steps are independent statistical samples.

## 1. Assignment quality: real defects, incomplete cost model

Current weights are descendant leaf-patch counts times `PATCH_SIZE**2`
(`src/domain.f90`, `subtree_weight_Domain`). They omit internal-patch work,
boundary/scaffold work, per-block/per-rank setup and retained Domain costs.
`compile_rhs_owner_map` assigns unresolved boundary/scaffold numerical targets
to the owner of the Domain's first catalog block. That extra burden is not an
explicit term in the leaf weight.

The catalog assignment in `src/parallel_block_mpi.f90` uses increasing-rank
next-fit against a fixed global average, relaxing tolerance until all blocks fit.
It does not update the target for remaining ranks/work or minimize the heaviest
rank's actual numerical/retained work. These are different defects: an optimal
partition of inaccurate weights can still be a poor runtime assignment.

| Ordinary-step work measure | 160 ranks | 256 ranks |
|---|---:|---:|
| Mean primitive events/rank/step | 37,033.5 | 23,145.94 |
| Mean of each step's maximum primitive count | 50,240 | 46,186.67 |
| Mean of each step's maximum/mean | 1.354 | 1.989 |
| Range of per-step maximum/mean | 1.255–1.388 | 1.525–2.120 |
| Maximum/mean after averaging work over steps for each rank | 1.265 | 1.734 |
| Correlation of rank ordinary-mean work with rank arithmetic time | 0.929 | 0.966 |

Averaging first hides changing critical ranks. The per-step maximum falls only
**8.07%**, while average work per rank falls37.5%. The previously reported15.9%
maximum-work fall used entire-trajectory rank totals; the ordinary-step statistic
above is the appropriate stronger evidence here.

At256 the160 Domain owners receive about77.22% of ordinary primitive work,
while96 zero-Domain ranks receive22.78% (means of step shares). Their mean
per-rank counts are28,684 versus13,916. Domain owners already have more RHS work
before considering retained Domain work. This contradicts an assumption of equal
work simply because every rank participates. It does not prove how much of the
full-step delay each owner causes.

The coarse logs also report leaf-load maximum/mean of1.068 then1.200 at160,
and1.139 then1.920 at256. These are two catalog snapshots, not per-step work
measurements. Their discrepancy with numerical-work statistics cannot be used
as an aligned quantitative proxy-error estimate. The full J5 per-block weight
and operation vectors are absent from this return.

Historical Stage272 independently found front-loading. For640 equal blocks on
256 ranks, old next-fit gives192 ranks three blocks and64 one block; a feasible
balanced assignment gives128 ranks three and128 two. **Both have the same peak
of three**: removing an underfilled tail need not shorten the critical path.
Stage273 already implemented a locally tested contiguous minimax partition;
reuse its reference/evidence when evaluating assignments, rather than presenting
it as a new idea. The reviewed Stage273 preparation report calls it unpromoted
pending cluster performance; this assessment has no cluster rejection result
for273. Stage277's small mixed timing result is a separate cache experiment.
Neither candidate is imported into v7 automatically.

## 2. Granularity: quantified limitation of assignment alone

All these J5 catalogs contain640 blocks, at most four root-child subtrees per
Domain. That gives four blocks/rank at160,2.5 at256 and at most one on average at640.
An adaptive set of unequal blocks cannot be flexibly balanced with one per rank.
More ranks require subdivision, work sharing, or another finer numerical unit;
changing the partition algorithm alone cannot create that flexibility.

The ordinary-step counters sharpen this conclusion. At256:

- On steps7,10 and11 the heaviest rank owns exactly one block, carrying51,600,
  30,360 and42,840 primitive events respectively.
- Across the nine ordinary steps, the largest single-block-rank work averages
  43,653.33 events; the actual rank maximum averages46,186.67.
- In every ordinary step, that single-block work exceeds global work/256:
  its ratio to the ideal mean ranges1.525–2.120.

For a model that preserves each existing block's numerical-work bundle,
including its currently attached boundary/scaffold targets, the indivisible-work
lower bound is `max(total_work/P, largest_block_work)`. A singleton rank exposes
one such bundle, so it provides a lower bound even without the entire catalog.
Consequently, mere reassignment can reduce the mean of per-step peak primitive
work by **at most5.49% at256** under this fixed-bundle model. The true largest
block could be heavier, making the remaining room smaller. At160 the analogous
bound permits at most16.53% reduction; it does not establish attainability.
These are bounds on an operation-count model, **not bounds or predictions for
solver speedup**. Reassignment also changes locality/communication and may change
other work. Removing or repartitioning scaffold work changes the model itself.

This is evidence against a broad assignment-only performance campaign. Target
heavy units for finer subdivision and explicitly decide how their first-block
scaffold burden is handled. Splitting only interior subtrees can leave the same
scaffold hotspot. Preserve compact direct kernels and numerical ordering; avoid
reintroducing a per-value execution interpreter. More blocks also increase halo,
metadata, route and migration costs, so subdivision needs a bounded cost test.

## 3. Strong-scaling costs: additional limits beyond imbalance

The bare ordinary-step median rises505.522→544.532ms from160 to256 (+7.72%).
V7 is still8.09%/8.00% faster than266 at these counts. Relative efficiency for
160→256 is58.02%; this intermediate test does not qualify the >=640-rank goal.

| Coarse ordinary-step quantity, ms | 160 | 256 |
|---|---:|---:|
| RHS arithmetic, mean across ranks/steps | 27.967 | 25.399 |
| RHS arithmetic, mean of per-step rank maxima | 34.797 | 36.442 |
| RHS preparation, same slowest whole-window rank | 58.565 | 72.516 |
| RHS exchange, same slowest whole-window rank | 130.661 | 145.051 |
| Inverse MPI wait, same slowest whole-window rank | 38.761 | 58.731 |
| Complete profiled window, mean slowest rank | 506.479 | 557.076 |
| Logged ordinary mean minus profiled-window mean | 9.556 | 29.227 |

The first two rows measure arithmetic distributions; the later phase rows use
the same whole-window critical rank for each step. **Do not sum phase maxima**
or mix those arithmetic maxima with the critical-rank phase budget.

Arithmetic averages fall only9.18%, and per-step maxima rise4.73%, despite37.5%
less mean primitive work/rank. Aggregate arithmetic self-wall rises4.475→6.502
rank-seconds per ordinary step (+45.31%). Fixed loops, copies, irregular work,
cache/memory effects and costs not represented by primitive counts can contribute;
this is not a measured FLOP-efficiency loss or proof of bandwidth saturation.
Arithmetic self time excludes nested exchange but is not a pure hardware compute
measurement. Work imbalance alone does not explain its poor average scaling.

The arithmetic mean-of-max minus all-rank mean is6.83ms at160 and11.04ms at256.
Those gaps describe available balance in the measured phase under idealized
redistribution; they are not achievable full-step savings. Delayed producers can
also affect other ranks' waits, while changed placement can worsen locality.
Even the entire arithmetic phase is much smaller than the remaining140.40ms
bare main gap at160. RHS rebalancing by itself is therefore not a credible plan
for legacy parity.

At256 zero-Domain ranks average155.85ms RHS exchange and68.94ms inverse wait,
versus139.96/43.17ms on Domain owners. This is consistent with less-loaded ranks
waiting for dependencies, but MPI waits also contain transport/synchronization.
Do not label them network cost or promise their removal from work balance alone.
The profile also shows growing adaptation, native-wavelet, migration and catalog
costs. Remaining Domain ownership is a concrete architectural constraint, not
proof that any one of these phases is now the dominant cause.

The coarse run perturbs ordinary mean timing by4.52% at160 and9.26% at256.
The logged/profile-window difference is retained explicitly above; its growth
is not assigned to a cause. Two rank counts on one node do not identify a fixed
serial fraction, a network-latency floor or a memory-bandwidth saturation point.
No Amdahl fit or multi-node scaling forecast is justified here.

### Memory is a separate constraint

Existing per-epoch counters, averaged over ordinary RK stages, show:

| Summed named live extent | 160 | 256 | Change |
|---|---:|---:|---:|
| Compact program nodes | 139,126 | 144,308 | +3.72% |
| Numerical buffer bytes | 1,631,451,872 | 1,664,415,515 | +2.02% |
| Route payload bytes | 633,721,184 | 721,768,212 | +13.89% |
| Rounding column bytes | 103,197,189 | 103,197,189 | unchanged |

These narrow categories exclude descriptors, geometry/accepted workspaces,
construction peaks and MPI storage. They are not RSS and need not form a complete
nonoverlapping allocation census. Their modest growth cannot be used to explain
summed rank peak RSS rising43.79→61.30GiB (+40.00%). RSS sums asynchronous rank
peaks. The unaccounted storage/construction/allocator/MPI components remain to be
identified before asserting a memory bottleneck or a shared-memory remedy.

## Repair priority and decision gates

1. **Evaluate assignment and granularity together, but change them separately.**
   Use actual block operation vectors (including internal patches and scaffold
   targets), retained rank-local work and locality constraints. Compare current
   assignment with the existing minimax reference under the same weights, then
   compare revised weights, then subdivided heavy units. Record largest-unit and
   total-work/P bounds. The current evidence is already sufficient to rule out
   promising a major gain from assignment alone.
2. **Prepare finer numerical ownership for heavy units.** Preserve stable geometry
   identities, exact ownership/epochs and compact stencils. Treat scaffold work
   explicitly. Bound added halo/import volume, routes, plan construction, migration
   and live memory; do not let finer subdivision recreate the overhead just removed.
   Do not commit to an arbitrary uniform eight/sixteen-block split without costs.
3. **Remove retained whole-step dependencies alongside this work.** Map inverse,
   RK boundary adapters and lifecycle readers to their real owners and move a
   complete producer/consumer slice with paired retirement. Seek reductions in
   both ordinary P<=D time and P>D waiting; balancing already migrated RHS work
   does not finish independent timestep ownership. Continue to include remaps.
4. **Close only the named evidence gap.** Before choosing an actual partition or
   subdivision, a validation-only catalog capture needs per-block identity, leaf
   weight, owned operation/target counts including scaffold association, owner,
   and boundary/import footprint, plus retained Domain work by rank/phase. Prefer
   existing compiled plan counts; no per-value timers or hot-path profiling.
   This resolves whether/where to split and whether alternative placement pays
   its traffic cost. Current55504 has rank aggregates, so it cannot reconstruct
   the full J5 counterfactual. Fold this bounded capture into implementation
   qualification; do not launch another broad timing campaign merely for it.
5. **Require complete-step evidence before promotion.** Preserve v7's absolute
   gains and compare proposed ownership, weights and subdivision changes
   independently. Qualify optimized/checked/oracle science, reload/remap and
   zero-Domain work, then paired bare ordinary timing at83/160/256, allocation
   and memory accounting. A better proxy distribution is insufficient. Proceed
   to matched160/320/640 on bb only when granularity and ownership are ready;
   do not splice bb13 and bb timings or claim final efficiency without its gate.

`MPI_Comm_split_type` identifies node-local groups; it does not itself redistribute
work or share buffers. The current256-rank test is already on one node. Explicit
shared buffers would require evidence of duplicated immutable storage and
synchronization/NUMA costs, with bounded communicator/window lifetimes. It cannot
split heavy numerical blocks or remove retained Domain dependencies, so it is
not the leading repair justified by this assessment.

The next substantial work should therefore combine a defensible finer-work-unit
plan with retirement of remaining critical Domain dependencies, measuring their
changes separately. There is evidence for those architectural needs; there is
not yet a defensible percentage promise for the next complete-timestep speedup.
