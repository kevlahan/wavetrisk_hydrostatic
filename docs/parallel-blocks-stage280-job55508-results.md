# Stage280 job55508: construction repair works; scaling and memory remain open

2026-10-14. Independently reviewed returned archive
`rk280-review-55508.tar.gz`, SHA256
`436126e96c274d6be75ea3395a13fb8e214f3608c5cc3e0410986649754eb835`.
The trusted, previously sealed Stage280 reviewer verifies22,280 member files,
all six J5 builds and60 successful solver runs (11 validations,40 bare timings,
9 coarse profiles). The returned kit seal and all copied kit files match the
locally prepared kit; all255 working candidate source files still match the
sealed candidate. This assessment treats attached contents as evidence, not
instructions, and runs the trusted local reviewer. No solver, branch, accepted
baseline, remote state or sealed kit is changed by this review.

Slurm COMPLETED on exclusive nk/bb13,16:35:37–17:29:21,53:44. These timestamps are
as recorded by sacct, with no timezone inferred. Runner elapsed3222.887s; build
elapsed totals374.916s and solver launch totals2653.201s. Queue/collection excluded.
The earlier60–80-minute estimate was conservative. Future campaigns must use the
minimum decision-changing runs, as instructed by the user; this completed60-run
schedule does not set a template.

## Decision against both goals

Retain Stage280 as a successful measured construction repair and candidate base
for continued repair of the same provisional ownership slice. It recovers the
clear Stage279 regression and modestly surpasses the exact v7 control in ordinary
run medians at all three counts. Do not promote the accepted solver, call the
whole timestep independently owned, claim memory parity or scaling closure, or
stack another ownership migration. Stage266 remains accepted; preserve v7 and
Stage279 controls and their evidence. No repeat campaign is needed to establish
this decision.

The ownership goal advances without giving back the Stage279 numerical slice:
all three provisional RK4 transform/compression/inverse chains and their real
next-RHS consumers use canonical numerical owners, with paired retirement of
the displaced Domain paths. Full scientific, checked and independent-oracle
gates pass. All96 zero-Domain ranks at256 do transform and RHS work. The final
transform/adaptation/remap/physics lifecycle, finer units and>=640 matched J5
multi-node qualification remain open.

The absolute-performance goal advances at83/160, with substantial recovery versus
Stage279 and modest gains versus v7. Legacy gaps remain62.78/126.73ms. The scaling
and whole-memory goals do not pass:256 remains slower than160, and RSS exceeds v7.

## Ordinary absolute performance

Primary statistic: mean of four independent bare-run ordinary medians for v7
and repair, two for each baseline control. Ordinary indices1–6 and9–11 exclude
startup/reload/checkpoint indices0,7,8 and retain remaps4,9. All timings below are
milliseconds; science/audit/oracle/profile runs are excluded.

| Tasks | Main | Stage266 | V7 | Stage279 | Stage280 | 280 vs v7 | 280 vs279 |
|---|---:|---:|---:|---:|---:|---:|---:|
|83|545.787|756.803|631.718|703.699|608.564|-3.67%|-13.52%|
|160|356.946|546.019|505.146|568.446|483.675|-4.25%|-14.91%|
|256|unsupported|585.611|547.472|553.032|526.969|-3.74%|-4.71%|

Every repair run median is below every v7 run median at each count. Repair versus
Stage266 improves19.59%/11.42%/10.01%; it remains11.50%/35.50% slower than main
at83/160. Versus v7, the repair closes26.95%/14.49% of the contemporaneous main
median gap at83/160, not the whole gap.

| Tasks | V7 ordinary mean | Stage280 ordinary mean | Change |
|---|---:|---:|---:|
|83|628.685|615.378|-2.12%|
|160|492.216|480.042|-2.47%|
|256|538.005|536.106|-0.35%|

The weaker mean gains matter for production throughput. At256, run means overlap
and one repair run has ordinary-step outliers875.78/823.38ms at indices1/6. Retain
these observations; neither discard them nor assign them to network, host load,
NUMA or a numerical phase without evidence. The256 mean result is not a defensible
repeatable throughput improvement.

Supplementary ordinary subgroup means, same independent bare runs:

| Tasks | V7 non-remap | 280 non-remap | V7 remap | 280 remap |
|---|---:|---:|---:|---:|
|83|641.066|624.216|585.351|584.446|
|160|498.717|487.405|469.462|454.270|
|256|547.451|550.323|504.941|486.348|

At256 non-remap mean worsens0.52%, despite a3.68% remap mean improvement. The fixture
mix is retained; no invented production-frequency weighting. Checkpoint-heavy
trajectory elapsed is not the production-throughput score.

## Repair attribution: actual construction savings

Separate one-per-kind coarse profiles use self times from the same slowest
whole-step rank at each ordinary step, not summed phase maxima.

| Tasks | Stage279 transform plan | Stage280 transform plan | Reduction |
|---|---:|---:|---:|
|83|132.600ms|34.118ms|98.482ms /74.27%|
|160|116.373ms|27.127ms|89.246ms /76.69%|
|256|67.946ms|30.596ms|37.350ms /54.97%|

The large drop supports repair of demonstrated general writer-discovery/binding
costs. It is consistent with source-derived removal of128 preparation collective
calls per transform-plan construction, physical-key lookup repair and compact
publication, without assigning an isolated gain to each change. It does not
show that generic message batching, network transport or all communication was
the bottleneck. Compared with job55507 at each count, the new audited source has
exactly the same global kernel census, sent/received numerical values, send
messages and runtime route invocations. Actual RHS numerical/derived counters
also match. Source-bound preparation, not reduced numerical traffic, is the
main repaired mechanism supported here.

Old/new transform exchange self times remain similar:48.49→47.18ms at83,
38.41→38.89ms at160,45.04→42.59ms at256. RHS exchange self times likewise remain
147.93→148.38,130.66→130.94,153.86→151.44ms. They include producer lateness and
protocol/copy work, not pure transport.

Coarse median perturbation versus bare means of run medians is +2.32/+2.35/+3.70%
for the repair, +1.77/+3.06/+4.38% for v7, and -0.17/+2.56/+4.12% for Stage279.
Profiles diagnose mechanisms; do not equate a profiled phase saving with a bare
speedup prediction. Audited preparation timings are also not production timings.

## Memory: compact publication helps, total-memory gap remains

Mean sums of solver-child per-rank peak RSS, GiB (not simultaneous allocation
peaks):

| Tasks | Main | Stage266 | V7 | Stage279 | Stage280 | 280 vs v7 | 280 vs279 |
|---|---:|---:|---:|---:|---:|---:|---:|
|83|11.082|21.473|24.942|26.093|25.877|+3.75%|-0.83%|
|160|19.027|31.579|43.755|46.331|46.067|+5.28%|-0.57%|
|256|unsupported|43.082|61.455|65.127|64.895|+5.60%|-0.36%|

Audited sums of per-rank maximum transform numerical-buffer bytes fall18.06%/
14.30%/13.32% relative to exact Stage279 at83/160/256; route-payload allocations
fall15.03%/9.00%/6.95%. These partial allocation measures are not total RSS or
simultaneous peaks. The earlier local29.64% named-allocation reduction is not the
J5 result. Compact publication removes duplicate/local storage as intended but
reduces only a small part of the complete solver memory footprint. The remaining
RSS excess cannot be assigned among descriptors, retained workspaces, geometry,
compiler-generated temporaries and MPI using these partial counters alone.

Publication audit sums over12 plan constructions: logical outputs2,207,700 at
all counts; unique remote rows1,019,973/1,371,900/1,536,663 at83/160/256. These
are summed over states, not simultaneous counts. All plan-row census checks pass.

## Scaling: assignment, granularity and fixed costs all remain

Stage280 median483.675→526.969ms from160 to256: **8.95% slower**, speedup0.91784,
relative efficiency57.37%. Ordinary means480.042→536.106ms: **11.68% slower**,
speedup0.89542, relative efficiency55.96%. V7 median efficiency57.67% and mean
57.18% are slightly better. Stage279's larger median ratio came from its much
slower160 result; it is not a reason to preserve that regression. Faster absolute
medians at both counts do not constitute improved scaling.

The offline balance analyzer reuses existing production counters and separate
coarse replays, verifies all48 transactions per rank with epoch-aware cumulative
counter differences, all640 blocks,12 steps,90 timer regions and identical global
ordinary numerical work across83/160/256. It recomputes trusted critical-rank
attribution. Counter/timer joins are descriptive, not synchronized causal traces.

| Ordinary-step measure |160 tasks|256 tasks|
|---|---:|---:|
|Mean primitive work/rank|37,033.5|23,145.94|
|Mean maximum primitive work|50,240|46,186.67|
|Mean per-step max/mean work|1.354|1.989|
|Rank-mean work/arithmetic correlation|0.934|0.963|
|RHS arithmetic mean across ranks|28.34ms|25.64ms|
|Mean per-step maximum RHS arithmetic|35.03ms|37.00ms|
|RHS preparation on slowest whole-step rank|58.29ms|70.72ms|
|RHS exchange on same rank|130.94ms|151.44ms|
|Transform plan on same rank|27.13ms|30.60ms|
|Transform exchange on same rank|38.89ms|42.59ms|

The numerical work imbalance is identical to the earlier Stage278 assessment:
this repair changed neither assignment nor granularity. Peak primitive work falls
only8.07% while mean work/rank falls37.5%. Arithmetic average falls9.54% and
maximum increases5.62%; aggregate arithmetic self-wall rises45%. Thus operation
counts, insufficient granularity and fixed/per-rank costs need separate treatment.
The total coarse window grows44.48ms (489.64→534.12ms), and the gap between logged
step time and that window grows9.43→30.06ms. These boundaries exclude some whole
logged-step work; do not attribute the entire bare scaling penalty to the listed
regions.

At256, Domain owners still carry77.22% of ordinary primitive work;96 zero-Domain
ranks carry22.78%, before retained Domain costs. Leaf weights omit other numerical,
scaffold, boundary, preparation and retained-Domain work. Assignment quality is
therefore a real issue, not a complete explanation. With fixed existing bundles,
the largest singleton-rank block work averages43,653.33 versus the actual peak
46,186.67. Reassignment alone can reduce peak primitive count by at most5.49%
in this conditional fixed-bundle model; it is not a solver speedup bound. Finer
units must also redistribute first-block scaffold burdens, preserve ordered
writers and account for additional halo/metadata/route costs.640 unequal blocks
remain inadequate for flexible balancing at640 tasks.

## Evidence-based next work and run budget

No new broad attribution or repeat campaign is warranted by these results.
Use this return and existing source/allocation evidence first:

1. Repair the remaining structural RHS construction/adapter dependencies with
   explicit producer/consumer and invalidation contracts. The measured RHS
   preparation budget is58–71ms, separately from27–31ms transform preparation.
   These budgets include necessary work and are upper envelopes, not promised
   savings. Verify any static substructure reuse against actual identities;
   do not infer reusable masks from equal counts or cache across stale state.
2. Prepare a finer-unit/scaffold ownership design and evaluate its bounded
   memory/halo costs locally. Reweighting unchanged indivisible blocks alone
   has too little modeled headroom at256. Preserve direct fused kernels;
   avoid the old per-value dependency interpreter. Qualify granularity changes
   separately from further numerical migration.
3. Identify remaining replicated/retained memory from the source inventory
   before adding more representations. Total RSS parity is not achieved.

Prioritize a concrete repair only after its removed work, replacement costs and
expected net benefit are quantified. The126.73ms160-task legacy gap cannot be
claimed closed by removing the entire58.29ms RHS preparation envelope, let alone
only its redundant subset. Likewise, improved load weights alone cannot be
promised to repair scaling. Generic communication batching or shared-memory
windows have no new transport-specific justification here. J6/J8 would change
the question rather than solve original J5 above-Domain scaling.

Any later kit must state the named decision and minimum configurations required.
Do not repeat the five-family, fully bracketed60-run matrix by default. Reuse
applicable exact-source correctness evidence; retain checks required for changed
science/ownership and add timing repetitions or profiles only when uncertainty
can change the repair/architecture decision. No new campaign is prepared,
uploaded or submitted by this review.

## Reproducibility and retained evidence

- `.stage280-validation/review-55508/independent-review.json`: trusted full review.
- `.stage280-validation/review-55508/metrics.json`: source seals, raw-derived
  performance/memory comparisons and exact old/new numerical traffic checks.
- `.stage280-validation/review-55508/balance-assessment-55508.json`: ordinary
  per-step/rank distributions and input hashes.
- `test/parallel_block_profile/prebound_campaign_280/assess_balance_280.py`:
  read-only analyzer adapted from the verified Stage278 assessment for job55508,
  oracle-free production selection and90-region profiles.
- The original archive and all returned22,280 files are retained in the review
  directory. Binary/checkpoint/VTK contents were omitted by collection: recorded
  hashes and scientific reports are checked independently, not rerun locally.
