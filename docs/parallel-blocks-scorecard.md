# Block refactor scorecard

Updated2026-10-14 after independent review of Stage280 job55508.
**Stage266 remains accepted. Stage280 repairs Stage279 construction costs and
improves ordinary medians3.67%/4.25%/3.74% versus v7 at83/160/256. Keep this repair
candidate and preserve v7/279 controls; no accepted-source promotion or further
ownership migration. Legacy parity, memory parity, whole-timestep independence
and useful scaling remain open.256 remains8.95% slower than160; RSS is still
3.75%/5.28%/5.60% above v7.**
Objectives: [authoritative reset](parallel-blocks-objectives.md).

Active development branch: `codex/stage278-structural-repair`, based directly on
accepted Stage266. [Branch roles](parallel-blocks-branches.md) preserve the separate
historical experiments. The qualified predecessor, v7 and reviewed Stage280 source are recorded as
candidate snapshots; production acceptance remains Stage266.

## Campaign budget rule — user instruction 2026-10-14

Future campaigns must use the minimum runs that can guide decisions on absolute
performance, scaling, repairs or architecture. Justify each configuration,
control, repetition and diagnostic by its decision value; reuse applicable
source-bound evidence and preserve necessary science/ownership gates. Add runs
only for a specific uncertainty that can change the decision. Avoid automatic
full repetition/profile matrices. State build/run counts and estimated cost
before delivery; predeclare bounded follow-up criteria where needed.

Job55508 completed its sealed60-run schedule in53:44 and is independently
reviewed below. This conservative matrix is not the default for future campaigns.
Earlier preparation/status statements are historical; no automatic promotion.

## Source and campaign decisions

| Item | Current state |
|---|---|
| Accepted solver | Stage 266 `ff740744608731a05d53771c32a3055309f3a586` |
| Legacy control | Main `6a245db7f1632b0205c2bebe95e5f0cd535c2582`, established checkpoint I/O hardening |
| Stage 273/277 | Separate unpromoted partition/cache candidates; excluded from 278 |
| Stage 278 v5 /55502 | Previous ownership candidate; science pass, performance superseded by v7 |
| Attribution 55503 | 12,099 files, 37 passing runs, 38:06; construction diagnosis informed v7 |
| Structural v7 /55504 | 15,672 files verified; five builds, 48 passing runs, 42:33; exact six-file repair delta |
| Fallback/cache proposals | Excluded; fallback lacks standalone case, cross-epoch route-cache probe had zero hits |
| Stage279 v3 /55507 | 19,118 files,49 passing runs,43:52; ownership/science pass,83/160 performance and memory regress; rejected for promotion |
| Stage280 v1 /55508 | 22,280 files,60 passing runs,six builds,53:44; construction cost repaired and absolute medians surpass v7; memory/main parity and scaling unresolved |

Implementation: `/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic`. Working source adds Stage279 and the two-file Stage280 repair to exact v7. Controls and all prior
sealed sources/evidence are preserved. Stage266 remains accepted; v7 remains the
cluster-qualified development control. User-submitted55508 is reviewed; no automatic promotion or new remote action.

## Latest decision — Stage280 job55508

[Full evidence and remaining repair decisions](parallel-blocks-stage280-job55508-results.md).
All60 runs pass science/checked/oracle and ownership/retirement gates. The exact
prepared kit/source seals match.96 zero-Domain ranks at256 perform transform
and RHS work; one-context lifetime and publication counts pass. Stage266 remains
accepted; retain280 for continued repair of the same slice, preserve v7/279
controls. No new campaign or migration is justified before a concrete cost case.

Mean independent bare ordinary run medians, ms:

|Tasks|Main|266|V7|Stage279|Stage280|280 vs v7|
|---|---:|---:|---:|---:|---:|---:|
|83|545.787|756.803|631.718|703.699|608.564|-3.67%|
|160|356.946|546.019|505.146|568.446|483.675|-4.25%|
|256|unsupported|585.611|547.472|553.032|526.969|-3.74%|

All four repair medians lie below all four v7 medians at each count. Recovery
versus279 is13.52%/14.91%/4.71%; versus26619.59%/11.42%/10.01%. Remaining main
gaps62.78/126.73ms (11.50%/35.50%). Ordinary mean gains versus v7 are smaller:
2.12%/2.47%/0.35%. At256 the means overlap, one repair run has retained ordinary
outliers and non-remap means worsen0.52%; do not claim repeatable throughput gain.

Coarse transform-plan self times fall132.60→34.12,116.37→27.13,67.95→30.60ms.
This supports the structural construction repair. The global kernel census,
numerical traffic, messages and runtime route counts match279 exactly. Exchange
self times are largely unchanged; no generic transport/batching conclusion.
Partial named transform buffers fall18.06%/14.30%/13.32%, but total summed rank
peak RSS remains25.877/46.067/64.895GiB,3.75%/5.28%/5.60% above v7 and only
0.83%/0.57%/0.36% below279. Peak sums are not simultaneous allocation peaks.

Scaling fails:280 medians483.68→526.97ms (+8.95%), speedup0.91784 and relative
efficiency57.37%; means480.04→536.11ms (+11.68%), efficiency55.96%. V7 scaling
is slightly better. No>=640 qualification or whole-timestep closure.

Existing ordinary counters/profile evidence confirms unchanged assignment and
granularity: mean step work max/mean1.354→1.989 at160→256, peak primitive work
falls only8.07% while average falls37.5%, and maximum RHS arithmetic rises
35.03→37.00ms. Fixed-bundle reassignment alone permits at most5.49% peak-count
reduction at256, not a runtime forecast. Domain owners retain77.22% of primitive
work before retained Domain costs. RHS preparation grows58.29→70.72ms and
exchange130.94→151.44ms; these are coarse costs, not pure transport attribution.

**Next work:** use existing evidence/source to design structural RHS
construction/adapter repair and finer numerical/scaffold units with explicit
invalidation, memory and halo costs. Identify retained/replicated memory before
adding representations. Quantify removed/replacement work before choosing the
next implementation; no broad repeat campaign, automatic full control/profile
matrix, generic shared-memory shortcut or J6/J8 substitution. Required final
lifecycle, memory/main parity and matched>=640 scaling remain open.

The following Stage280 preparation/local notes and Stage279 rejection are
historical evidence, superseded by the cluster decision above.

## Stage280 preparation and local evidence — Stage280

[Implementation, contract and quantitative evidence](parallel-blocks-stage280-prebound-repair.md).
Stage280 repairs the same provisional ownership slice:16 ordinary routes bind
sole canonical providers directly, avoiding128 redundant preparation collectives
per epoch. Outer-vector and fixed-input writer proofs remain. Publication uses
local canonical rows or one unique remote row, removing duplicate numerical
storage/staging. Physical key lookup is mixed; topology/mask/state invalidation,
within-step reuse and one-context lifetime remain. Only two solver files change.

All11 local full replays and48 MPI exchange cases pass, including checked and
independent oracle. Exact-source J5 GCC13/MVAPICH checked/optimized builds pass.
The source and prior candidate/control evidence are sealed and preserved.

Local4-rank J4/J6 A/B/C/C/B/A bare comparison(two runs each): mean ordinary run
medians v7 3.866867s, Stage279 3.890662s, repair3.681480s. Repair is4.79%/5.38%
faster than v7/279 locally; both repair medians are below every control median.
Ordinary means improve5.55%/6.39%. This is local evidence, not cluster qualification.
Named transform numerical/route allocations601.55→423.26MiB(-29.64%); total RSS
and memory parity remain unqualified. Global kernel census, remote numerical
payload and runtime route invocation counts are exactly unchanged.

One local coarse pair shows plan165.24→157.13ms and exchange212.19→178.82ms,
but RHS exchange increases513.54→585.74ms and whole profiled ordinary means change
only3.81676→3.80259s. Do not dismiss shifted waiting or attribute the complete bare
gain to one mechanism. Source-derived128-call removal is a cost case, not a
latency-derived speedup forecast. Preserve v7/279/266/main controls and require
83/160/256 ordinary-time, RSS, production/oracle and scaling validation before
promotion. No cluster solver or campaign was run for this repair; only an isolated
two-source compiler check was uploaded. Final lifecycle, granularity and>=640
qualification remain required. Stage266 remains accepted.

The Stage280 bbserv kit is sealed and ready for user submission: six serial
builds,11 validations,49 measurements,60 total runs on exclusive bb13. Forty bare
runs retain v7/279/266/main controls; nine separate coarse runs test attribution.
42 current helper tests and fresh-extraction/source checks pass. Ordinary time,
RSS and160-to256 scaling remain pending; no upload/submission or promotion.
Planning60-80 minutes within90, excluding queue/collection.

The job55507 rejection below is the established cluster result; it is not erased
by local repair evidence. No additional ownership migration was stacked.

## Latest decision — Stage279 job55507

[Full review and attribution](parallel-blocks-stage279-job55507-results.md).
All49 runs pass science, exact-source/fixture checks and ownership/retirement gates.
At256 all96 zero-Domain ranks have transform records and real RHS work. This is
provisional ownership progress; final lifecycle and whole-timestep closure remain open.

Mean of bare run ordinary medians, ms:

| Tasks | Main | 266 | V7 | Stage279 | Stage279 vs V7 |
|---|---:|---:|---:|---:|---:|
| 83 | 540.209 | 754.080 | 634.078 | 683.278 | +7.76% |
| 160 | 359.166 | 547.708 | 504.941 | 569.310 | +12.75% |
| 256 | — | 593.329 | 550.577 | 555.903 | +0.97% |

Every candidate median is slower than every v7 median at83/160;256 ranges overlap.
Ordinary means regress13.70%/14.62%/6.00%, with slower remap and non-remap means.
Candidate is now3.94% slower than266 at160 and26.48%/58.51% slower than main83/160.
The local slowdown warning is confirmed. No checkpoint-frequency explanation.

Median160→256 speedup improves0.9171→1.0241, relative efficiency57.32%→64.01%,
largely because160 worsens64.37ms while256 worsens5.33ms. Candidate ordinary means
still worsen564.709→570.450ms. No absolute256 gain or scaling closure;>=640 remains
untested. Mean summed rank peak RSS rises24.92→26.09,43.71→46.32,61.21→65.10GiB
(+4.71%/+5.97%/+6.35%); these are not simultaneous memory peaks.

Separate coarse critical-rank self times: new transform total187.66/164.02/119.16ms,
old related envelope reductions115.96/85.57/108.33ms. New preparation alone costs
128.83/117.23/68.13ms. Audited segments locate most preparation in owner-home
resolution and route/adapter binding, not geometry arithmetic. All12 trajectory
states trigger new topology/mask epochs; within-step reuse already works. No
unsafe count-only reuse or generic transport/shared-memory conclusion follows.
Profiles and bare medians are different statistics; the report preserves that distinction.

**Next gate:** address measured construction and duplicate-storage costs before
another campaign. Require a defensible net cost reduction; do not stack another
ownership migration or repeat merely to confirm the clear regression. Regaining
v7 requires49.20/64.37ms at83/160; main requires143.07/210.14ms. These are target
gaps, not projected savings. Preserve Stage279 numerical/oracle evidence; retain
v7 as development control and266 as accepted solver. No solver change or remote
action was made during review.

The following preparation notes and job55504 tables are historical evidence.

## Stage279 local implementation status (before job55507)

[Implementation and evidence](parallel-blocks-stage279-resident-transform.md).
Three provisional resident RK → wavelet/compression → inverse → next-RHS chains
now run on numerical owners, retiring the corresponding old WT paths and canonical
Domain interior publication. Raw-RK boundary/fixed-scaffold adapters and final
lifecycle work remain. Direct canonical reads replace a redundant provider-plan
construction; sequential exchanges share maximum-sized scratch and one epoch
communicator. The640-block catalogue and assignment policy are unchanged.

All nine final local four-rank J4/J6 replays pass, including production optimized,
checked and independent oracle;24 MPI exchange cases and32 helper tests pass.
Local tests do not cover zero-Domain ranks. The cluster kit includes a dedicated
256-rank production audit in addition to oracle validation.

Local A/B/B/A mean run medians: v7 6.802392s, candidate6.945858s (+2.11%).
Ordinary means:6.879176s versus7.190612s (+4.53%). Run ranges overlap amid heavy
host load; no demonstrated speedup. Increased coarse RHS wait could include
new pipeline producer lateness; do not dismiss it as unrelated or infer pure
transport. Candidate-specific named numerical/route allocation falls58.52%
versus its initial implementation, not versus v7 and not total RSS. The retained
final workspace and compiled metadata still contribute memory.

The sealed J5 kit runs49 cases on exclusive bb13: nine validations,34 bare and
six coarse runs, five serial builds. Stage279 performance, memory,256 scaling
and legacy parity await this comparison. The historical cluster gates below
remain the latest established results. No further cutover before review.

## Job55505 build failure and v2 correction

Job55505 failed after4:32 in the GCC13 checked build, before any solver runs.
This is a portability failure, not performance evidence. GCC13 flagged the
allocatable index-selection self-reassignment; v2 replaces it with explicit
allocation and stable ordered filling, including empty selections. Checks remain
fully enabled. [Correction and evidence](parallel-blocks-stage279-gcc13-fix.md).

The exact v2 candidate passes full isolated GCC13.2/MVAPICH4 J5 checked and
optimized builds, three fresh local4-rank optimized/checked/oracle replays and34
helper tests. V1's local timing observations remain historical evidence only;
no new speedup or scaling result is available. V2 retains the same49-run campaign
and unchanged v7/266/main controls. The failed kit is preserved. The full v2 kit
is ready for user transfer/submission; only an isolated one-source compiler
check was uploaded and built remotely. No solver or Slurm job was launched.

## Job55506 validation-checker failure and v3

Job55506 completed five builds and six solver validation runs in10:28. The first
five passed in the runner; the sixth (83-rank production) passed science then
stopped on an incorrect requirement that every canonical owner have all three
transform kernel classes. Four ranks legitimately had outer-only descriptors.
No measurement,160-rank or256-rank run began. This is not a performance result.

[Checker correction and evidence](parallel-blocks-stage279-checker-fix.md).
V3 requires positive total records per canonical owner and preserves all existing
science/epoch/traffic/retirement gates. It adds global class-census equality across
all36 transitions and validation counts. All six completed cluster runs and three
local audits pass reanalysis; checkpoint/binary/VTK hashes pass read-only checks.
40 helper tests include the actual83-rank regression. Solver source, flags and
fixture remain exact v2, retaining its GCC13/full-local qualification. V3 awaits
user transfer/submission with the unchanged49-run schedule. Failed originals are
preserved; no new solver or Slurm run was performed during this correction.

## Historical v7 dual-goal gates — job55504

| Gate | Evidence / remaining work |
|---|---|
| Connected RHS and real RK ownership | Pass in candidate; all 96 zero-Domain ranks do numerical work, 16,212,000 primitive events at 256 |
| Paired retirement | Domain compatibility RHS, scalar capture/transport, velocity remainder export, full temperature-stage publication and displaced interior RK remain retired; v7 adds no numerical retirement |
| Science/lifecycle | All 48 runs pass; exact accepted-family/main policy, 12 states, 2 reloads, 2 remaps, 62 VTK members; checked/oracles retained |
| Whole-timestep independence | Open: boundary/scaffold/RK adapters, inverse/lifecycle publication, split Simple physics, geometry/workspaces |
| Absolute performance vs 266 | Intermediate gate surpassed on this fixture: v7 is 16.32%/8.09%/8.00% faster at 83/160/256; every candidate median below every 266 control median |
| Legacy P<=D performance | Open: v7 remains 17.03%/38.45% slower than main at 83/160; further 14.55%/27.77% candidate time reductions required |
| Intermediate P>D scaling | Unresolved:160→256 slowdown 7.72%, speedup 0.9284, relative efficiency 58.02%; previous 61.48%, 26658.08% |
| Memory | Improved 9.84%/18.97%/24.22% vs previous, but 16.00%/38.67%/42.42% above 266 summed rank peak RSS |
| Useful P>=4 D | Unqualified; J5 requires >=640 matched multi-node ranks and sufficient subdivision/load balance; no invented efficiency threshold |

Scope: restarted compressible non-mode-split RK4 climate/Simple physics, no
NCAR/SSO, established short J5 fixture on nk/bb13. Ownership credit is for the
qualified slice, not an independent whole timestep or all solver configurations.

## Historical matched ordinary-step performance — job55504

Milliseconds, mean of four previous/v7 and two main/266 run medians. Ordinary
steps exclude startup/reload/saves and retain remaps; steps are not independent
repetitions. Production saves remain rare(O1000 steps), remaps about every five.

| Ranks | Main | 266 | Previous | V7 | V7 vs previous | V7 vs 266 | V7 vs main |
|---|---:|---:|---:|---:|---:|---:|---:|
| 83 | 540.900 | 756.425 | 790.529 | 632.989 | -19.93% | -16.32% | +17.03% |
| 160 | 365.127 | 550.008 | 598.310 | 505.522 | -15.51% | -8.09% | +38.45% |
| 256 | — | 591.899 | 608.270 | 544.532 | -10.48% | -8.00% | — |

Every adjacent bare A/B pair favors v7. Ordinary means improve 19.04%/14.62%/
10.96%; remap and non-remap groups improve at each count. Gains are supported
within this allocation; long representative workloads and multi-node performance
are not qualified. Remaining main gaps 92.09/140.40 ms cannot be declared solved
by the now much smaller preparation budget.

## Attribution, scaling and memory

Route construction falls 151.51→37.33, 92.41→33.57, 79.95→40.72 ms at 83/160/256;
input plans 33.83→14.28, 27.15→14.70, 31.02→19.15 ms. Named preparation self subtotal
falls 216.94→62.74,141.50→58.56,133.12→72.52 ms. RHS exchange is almost unchanged
at 160/256. This supports structural construction/hash repair rather than a
transport/batching explanation. Individual mechanisms were not ablated.

V7 is faster at 256 than before, but gains more at 160:92.79 vs 63.74 ms saved.
This increases the 160→256 median gap from 9.96 to 39.01 ms. The scaling failure
persists; neither better absolute time nor correct zero-Domain work closes it.
Candidate ordinary-mean scaling also fails(speedup 0.9201). Coarse window growth
includes inverse wait, RHS exchange, adaptation, native wavelet, migration and
catalog/route construction. Wait time includes late producers; causality is not
isolated. Coarse mean perturbation reaches 9.26% at 256; do not partition bare
median growth with these profiled means or omit unprofiled logged-step time.

| Ranks | Previous RSS GiB | V7 RSS GiB | 266 RSS GiB | V7 vs previous | V7 vs 266 |
|---|---:|---:|---:|---:|---:|
| 83 | 27.65 | 24.93 | 21.49 | -9.84% | +16.00% |
| 160 | 54.04 | 43.79 | 31.58 | -18.97% | +38.67% |
| 256 | 80.89 | 61.30 | 43.04 | -24.22% | +42.42% |

RSS is summed individual rank peaks, not concurrent job peak. Lower memory is
measured; MPI/solver allocation attribution and long-run growth are unresolved.
All 1,094 integer work/traffic/narrow-storage rank files exactly match 55503.
Existing v7 source lifetime limits and metadata payload tests remain valid.

## Load balance and granularity evidence — job 55504 follow-up

[Detailed derivation](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/review-55504/load-balance-followup.md),
[rank statistics](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/review-55504/load-balance-followup.json). The same 71,052,000 trajectory
primitive events are spread more unevenly at256: maximum/mean rises1.24→1.67.
Only640 blocks exist; average blocks/rank falls4→2.5, and122 of256 ranks have one
final block. The maximum rank ordinary RHS arithmetic mean barely changes:
33.78→33.51 ms. At256 Domain owners average342,750 primitive events per rank,
zero-Domain owners168,875; the latter wait more in RHS exchange and inverse.
These are explicit load-distribution/granularity signals, not just transport cost.
Counters cover the validation trajectory; coarse timings cover ordinary steps.

The [completed ordinary-step assessment](parallel-blocks-stage278-balance-assessment.md)
separates assignment quality, granularity and strong-scaling costs. All nine
ordinary steps have identical global work at83/160/256. Mean per-step maximum/mean
primitive work rises1.354→1.989 at160→256; maximum work falls only8.07%, versus
37.5% less mean work/rank. This supersedes trajectory aggregates for ordinary-step
imbalance: averaging before taking the maximum hides changing critical ranks.

At256 the busiest rank owns only one block on three ordinary steps. With existing
block/scaffold work bundles intact, single-block observations bound the possible
reduction in mean per-step peak primitive work from reassignment alone to at most
5.49%. This is a conditional work-model bound, not a solver speedup forecast.
Heavy-unit subdivision and explicit treatment of first-block scaffold work are
needed for substantial flexibility; extra metadata/halo/memory costs must be bounded.

Mean per-step maximum RHS arithmetic time rises34.80→36.44ms, while all-rank mean
falls only27.97→25.40ms. Preparation and waits also grow. The profile cannot isolate
network transport or identify a serial fraction. Narrow numerical buffers grow
only2.02% and route payload13.89%, whereas summed rank peak RSS grows40.00%; memory
attribution remains incomplete. Arithmetic balance alone is insufficient for the
140.40ms remaining main gap at160.

Next evaluate actual work weights, assignment and finer units separately, alongside
paired retirement of retained inverse/lifecycle dependencies. A bounded validation
capture of per-block work/scaffold/import costs is needed to choose actual J5 cuts;
existing rank aggregates cannot reconstruct that counterfactual. No broad diagnostic
campaign or automatic partition/cache adoption follows from this assessment.
Stage273 was locally qualified and unpromoted in the reviewed preparation record;
this assessment has no273 cluster rejection evidence. Stage277 is a separate,
unpromoted cache experiment with small mixed cluster timing changes.

## Earlier repair specification — superseded by job55507 decision

- Retain v7's demonstrated gain and ownership retirements; do not repeat the same
  campaign merely to reconfirm a large benefit.
- Validate the implemented [resident transform repair](parallel-blocks-stage279-resident-transform.md)
  on bbserv. The implemented slice closes resident RK → forward wavelets/compression → inverse → next RHS, including
  actual boundary/scaffold versions at the three provisional RK4
  transitions as a complete producer/consumer slice; final lifecycle closure
  remains a separate required gate. Preserve direct compact column execution.
- Verify the selected retirement of geometry-owner vector/inverse work and its coefficient,
  initial-gather and final-installation transfers together. Do not credit already
  removed in-phase gathers or the deferred first installation; current canonical
  Domain publication is local and full SOL network returns are already zero.
- Evaluate finer, cost-aware numerical units separately from this ownership cutover.
  Include internal/coarse patches and first-block scaffold work in the subdivision
  contract; a partition-only repair has insufficient modeled room at256.
- The broad transform/boundary/publication profile envelope is208.11/142.25/152.56ms
  at83/160/256, including necessary work and callers beyond the selected slice.
  It is not a savings forecast. At160 it barely exceeds the140.40ms main gap, so
  this slice alone cannot defensibly be promised to close the gap. Retain RHS
  preparation/execution as a separate residual target; no transport-dominance claim.
- Avoid another metadata-only plan as the whole answer to main parity. At 160 the
  remaining preparation subtotal 58.56 ms is much smaller than the 140.40 ms bare
  main gap (different statistics; not an additive forecast).
- Require an architecture repair that addresses 256 scaling and memory while
  preserving ordinary absolute performance. Generic batching/shared memory and
  count-only cross-epoch reuse remain unsupported causal shortcuts. Any extra
  instrumentation needs a named unresolved implementation decision.
- Finish whole-timestep ownership, granularity/load balance and matched >=640
  qualification as separate required gates. Stage 266 remains accepted.

## Workload continuity — user clarification 2026-10-14

J5/D=160 remains the deciding83/160/256 performance track and the >=640-rank
block-independence target. A possible future J6-coarsest/J8-finest test is additional,
not a substitute: D would become640, so160/256 would both be below the Domain
count,640 would only reach it, and P>=4D would require2,560 tasks. No larger-workload
campaign is scheduled. Increasing work may reveal strong-scaling headroom but
must not hide the original Domain barrier.

## Evidence and history

[55504 full report](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/docs/parallel-blocks-stage278-job55504-results.md),
[verified return](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/review-55504/verified.json), [metrics](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/review-55504/metrics.json),
[prior scorecard](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/review-55504/scorecard-before-review.md),
[runtime ledger](parallel-blocks-campaign-runtimes.json).
[55503 attribution report](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/docs/parallel-blocks-stage278-job55503-results.md)
and [v7 local qualification](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/docs/parallel-blocks-rhs-structural-repair-278-v7.md)
retain prior evidence and limitations. Original returned payloads are not rerun;
raw reports/hashes are independently verified. No remote actions during review.
