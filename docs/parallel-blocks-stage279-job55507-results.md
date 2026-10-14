# Stage279 job55507 — ownership passes, performance repair rejected

Reviewed2026-10-14. Stage266 remains accepted. Retain Stage278/v7 as the
cluster-qualified development control; do not promote or stack another migration
on Stage279. Preserve the experimental implementation and its numerical evidence.
This review changes no solver source and submits no remote work.

## Verification and scope

19,118 archive files verify; source, fixture and helper identities match the
locally sealed v3 kit. Five builds and all49 scheduled runs pass: nine validations,
34 bare runs and six coarse profiles. Slurm COMPLETED,14:55:26–15:39:18 EDT,
43:52 (2632s); runner2630.446s. Queue and collection are excluded. Exclusive
nk/bb13, J5/D160; main has no supported256-rank control. No J6/J8 substitution.

All12 states,2 reloads,2 remaps,62 VTK members and exact accepted-family/
main-relative science policies pass. Checked and independent oracles pass.
The three provisional transforms are executed on numerical owners, and the
corresponding old Domain paths remain retired. All96 zero-Domain ranks have
transform descriptors at256 (956,031 across the trajectory), alongside16,212,000
RHS primitive events. Transform descriptors include structural/inactive records;
they are not active FLOPs or calibrated work weights. Whole-timestep independence,
final lifecycle closure, finer subdivision and efficient>=640 remain unfinished.

## Matched bare ordinary performance

Milliseconds, mean of four independent v7/candidate run medians and two main/266
controls per supported count. Nine ordinary steps exclude0,7,8 and retain remaps4,9.
Timesteps are not independent repetitions. Profiling/oracles/audits are excluded.

| Tasks | Main | 266 | V7 | Stage279 | vs V7 | vs main |
|---|---:|---:|---:|---:|---:|---:|
| 83 | 540.209 | 754.080 | 634.078 | 683.278 | +7.76% | +26.48% |
| 160 | 359.166 | 547.708 | 504.941 | 569.310 | +12.75% | +58.51% |
| 256 | — | 593.329 | 550.577 | 555.903 | +0.97% | — |

All four candidate medians exceed all four v7 medians at83 and160: the regressions
are clear within this allocation. At256 the ranges overlap; median+0.97% does not
establish a repeatable small regression, but it certainly demonstrates no gain.
Candidate ordinary means are worse by13.70%/14.62%/6.00% at83/160/256. Both remap
and non-remap subgroup means are slower at every count. The regression cannot
be dismissed as rare checkpoint overhead. Versus266, medians are9.39% faster
at83,3.94% slower at160 and6.31% faster at256. Main gaps widen to143.07/210.14ms.

| Tasks | V7 median range (ms) | Stage279 median range (ms) | V7 mean (ms) | Stage279 mean (ms) |
|---|---:|---:|---:|---:|
| 83 | 629.404–640.014 | 664.214–719.889 | 627.878 | 713.896 |
| 160 | 500.613–506.969 | 565.752–573.283 | 492.689 | 564.709 |
| 256 | 545.063–556.919 | 547.419–562.445 | 538.180 | 570.450 |

## Scaling and memory

Median160→256 speedup changes from0.9171(v7) to1.0241(candidate); relative
strong-scaling efficiency rises57.32%→64.01%. This is largely a worse160 denominator:
160 becomes64.37ms slower while256 becomes5.33ms slower. It is not an absolute
256-task improvement. Candidate ordinary means actually rise564.709→570.450ms,
speedup0.9899, relative efficiency61.87%. Do not declare the256 scaling issue solved.
No>=640 or multi-node qualification was performed.

Mean summed rank peak RSS, GiB; sums are not simultaneous job peaks:

| Tasks | Main | 266 | V7 | Stage279 | Stage279 vs V7 |
|---|---:|---:|---:|---:|---:|
| 83 | 11.07 | 21.48 | 24.92 | 26.09 | +4.71% |
| 160 | 19.04 | 31.58 | 43.71 | 46.32 | +5.97% |
| 256 | — | 43.05 | 61.21 | 65.10 | +6.35% |

Reported transform context count stays1, and plans are built12 times for12
trajectory states, reused across the three provisional transitions within each
step. These observations do not establish long-run memory stability. Pooling the
new buffers reduced the initial candidate's allocation, but actual RSS still
increases relative to v7. Retained final-transform workspace and metadata overlap
remain costs.

## Attribution: replacement construction overwhelms retirement

Separate coarse runs, ordinary mean self-time on each run's same slowest whole-step
rank. These figures are not bare medians and cannot be subtracted from them as a
speedup forecast. Profile median perturbation is1.44–5.04% for candidate and
2.96–5.02% for v7. One profile per count is diagnostic, not replicated attribution.

| Tasks | New plan85 | New exchange90 | New total85–90 | Old envelope before | Old envelope after | Net envelope increase |
|---|---:|---:|---:|---:|---:|---:|
| 83 | 128.82 | 47.34 | 187.66 | 201.49 | 85.54 | 71.70 |
| 160 | 117.23 | 39.73 | 164.02 | 140.16 | 54.58 | 78.45 |
| 256 | 68.13 | 45.60 | 119.16 | 159.68 | 51.35 | 10.83 |

Old envelope: inverse14–21, publication/boundary32,36,38,40, scaffold42,
wavelet43 and compression44. It includes necessary work and final-stage callers;
its decrease is observed, not a claim that every included cost was retired.
At160 the new164.02ms replaces an85.57ms reduction, leaving+78.45ms in this
profiled envelope. At83 the corresponding excess is71.70ms; at256 it is10.83ms.
RHS exchange also grows21.99/10.90ms at83/160. MPI waits include producer lateness;
this does not isolate network transport or justify generic communication batching.

Audited validation preparation segments independently identify where to inspect.
On the rank slowest in preparation for each ordinary event, physical-owner home
resolution averages28.83/26.91/12.72ms at83/160/256, and route/adapter binding,
ordering and scratch setup78.01/71.58/35.54ms. These two segments are roughly
79%/81%/63% of audited preparation. Geometry compilation itself is only
0.48/0.22/0.04ms. Audit timings are a different run and rank-selection statistic;
do not splice them into the bare or coarse totals. These observations support
investigating repeated owner/route construction, not blaming numerical arithmetic.

Every plan's topology and mask epoch changes across successive states. The
implementation already reuses within each step; this is not three complete builds
per step. Reuse across steps requires a demonstrated stable substructure contract.
Ignoring epochs, matching counts or caching old routes is not justified.

## Decision and next gate

Reject Stage279 as a performance repair in its present form. The cluster validates
the local slowdown warning and shows a larger penalty at83/160. Ownership progress
is real, but it does not satisfy both goals. Keep accepted266 and the exact v7
control; retain this source/oracle evidence for a bounded redesign, without
promoting it or stacking another ownership migration.

Before another campaign, address the measured home-resolution/route-construction
and duplicate-storage costs. Show how stable geometry/ownership information can
be bound without reconstructing equivalent structures, while preserving changing
topology, masks and exact boundary versions. Quantify replacement costs first.
Recovering v7 requires about49.20/64.37ms net ordinary-median savings at83/160;
legacy parity needs143.07/210.14ms. These are measured target gaps, not promised
savings from any profile region. Even the entire new plan budget is smaller than
the main gap, so removing it cannot be claimed as full legacy closure.
Finer units/load balance and final lifecycle remain separate architecture problems.
Do not substitute larger J6/J8 workloads or a generic shared-memory campaign.
No repeat campaign is needed merely to establish the83/160 regression.

## Reproducibility note

The trusted local reviewer verifies every raw member and scientific gate. Its
Python3.14 run stopped only at a final derived percentage:13.699784523169267 versus
cluster13.699784523169246, a difference2.13e-14 percentage points(12 ulps after
subtraction). Independent Python3.11 recomputation matches the cluster comparison
exactly. Full verification then passes with Python3.11 sequential summation applied
only to diagnostic.py under3.14; existing raw/science/eight-ulp checks are unchanged.
The compatibility wrapper, exact raw archive, regenerated metrics and trusted-kit
identity check are retained alongside this report. This rounding issue has no
scientific or performance-decision consequence.
