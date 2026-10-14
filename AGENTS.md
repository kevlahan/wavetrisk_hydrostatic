# Block refactoring: current objectives and operating rules

User-authorized objectives reset: 2026-10-09, following the project review and
Stage277 results. Read `docs/parallel-blocks-objectives.md` and maintain
`docs/parallel-blocks-scorecard.md` after every reviewed campaign or accepted
production change. Those documents supersede conflicting historical next-step
recommendations. `docs/parallel-blocks-measured-restart.md` retains the history;
`docs/parallel-blocks-handoff-277.md` supplies the continuation context.

## Outcome and next substantial deliverable

- Complete independent block ownership of the supported production timestep,
  removing non-essential Domain numerical pathways and redundant transfers.
- Achieve legacy-comparable or better performance at P<=D and useful efficient
  execution at P>=4D. J5 D=160 requires >=640 ranks; 256 is an intermediate test.
  No numerical efficiency threshold has been agreed; do not invent one or claim
  that correctness on zero-Domain ranks establishes efficient scaling.
- Next deliver a connected block-owned production RHS chain and its real RK
  consumers, including necessary shared pressure/Bernoulli/physics/diffusion and
  restriction dependencies. Retire displaced Domain computation and transfers
  in the same qualified slice. Do not substitute another geometry-owner
  optimization, diagnostic prototype, or bulk primitive-return path.
- Record producer/reader ownership, topology/mask/state/halo epochs, replacement
  costs and retired work before changing placement. Use the existing evidence;
  new diagnostics need a named unresolved implementation decision.
- Preserve direct fused numerical loops, compact epoch-valid stencils and bounded
  communication. Learn from the Stage251–255 regression: do not recreate the
  per-value/dependency-depth owner interpreter and defer optimization.
- Geometry identifiers, justified external adapters and independent oracles may
  remain. Counting Domain references or renaming storage is not a completion metric.

## Accepted baseline and latest decision

Accepted solver: Stage266 `ff740744608731a05d53771c32a3055309f3a586`.
Legacy main control: `6a245db7f1632b0205c2bebe95e5f0cd535c2582`, with established
checkpoint I/O hardening. This worktree follows codex/stage278-structural-repair,
branched directly from Stage266. Its history records the job55504-qualified v7 candidate followed by the complete
job55508-reviewed Stage280 source, including Stage279 and the two-file binding/storage repair;
Stage273/277 partition/cache changes are excluded. Do not confuse a candidate
commit with promotion of the accepted solver. See docs/parallel-blocks-branches.md.
Build further candidates from the qualified source and record their deltas.

Stage277 job55497 is independently reviewed: 4,323 files, sixteen successful
runs, all assignment/science gates pass, elapsed18:45. It remains unpromoted.
Candidate timing changes versus266 are -0.415%/+1.094%/+0.452% at83/160/256;
single runs establish no clear benefit or repeatable small regression. At256 all
96 zero-Domain ranks perform zero primitive nodes. Retain its numerical kernels
and assignment evidence for reuse, not its cache automatically. No further
standalone Stage277 tuning/repeat campaign is planned. Stage266 remains accepted.

Stages274/275 support local compact scheduling and sampled input residency;
Domain-fed velocity inputs still cost over3x the existing tendency export under
the observed budgets. Stage276 validates sampled block-owned primitive execution,
including ghosts, but retired no production work. Stage277 validates all-call
primitive arithmetic on Domain owners. These are distinct evidence scopes.
Stage267 canonical publication/boundary dependencies remain until readers move.

## Current Stage278 continuation

Job55504 verifies all48 runs and the exact v7 source. Ordinary-step time improves
19.9%/15.5%/10.5% versus its predecessor at83/160/256 and exceeds266 performance.
Main parity, whole-timestep ownership, memory parity and efficient scaling remain
open. Retain v7 as the development base. Inspect load imbalance, insufficient
block granularity, fixed per-rank overhead and inherited Domain dependencies;
do not infer transport dominance from MPI waits or automatically reuse unpromoted
partition/cache proposals. The current scorecard contains evidence and next gates. The ordinary-step follow-up
in docs/parallel-blocks-stage278-balance-assessment.md finds mean per-step work
max/mean1.354→1.989 at160→256 and a conditional5.49% upper bound on peak primitive
work reduction from reassignment of unchanged blocks at256. Distinguish that
operation-count bound from runtime gains; target granularity/scaffold work and
retained dependencies separately. That assessment itself made no production change; the subsequent Stage279 implementation is recorded below.

## Selected next repair specification

See docs/parallel-blocks-resident-transform-repair-plan.md for the concrete
producer/reader and retirement contract. The selected next ownership slice closes
resident RK through vector/scalar wavelets, compression and inverse into the next
RHS, initially at the three provisional RK4 transitions. The three-provisional-transition slice is implemented and locally science-qualified.
See docs/parallel-blocks-stage279-resident-transform.md: nine four-rank replays,
24 MPI exchange cases and32 helper tests pass; the sealed49-run bbserv kit awaits
user submission. Local bare point estimates are +2.11% median/+4.53% mean versus
v7 under large variation: no speedup or performance acceptance is established.
Job55505 failed during GCC13 checked compilation before any solver run. V2 fixes
one allocatable selection expression; full isolated remote GCC13 checked/optimized
builds and three fresh local optimized/checked/oracle replays pass, with34 helper
tests. The v2 kit awaits user submission; v1 and failure evidence are preserved.
See docs/parallel-blocks-stage279-gcc13-fix.md.
Job55506 completed six validation solvers but stopped at83 in an incorrect
per-rank kernel-class checker. V3 changes helpers only: positive total per owner,
exact global class census across counts, and shared current revision identity.
Six cluster runs and three local audits pass reanalysis;40 helper tests pass.
Sources/flags remain exact v2. No performance measurements or160/256 validation
yet. See docs/parallel-blocks-stage279-checker-fix.md; v3 awaits user submission.
Resolve the cluster result before stacking another cutover. Evaluate finer,
cost-aware units separately; include internal/coarse and scaffold work. Keep J5
as the deciding track. Possible later J6/J8 tests do not replace the original
above-Domain target. The current scorecard records remaining gates and budgets.

## Latest campaign decision — job55507 supersedes preparation status above

Independent review verifies19,118 files, five builds and all49 runs; elapsed43:52.
Stage279 passes science, provisional transform ownership/paired retirement and
96 zero-Domain participation at256. It FAILS the performance gate. Ordinary
median times683.278/569.310/555.903ms at83/160/256 versus v7
634.078/504.941/550.577ms: +7.76%/+12.75%/+0.97%. All83/160 candidate medians are
slower than all v7 medians;256 overlaps. Ordinary means regress13.70%/14.62%/6.00%.
Memory rises4.71%/5.97%/6.35% versus v7. Better median160→256 ratio mainly reflects
a worse160 denominator; candidate means still worsen with256. No scaling closure.

Reject promotion; keep266 accepted and v7 the exact development control. Preserve
Stage279 source/evidence but do not stack more migration or repeat the regression
campaign. The measured next investigation is expensive per-step owner-home and
route/adapter construction plus duplicate storage: new preparation128.83/117.23/
68.13ms in coarse critical-rank ordinary means. Every step has a new recorded
mask/topology epoch; existing within-step reuse is not the issue. Require stable
substructure identities and a quantitative net cost case before another repair
campaign. No transport-dominance inference, generic shared-memory shortcut or
J6/J8 substitution. Main gap is now143.07/210.14ms at83/160; full ownership,
finer units/load balance, memory parity and matched>=640 remain required.
See docs/parallel-blocks-stage279-job55507-results.md and current scorecard.
This review changed no solver source, baseline, branch or remote state.

## Latest campaign decision — job55508 supersedes pending/preparation notes

Independent review verifies22,280 files, six builds and60 successful runs in53:44.
Stage280 preserves science/checked/oracle, provisional numerical ownership and
paired retirement; all96 zero-Domain ranks do transform/RHS work at256. All
source/kit identities match. Ordinary mean run medians improve3.67%/4.25%/3.74%
versus exact v7 at83/160/256; recovery versus279 is13.52%/14.91%/4.71%. Every
repair median is below every v7 median. Remaining main gaps62.78/126.73ms.
Ordinary mean gains are2.12%/2.47%/0.35%;256 means overlap with retained outliers
and non-remap mean regresses0.52%. Do not claim repeatable256 throughput gain.

Retain280 as a measured construction repair candidate for continued repair of
this slice; preserve v7/279 controls, keep266 accepted and main untouched. No
accepted-source promotion, additional migration, repeat or new cluster campaign.
Coarse transform plan falls132.60→34.12,116.37→27.13,67.95→30.60ms. Numerical
traffic, messages, route calls and kernel census are unchanged versus279: this
supports the construction repair, not a generic network/batching conclusion.
RSS still exceeds v7 by3.75%/5.28%/5.60%; whole-memory parity fails.

256 remains8.95% slower than160 by medians,11.68% by ordinary means. Median
relative efficiency57.37%, mean55.96%; no scaling closure or>=640 qualification.
Offline existing evidence confirms numerical imbalance1.354→1.989, max work
falls only8.07%, maximum arithmetic rises35.03→37.00ms. Fixed-bundle assignment
alone permits at most5.49% peak primitive reduction at256, not a speedup forecast.
RHS preparation grows58.29→70.72ms; exchange130.94→151.44ms includes lateness.
Next quantify structural RHS construction/adapter repair, finer numerical and
scaffold units and retained memory using existing source/evidence, before another
implementation/campaign. Preserve dual goals and minimum decision-changing run
policy; no generic shared memory, unsafe mask reuse, assignment-only scaling
promise or J6/J8 substitution. Final lifecycle, memory/main parity remain open.
See docs/parallel-blocks-stage280-job55508-results.md and current scorecard.
Review changed only evidence/docs/offline analyzer; no solver or remote actions.

## Current implementation — Stage280, cluster performance pending

The user authorized repair of job55507's demonstrated construction/storage costs.
Stage280 preserves the provisional numerical ownership slice.16 ordinary routes
use current unique canonical providers directly, eliminating128 redundant
preparation collectives per epoch. Outer and fixed-input writer proofs remain.
Local publication reads canonical rows; unique remote imports share storage.
All epoch/state checks, one-context lifetime and within-step reuse remain.
No partition or additional numerical migration; two solver files change.

11 full local solver replays and48 optimized/checked MPI cases pass. Exact J5
GCC13/MVAPICH checked/optimized builds pass in an isolated remote compiler copy.
Local4-rank J4/J6 matched medians improve4.79% vs v7 and5.38% vs279(two runs each),
with nonoverlapping repair/control ranges. Named extra-transform numerical/route
allocation falls29.64%; total RSS remains unqualified. Remote numerical payload,
route invocations and global work census are unchanged. Local coarse RHS waiting
grows while the whole coarse benefit is small; retain that limitation. No cluster
speedup/scaling/main-parity claim, no promotion and no further migration before
83/160/256 validation. Exact controls and source snapshots are retained.
See docs/parallel-blocks-stage280-prebound-repair.md and .stage280-validation.
The sealed rk280-prebound-kit-v1 is prepared in .stage280-validation for user
transfer/submission: six serial builds,11 validations,49 measurements,60 runs.
42 current helper tests and fresh extraction of all five source families pass.
Tracked helpers: test/parallel_block_profile/prebound_campaign_280. Planning
60-80 minutes within90. No new kit upload/submission. Use scorecard ordinary
time/RSS/ownership/scaling gates; maintain v7/279/266/main controls.
Only the two source files were uploaded for compiler checks; no remote solver
or Slurm campaign was run. Legacy main and accepted266 remain untouched.

## Campaign size constraint — user instruction 2026-10-14

Use the minimum runs needed to decide performance improvements, scaling
improvements, repairs or architecture choices. Before sealing a campaign, name
its unresolved decisions and justify each configuration, control, repetition
and diagnostic against those decisions. Reuse applicable source-bound evidence;
retain validation needed for the changed source and above-Domain ownership.
Do not automatically multiply every control, rank count and profile by the same
repetition count. Extra repetitions require a specific decision-changing
uncertainty; extra profiles require a named attribution or architecture question.
Predeclare any bounded follow-up criteria rather than inflate the initial kit.
If a small effect is below variation, decide whether it is worth pursuing before
requesting more runs. Report build/run counts and estimated total cost up front.

Job55508 completed its sealed60-run schedule and is independently reviewed
above. This constraint governs future campaigns; the historical60-run matrix
must not become a default. No additional campaign submission authorized.

## Qualification and scope

- Gate each independently activatable connected production change immediately
  against accepted source and main where supported, including above-Domain
  timing. Resolve repeatable regressions before stacking another migration.
  Repeat only when uncertainty changes the decision. Current baseline parity
  is an intermediate gate, not closure of the existing legacy gap.
- Preserve optimized/checked/independent-oracle tests, accepted-family exact
  science, established main-relative tolerance, topology/masks, reload/remap,
  boundary/scaffold epochs and zero-Domain-rank coverage. Keep duplicate
  arithmetic/reference packets out of production timings.
- Current root-child subdivision supplies at most4D blocks; finer granularity
  and useful load balance are necessary for efficient >=4D execution. Test
  partition/subdivision separately from numerical ownership changes.
- Current P>D support is restarted compressible, non-mode-split RK4 climate with
  Simple physics, without NCAR topography/SSO. Other configurations, cold starts
  and ranks with no blocks are not generally qualified. Do not overstate scope.
- Explicit MPI shared memory is deferred. It does not distribute Domain-owner
  arithmetic or fix inter-node scaling. Any later test must distinguish within-
  NUMA, cross-NUMA and inter-node paths and preserve local first-touch placement.

## Operational constraints

- Local machine: four usable cores; serial builds and at most four MPI ranks,
  without concurrent build/replay work. Local timings are not cluster evidence.
- nk contains only bb13 (256 effective CPUs); no bb04 exclusion needed. bb has
  bb01–bb12 (64 effective CPUs each); exclude bb04 on bb/all. Established jobs
  use exclusive allocation and --mem=0; small explicit memory requests have failed.
- Verified bb13 topology: two EPYC9754 sockets,128 physical cores/socket, one
  online thread/core, two NUMA nodes and32 L3 groups. bb-node topology and page
  residency remain unverified; capture actual allocation topology before binding.
  Compare on identical nodes and placement; do not combine nk/bb timings as one
  scaling curve. Ten64-CPU bb nodes can supply640 ranks, subject to allocation.
- Record real Slurm start/end/elapsed and build/run schedule in
  `docs/parallel-blocks-campaign-runtimes.json`. Separate queue/collection time.
  Stage270 took15:38, Stage27616:11, Stage27718:45 on their respective schedules;
  new ownership/multi-node campaigns need their own estimates.
- Preserve sealed source/fixture identities and original failed-job evidence.
  Never use unrestricted stage-number replacement: it previously corrupted a
  grid checksum. Verify packaged bytes and fresh extraction before delivery.
- User commands must be fenced and copyable, without # comment lines or editable
  placeholders. Use the proven launch/status/collection structure.
