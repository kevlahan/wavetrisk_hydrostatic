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
branched directly from Stage266. Its source is the job55504-qualified v7 candidate;
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
do not infer transport dominance from MPI waits or automatically reuse rejected
partition/cache proposals. The current scorecard contains evidence and next gates.

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
