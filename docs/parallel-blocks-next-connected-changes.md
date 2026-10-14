# Next substantial changes: connected RHS production cutover

Current completion status: the connected production selector and final local
qualification are complete; the sealed bbserv kit is ready for user submission.
See [the production result and run commands](parallel-blocks-rhs-production-278-results.md)
for current evidence and remaining goal gaps. The design/reassessment text below
is the chronological pre-implementation record; its draft-readiness statements
are superseded by that result. Accepted source remains266; no cluster promotion.


Prepared 2026-10-09. This is the concrete continuation of the Stage277 handoff
and the revised objectives. The accepted source remains Stage266
`ff740744608731a05d53771c32a3055309f3a586`; Stage278 numerical development is
inactive and has not earned ownership or performance credit.

Implementation checkout:
`/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic`.
Single authoritative scorecard:
`/Users/kevlahan/.codex/worktrees/26d8/wavetrisk_hydrostatic/docs/parallel-blocks-scorecard.md`.
[Field/reader/cost contract](parallel-blocks-connected-rhs-278.md).
Preserve the original checkout's Stage273 partition candidate and uncommitted
review documents. Do not fold subdivision, partition tuning or Stage277 caching
into this numerical comparison.

## Quantitative revision, 2026-10-09

The [quantitative performance case](parallel-blocks-rhs-performance-case-278.md)
now controls the implementation/layout decision and the campaign entry gate.
It verifies998 archived raw observer files; no new solver timings are claimed.
The83-rank candidate cost pool is49.196ms of mixed capture/remainder work, not
193.987ms of freely removable producer arithmetic. Only an identified removed
fraction earns savings; charge replacement halos/preparation and producer cost
changes once. The scenario with50% removal, unchanged producer cost and10ms new
cost saves14.60ms in a mean-rank screen; a10% producer cost increase makes it a
4.80ms regression. These assumptions are not a forecast.

Revised production design: direct writes to existing resident consumer storage,
local RHS diffusion/physics, shared epoch-valid geometry, patch-contiguous fused
loops and phase-bounded scratch. Retain the separate Simple split-physics adapter.
Compile complete route bytes/peers/phases and identify exact removed calls before
sealing; a reused Domain RHS-physics residual stream has no current performance
justification. Preserve all scientific ordering and the full connected RHS/RK
boundary. No standalone module or producer benchmark is the next campaign.

The raw160/256 block work maxima are unchanged in five of nine ordinary steps;
current ownership alone does not justify1.6x scaling. Keep partition/subdivision
separate. The campaign below is conditional on this revised implementation and
cost gate; no positive whole-step forecast or launch readiness has been established.

## Decision from the initial evidence

The512-case numerical development gates establish the draft's tested arithmetic
only. They do not establish distributed producer placement, epoch correctness,
full-step science, retired work, or speed. The next deliverable therefore is an
**independently selectable connected production RHS transaction**, not another
numerical helper collection, producer benchmark or geometry-owner cache.

Two goals remain mandatory and are evaluated separately:

| Goal | Current scorecard evidence | Next candidate must establish |
|---|---|---|
| Independent ownership | Mass/velocity shared production and RK compatibility still geometry-owned. Stage278 has no production caller and retires nothing. | Authoritative RHS and actual RK consumers on final block owners, including zero-Domain ranks, with matched interior Domain computation and transports absent. |
| Performance and scaling | Accepted266 is39.80% behind main at83 and54.72% at160 in job55497;160-to-256 speedup0.9333, relative efficiency58.33%. | Bare complete-step comparison including every new cost, no repeatable regression to accepted266, and explicit remaining gap to main. Subsequent matched bb160/320/640 evidence is needed for>=4D acceptance. |

Baseline non-regression is an intermediate qualification, not achievement of
legacy parity. Correct zero-Domain execution is an ownership check, not proof
of efficient scaling. No efficiency threshold is invented here.

## Reassessment: what to retain and what to change

Reviewed evidence: [Stage275 input budgets](/Users/kevlahan/.codex/worktrees/26d8/wavetrisk_hydrostatic/docs/parallel-blocks-stage275-results.md),
[Stage276 resident execution](/Users/kevlahan/.codex/worktrees/26d8/wavetrisk_hydrostatic/docs/parallel-blocks-stage276-results.md),
[Stage277 correctness and matched timing](/Users/kevlahan/.codex/worktrees/26d8/wavetrisk_hydrostatic/docs/parallel-blocks-stage277-results.md),
[Stage278 boundary-test scope](parallel-blocks-rhs-boundary-278-results.md).
The run matrix was checked against `primitive_277/common.py` and `run.py` in the
handoff checkout, and its measured runtime against the campaign ledger there.

Reassessed 2026-10-09 against the actual accepted call path and reviewed
Stages274–277. **Do not launch the current Stage278 draft.** Adding its module to
Makefile did not connect it to production: only the unit driver imports it.
Recent preparation made numerical progress but stopped short of the requested
ownership change. More standalone arithmetic additions would continue that drift.

| Evidence or draft component | Disposition and implementation consequence |
|---|---|
| Stage275 minimal Domain-fed velocity inputs cost3.27–3.46x existing tendency export at83,3.12–3.38x at160 and3.21–3.48x at256 | Reject velocity-only relocation fed by Domain primitives. These are sampled per-layer budgets, not measured whole-step traffic or speedups. Flux/PV account for52.84–52.91% of the sampled256-rank input budget: place their shared producer with resident consumers. |
| Stage276 exact resident primitive inputs/outputs, including sampled ghosts; all96 zero-Domain ranks execute at256 | Reuse address/residency evidence and fused arithmetic. Extend validation to every actual call/layer and the complete consumer chain. Its reference Alltoall stream remains validation-only. |
| Stage277 all-call assignment checks pass but zero-Domain ranks do zero primitive work; timing changes are small and mixed | Retain arithmetic/oracle coverage. Do not inherit the geometry-owner cache automatically or run another standalone277 campaign. |
| Stage278 direct/SW/pentagon numerical tests | Retain as development regression tests. Synthetic slot maps and manually supplied physics residuals do not prove real traversal, final-writer delivery, halos or full science. |
| Draft per-event geometry, gather arrays and hot-loop address checks | Revise around actual block patch storage. The primitive descriptor embeds49 doubles of geometry (392 bytes with8-byte dp), plus addresses and flags, including SW storage on ordinary nodes. Account for replication and workspace size; share immutable metrics by valid topology/layout epoch. Validate addresses when installing plans, retain transaction epoch guards and checked-build bounds checks. Do not claim a measured bottleneck yet. |
| Draft whole-column pressure helper | Preserve accepted summation and vertical arithmetic, but design storage and schedule around actual pressure/Exner consumers. Raw column Exner and fine-to-coarse restricted Exner must not alias prematurely. Avoid assuming three full node-by-layer output arrays are the necessary production layout. |
| Draft pointwise diffusion helpers | Insufficient: implement the actual scalar/divu/rotational hierarchies, halos, pentagon viscosity and post_vort ordering. Synthetic nonzero physics inputs cannot qualify that pipeline. |
| Draft rhs_rk helper | Keep as an arithmetic test only. Reuse the existing saved-base accumulator/trial/commit path in parallel_block_mpi; do not create a second resident RK state machine. |

Concrete accepted-source anchors are `time_integr::prepare_block_multistage_tendency`,
`multi_level::block_tendency_compatibility_ml`, `basic_operators`,
`prepare_native_velocity_programs`, and
`parallel_block_mpi::begin_block_domain_multistage_candidate_stage`.
The latter already assembles resident stages from a retained timestep-start base.
`RK_sub_step_compatibility` still performs geometry-owner mass/velocity closure;
`retain_block_native_multistage_candidate` still has publication expectations.
Change those readers and assertions together, including the final-stage alias.

The implementation decision is one supported-mode RHS/RK selector spanning the
shared producers, hierarchy dependencies and actual resident consumers. Keep the
accepted selector as the control. Do not stage multiple partial production
activations unless each independently closes its readers and retires paired work.

Four integration questions must be answered in source before packaging, using
the real compiler and local replay rather than a fresh broad cluster survey:

1. Can every active mass restriction, mixed direct/restrict velocity action and
   boundary alias be mapped to an owned or phase-complete halo slot? Compile
   actual dependencies; Stage274/275 samples do not prove universal locality.
2. Does every shared-vorticity/pentagon output preserve the last legacy writer
   and each intermediate read? Encode ordered phases and explicit deliveries;
   accepting only the final field can miss an incorrectly consumed intermediate.
3. Are topology/ownership/layout, actual mask predicates, stage state and halo
   epochs distinct? Accepted velocity plans already detect mask-only changes
   without a generation change; retain that property for mass and velocity.
4. Which Domain materializations remain required by inverse/scaffold, split
   Simple physics and lifecycle readers? Record their exact field/phase/byte
   contract. Remove duplicate RK arithmetic, but do not promise removal of an
   interior publication still required by a retained reader. Count its cost.

## Evidence-based benefit and its limit

Job55497 gives the latest same-allocation baseline: at83 ranks accepted266 is
760.845ms versus main544.225ms; at160 it is551.163ms versus356.235ms. Thus the
observed gaps are216.620ms and194.929ms. Closing them would require approximately
28.47% and35.37% reductions in accepted266 time respectively. These are required
reductions to parity, **not predicted gains**. The often-quoted39.80%/54.72%
figures instead use main as the denominator.

The proposed change can distribute currently geometry-owned RHS/RK work and
remove its tendency-export/compatibility costs. That is the causal opportunity.
Stage270's roughly194ms shared-producer self time is neither wholly removed nor
a portable speedup budget. Stage276's sampled gather/reference timings are not
production-step measurements. New halo, preparation, retained-adapter and memory
costs may offset savings; no defensible speedup percentage exists yet.

For matched run-level medians let T266, Tc and Tmain denote accepted, candidate
and main times. Report saved time T266-Tc, candidate change100*(Tc/T266-1),
remaining legacy gap Tc-Tmain, and fraction of the matched legacy gap closed
(T266-Tc)/(T266-Tmain), when the denominator is positive. For160-to256 report
speedup T160/T256 and relative efficiency speedup/1.6. Latest accepted values
are0.9333 and58.33%; useful work beyond D must improve measured outcomes, not
just increase the count of active ranks. No numeric acceptance tolerance or
final efficiency target is invented.

## Reviewable change sets, one connected activation

These changes may be developed separately, but no intermediate producer-only
path is a candidate for promotion. Enable the complete transaction once its
readers and replacement costs are connected. If a genuinely independently
activatable smaller slice is found, document its complete reader closure and
retirement before qualifying it; do not assume each change set is such a slice.

### 1. Compile actual block workspaces and routes

Extend the numerical draft with per-local-block RHS workspaces and a compiler
in the MPI adapter. Reuse `Block_Data` resident scalar/vector and mean families,
`Block_Stencil_Address`, boundary/ghost storage descriptors and established
catalogue/writeback identities. Geometry owners may compile geometry metadata
on topology change; ordinary RHS arithmetic and state operands must execute
on destination block owners.

Remap every `Mass_Restriction` address (including its fine flux, partial-node
and coarse-node operands), every ordered `Velocity_Level_Program` action and
`Velocity_Gradient_Stencil` into owned/halo slots. Retain direct numerical loops;
resolve storage class/producer once at compilation, not at each arithmetic read.
Install metrics/stencils only when their topology/ownership/layout epoch changes.

Explicitly preserve mask-only rebuilds: `prepare_native_velocity_programs`
compares actual TRSK/ADJZONE predicates even if its generation is unchanged;
`prepare_native_mass` separately compares TRSK/RESTRCT predicates. A single
address-generation test is insufficient. Keep mask and metric lifetimes separate
so unchanged geometry is not resent on a mask-only adaptation.

Before installing a plan, verify each nonlocal operand has a producing phase,
all active addresses are mapped, and target writes are unambiguous. Confirm
source and mass-restriction dependencies against actual compiled topology,
not only Stage275 sampled plans. Reject an unresolved dependency before entering
production; do not create an implicit Domain-owner fallback or a dynamic
per-value dependency interpreter. A real same-level cross-owner dependency
requires a documented bounded phase split preserving action order.

### 2. Complete shared producers and hierarchy dependencies

The mapped `step1::comp_SW` and all four `post_step1` pentagon arithmetic cases
are now implemented and independently unit-checked. Connect them, positive
Domain-boundary mapping and rotational `post_vort` to the actual fused producer schedule. Preserve the original node traversal and the final
shared-vorticity writer, including `post_step1` and `post_vort`. Tests must use
accepted legacy expressions for these cases, with unequal/nonuniform inputs,
orientation changes, shared aliases and all relevant masks.

Connect raw pressure/Exner and Bernoulli production to ordered fine-to-coarse
injection and FINAL restricted gradient operands. Keep raw column Exner separate
from restricted Exner; do not reuse the existing resident hydrostatic cache's
alternative surface-pressure summation order as an exact RHS replacement.

Implement the full scalar/divu/rotational diffusion sweeps and their node/edge
halos. The draft's pointwise diffusion functions alone are insufficient.
Preserve state versus Laplacian inputs, perturbation versus state-plus-mean,
pentagon viscosity reduction, and the order2 vorticity overwrite. All prepare,
restriction, exchange, allocation and physics costs belong to the transaction.

### 3. Execute the authoritative mass/velocity RHS

At the actual tendency entry, finish the stage's state/mean halo and bind its
state transaction. Evaluate columns and physical layers, with explicit fine-to-
coarse edge/node/source phases from the field contract. Use final restricted
mass flux for both mass divergence and velocity source. Mark a layer ready only
after the final B/Exner gradient and required communications complete.

Install results directly into the existing local tendency/accumulator families.
Keep the native temperature path connected to the same producer epochs. Remove
the production substitution of completed Domain mass RHS into the scalar physics
channel (`publish_native_mass_rhs`) for moved interiors. Remove the velocity
RHS export/remainder path for moved interiors. No full primitive return is added.

### 4. Connect real RK and retention consumers, retire paired work

Change `time_integr::prepare_block_multistage_tendency` to use the connected
block RHS entry in the new selector. The accepted selector remains available
as the exact control. The oracle continues to run `trend_ml` independently in
validation only; bare production must not run duplicate reference arithmetic.

Reuse `begin_block_domain_multistage_candidate_stage`'s saved timestep-base,
accumulator and reversible trial semantics. `rhs_rk` is numerical development,
not a substitute for wiring these existing production registers. Respect
`production_multistage_captured_tendency_stage`, stage count, checkpoint and
trial lifetimes. The final RK call currently aliases q with its destination;
resident stage assembly must still use the retained timestep-start base once.

Replace the interior dependency on `RK_sub_step_compatibility` and update
`retain_block_native_multistage_candidate` publication expectations together.
Stage266's temperature-only publication count assertion cannot simply be
removed. Separate the necessary boundary/scaffold/inverse adapter from the
newly unnecessary interior materialization, and assert the resulting contract.

Activation must prove these paired removals for moved interiors:

| Displaced computation | Displaced data movement / materialization |
|---|---|
| RHS surface pressure, pressure integration, shared primitive assignments | RHS-only interior state/primitive preparation for Domain readers |
| B/Exner injection and scalar/divu/rotational diffusion | Their geometry-owner numerical boundary/operand streams, replaced by block phase halos |
| Native mass restriction/divergence on geometry owners | Completed mass-RHS substitution/publication and its former mass boundary work |
| Native velocity source/restriction/gradient on geometry owners | Native tendency export and compatibility remainder capture |
| Geometry-owner interior mass/velocity RK closure | Duplicate interior stage copy/materialization for the replaced readers |

Every retained scaffold or adapter must name its fields, reader, stage frequency,
bytes and remaining Domain arithmetic. Simple split physics is a separate
`main::apply_domain_physics_consumers` transaction after RK: retaining it must
be explicit in the scorecard and full-step timings. Its later migration,
remaining inverse and lifecycle readers are necessary for complete independent
timestep ownership; RHS acceptance alone is not project completion.

### 5. Qualify before promotion or another cutover

Use serial local builds and<=4 MPI ranks, with no concurrent build/replay.
Retain optimized, checked and independent-oracle fixture runs against sealed
accepted266 and supported main controls. Cover exact accepted-family science,
main-relative tolerances, mask-only adaptation, topology changes, reload/remap,
shared boundary epochs and final-stage aliasing. Unit tests remain development
evidence, never a replacement for these complete-step gates.

## Proposed next bbserv campaign: connected production qualification

**Initial proposal status (superseded): no sealed kit and no cluster launch.** Base the candidate
on accepted266 `ff740744608731a05d53771c32a3055309f3a586`, with an explicit
source overlay. Main remains `6a245db7f1632b0205c2bebe95e5f0cd535c2582` plus the
established checkpoint I/O hardening. Exclude Stage273 and Stage277 caching.
Use the immutable established J5 restart fixture and its recorded hashes.
Scope: restarted compressible non-mode-split RK4 climate with Simple physics,
without NCAR topography/SSO. Other modes, cold starts and ranks with no blocks
remain outside this qualification.

### Entry gate: complete locally before delivery

Require a full optimized and checked climate build, independent-oracle full-step
replays on at most4 local MPI ranks, and actual activation of the new selector.
The current module-only build and synthetic tests do not meet this gate.
Exercise topology and mask-only changes, reload/remap, final-stage aliasing,
boundary/scaffold epochs and supported diffusion branches. Record exactly which
branches the immutable J5 campaign exercises; branch-specific local fixtures
must cover any supported branch absent there. No claims for unexercised modes.

Validation must show the displaced Domain numerical entry points and superseded
exports are unreachable for migrated work, using counters plus validation-only
poison/absence checks at actual readers. Distinguish independent oracle calls
from production work. Unsupported or unresolved routes must fail closed before
output mutation, without an implicit Domain fallback. A dry-run/fresh extraction
must verify source/fixture/binary identities and review contracts; update old
publication-count assertions to the new reader contract rather than disabling
them. Only then prepare proven launch/status/collection scripts for user launch.

### Leg A: nk/bb13, one exclusive256-CPU node, mem=0

Keep four serial builds: optimized main, optimized266, optimized connected
candidate, checked connected candidate. Probe and verify existing placement at
16/83/160/256. No concurrent builds or solver runs. Execute this fixed order:

| Phase | Ranks | Runs in order | Count |
|---|---:|---|---:|
| Validation |16| main;266; candidate with full contract checks; checked candidate with checks; candidate with checks and independent oracles |5|
| Validation |83| candidate with full contract checks |1|
| Validation |160| candidate with full contract checks |1|
| Validation |256| candidate with full contract checks and independent oracles |1|
| Bare timing |83| main;266; candidate |3|
| Bare timing |160| main;266; candidate |3|
| Bare timing |256|266; candidate |2|

Total: **eight validation runs plus eight bare runs**. This reuses the Stage277
schedule but replaces the obsolete assignment-only contract with the connected
ownership/retirement contract. All runs retain12 steps,2 reloads,2 remaps and62
VTK members, exact accepted-family science and established main-relative gates.
Keep the ordinary-step exclusions0/7/8 unchanged. Also report complete trajectory
elapsed, rebuild/adaptation costs and memory so excluded exceptional steps cannot
hide a new regression. All communications and retained adapters remain inside
bare solver execution. Disable duplicate arithmetic, reference packets and detailed
validation instrumentation for bare timing; verify that they are absent.

In validation, record actual Domain/block ownership and authoritative work by
phase for every rank, all48 RHS calls and30 layers, including all96 zero-Domain
ranks at256 with relevant owned work. Verify consumers use those results, not
shadow copies. Reconcile expected work with the compiled catalogue, rather than
hard-coding Stage277 node counts. Check no missing/duplicate work, final-writer
semantics and real RK stage installation. Record plans/rebuilds, workspace and
peak memory, halo/adapter bytes and messages, replaced exports and retained
publications; collect bounded phase costs only in instrumented runs. These
measurements answer the candidate's replacement-cost question, not a new audit.

Fail science/ownership immediately and preserve raw failure evidence. Once those
gates pass, evaluate bare full-step costs at all counts. If a small or unexpected
timing difference could change acceptance, allow **at most one additional reversed
control group per affected count** in a separately reviewed follow-up schedule: candidate/266/main at
83 or160; candidate/266 at256. Maximum eight additional bare runs,24 solver runs
overall. Timesteps are not independent repeats. Report each run median; if this
bounded repeat remains inconclusive, hold promotion and name the unresolved
measurement instead of launching another generic survey. A repeatable regression
requires recovery within this slice before stacking another ownership change.

Planning estimate for the base schedule: **25–40minutes elapsed**, excluding queue
and collection. This is a provisional engineering allowance above Stage277's
measured18:45 for the same four-build/sixteen-run schedule; the new full-chain
oracle and route checks have no measured cost yet. Use a60-minute wall limit,
not a predicted runtime. Before sealing, revise the estimate from measured build
and local validation costs; retain the60-minute cap for the fixed initial schedule; budget any reviewed repeat separately.
Record actual Slurm accounting and the exact schedule in the runtime ledger.

### Leg B: conditional inter-node readiness, same candidate source

Run only after Leg A ownership/science passes and its timing decision is resolved.
This tests a specific remaining question: whether the new phase/halo route works
and remains competitive across nodes. It is not a subdivision experiment or a
claim of efficient4D execution. Do not include it automatically in the first kit.

Request four exclusive bb nodes, exclude bb04, mem=0. Capture their actual
CPU/NUMA topology and verify compatible hardware, compiler flags and bindings
before building/running; do not reuse a bb13-native binary blindly. Distribute
160 ranks as40/node and256 as64/node, using the same allocated four nodes and
placement policy for each control/candidate pair. Build three optimized binaries
serially (main,266,candidate). Run two candidate validations (160 with contract
checks;256 with checks and oracles), then five bare runs: at160 main/266/candidate;
at256266/candidate. Total **seven runs, three builds**. Apply the same science,
retirement and failure gates as Leg A. No automatic repeat extension for this leg;
a timing ambiguity holds its conclusion pending a specifically justified follow-up.

Do not infer inter-node cost by subtracting nk times from bb times, and do not
join the two hardware families into a scaling curve. There is no measured bb
runtime estimate for this candidate: a provisional60-minute wall cap is a resource
limit, not a forecast. Reassess feasibility before sealing, record actual elapsed,
and retain an incomplete result honestly if the limit is reached.

## Scorecard evaluation for each returned campaign

Record exact source/overlay/fixture/binary identities and reviewed completeness.
For each selected production field, record executing owner, zero-Domain-rank
coverage, consumer and epoch evidence, retired computation and removed transport,
new bytes/messages/storage/rebuild costs, and any retained adapter. If those
entries are absent, ownership qualification is incomplete even with exact science.

For performance, record run-level ordinary-step medians under established
exclusions, matched main and accepted gaps, independent repeats when justified,
node/placement/topology, Slurm elapsed and failed runs. Compute
`speedup=T(reference ranks)/T(target ranks)` and
`relative efficiency=speedup/(target ranks/reference ranks)` only on matched
hardware. Do not treat timesteps as repetitions or join nk and bb scaling data.

Classify the result as development-only, correctness-qualified pending timing,
production accepted, regression requiring recovery, or unpromoted/inconclusive.
Acceptance requires BOTH ownership/science and the performance gate; report
legacy gap and>=4D qualification separately even after an intermediate acceptance.
Maintain the original scorecard in place; do not create a competing current copy.

## Following milestones, kept separate from this source comparison

Leg B above is the bounded multi-node readiness check. It does not close the
final scaling gate or justify mixing a partition change into this candidate.

Develop finer subdivision/balance separately from numerical ownership; current
at-most4D blocks provide no balancing slack at4D ranks. After independent
qualification, measure J5 on the same ten bb nodes at160/320/640 ranks
(16/32/64 ranks per node, D=160), with
longer adaptive workload and memory/halo evidence.640 ranks is the minimum
fourfold target;256 ranks remains intermediate. Agree a useful-efficiency
criterion before final acceptance. Explicit MPI shared memory remains deferred.
