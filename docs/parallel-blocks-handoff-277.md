# Handoff after Stage277 and objectives reset

Prepared 2026-10-09. The user requested revised objectives that eliminate drift,
substantial progress toward complete block ownership, legacy-comparable timing
and scaling beyond the Domain limit. They explicitly requested maintaining the
new scorecard and then authorized creating a handoff chat in this project.
The new local project chat is `01a12152-11da-7bf0-9580-8381731e6dd7`, created
2026-10-09 and confirmed active. Its opening instructions explicitly identify
this existing checkout, the uncommitted documents and the accepted source;
the project's default checkout must not be mistaken for that source.

## Read first

1. `AGENTS.md` — current constraints and latest decision.
2. `docs/parallel-blocks-objectives.md` — authorized objectives and acceptance.
3. `docs/parallel-blocks-scorecard.md` — maintain this current scorecard in place.
4. `docs/parallel-blocks-stage277-results.md` and `.json` — verified latest return.
5. `docs/parallel-blocks-project-review-2026-10-09.md` — historical project assessment;
   its Stage277-pending statements reflect its creation time, not current status.

Consult detailed history only when needed. Do not repeat the broad project audit
or start a broad profiling campaign. The next substantive work is the connected
block-owned production implementation; no Stage278 kit has been prepared.

## Workspace and source identity

Workspace: `/Users/kevlahan/.codex/worktrees/26d8/wavetrisk_hydrostatic`.
Branch: `codex/domain-migration-measured`.
Accepted solver: `ff740744608731a05d53771c32a3055309f3a586` (Stage266).
Main control: `6a245db7f1632b0205c2bebe95e5f0cd535c2582` plus the established
checkpoint I/O hardening. Latest pre-review commit is `b5f03c8e` (Stage277 kit).

Working `src/parallel_block_mpi.f90` includes an unpromoted Stage273 partition
candidate. Do not stack it into the next solver implicitly or destroy it.
Stages274–277 use sealed Stage266 overlays. Recent review/objectives/scorecard
documentation may be uncommitted; inspect status and preserve it. The pre-existing
untracked `.stage261-validation/`, `.stage262-validation/`, `.stage263-validation/`
directories must be preserved. No solver source was changed by this objectives
reset and no new cluster campaign was launched.

## Latest reviewed result

`/Users/kevlahan/Desktop/rk277-review-55497.tar.gz`:
SHA256 `aa268d0d5e239ea4bd945ade6397e7921bcc323d99d4824cb96dcd9a5cb1fbb8`.
Reviewer: `test/parallel_block_profile/primitive_277/review.py`.
Saved raw local review: `.stage277-validation/review-55497.json`.
4,323 file hashes and all sixteen runs verified; delivered kit manifests match.
Slurm COMPLETED,18:45, four serial builds/eight validations/eight single timings.

Stage277 is **not promoted**: -0.415% at83,+1.094% at160,+0.452% at256 versus266;
small single-run changes, no established benefit. All96 zero-Domain ranks at256
still execute zero primitive production nodes. Preserve its fused arithmetic
and 19-value assignment checks as reuse material, not the Domain-owner cache as
an accepted optimization. No standalone repeat/tuning campaign is planned.

## Known implementation scope, not another research question

`src/multi_level.f90::block_tendency_compatibility_ml` still prepares shared
Domain inputs and `native_mass/native_velocity(size(grid))`. Source/restriction/
gradient arithmetic is direct, but geometry-owned. `basic_operators` and
`src/ops.f90::step1` supply physics, pressure/Bernoulli/Exner and shared primitives.
`src/time_integr.f90::RK_sub_step_compatibility` materializes the retained mass/
velocity stage. `parallel_block_mpi.f90` exports native velocity tendency to final
owners. The accepted Stage203 inverse workspace is also geometry-owned.

Stages274/275 show moving velocity alone with Domain-produced inputs requires
over3x the current tendency-export volume under observed budgets. Stage276
validates resident primitive operands on block owners at sampled phases; its
reference Alltoall packets are validation only. Stage277 validates production
assignment arithmetic across all calls/layers but changes no execution owner.

Complete the field/reader/epoch/new-cost/retired-work contract and implement a
connected production RHS ownership slice, including dependencies that otherwise
force bulk primitives back to Domain owners. Preserve fusion, exact restriction
ordering, column/physics sequencing, final shared-vorticity writes and boundary
epochs. Retire displaced work together. Do not recreate the old per-value,
dependency-depth owner interpreter or accept a slower intermediate cutover.

Then remove newly unnecessary representation transitions and tackle the remaining
inverse/lifecycle ownership. Granularity/balance is a separate necessary workstream:
current maximum4D blocks gives no balancing slack at4D ranks. Final J5 target is
>=640 ranks on multiple bb nodes, not256 on bb13. No numerical efficiency threshold
has been agreed. Explicit MPI shared memory remains deferred.

## Suggested opening instruction for the new chat

Continue the block refactor in this workspace using AGENTS.md and
docs/parallel-blocks-handoff-277.md. Follow the revised objectives and maintain
docs/parallel-blocks-scorecard.md. Start from accepted Stage266, preserve the
separate Stage273 candidate, and prepare the next substantial connected
block-owned production RHS change with its real consumers and immediate science
and performance gates. Use the existing Stage274–277 evidence; do not substitute
another standalone geometry-owner optimization or broad diagnostic campaign.
