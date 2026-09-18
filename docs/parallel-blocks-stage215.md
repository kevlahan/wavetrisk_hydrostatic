# Stage215 — compile the remaining generic Domain halo routes

Promotion: the user authorized committing the validated candidate after review of
job54714. The exact validated `comm_mpi.f90` is now applied on
`codex/performance-investigation`. Main is unchanged. Use Stage202 on, Stage203
off and FP contraction off for the accepted performance baseline. The campaign
kit and historical evidence below remain preserved as originally recorded.

Stage215 is an integrated solver candidate for one bounded performance screen.
It targets generic `comm_mpi_mod` field boundary exchanges, not the already
compiled native block/inverse plans. Baseline is the accepted Stage202 path,
Stage203 off, FP contraction off. Nothing is promoted by preparing or running
this kit; main is untouched.

## Why this candidate, and why stop quickly

Stages187,193 and202 produced useful gains by amortizing horizontal preparation
across fields/layers and bypassing staging. The accepted generic Domain exchange
still repeats destination scans, level lookup/filtering and element-wise buffer
append for every field on every call. Receive counting repeats the same scan;
finish scans again to install. Those operations survive independently of the
native block plans reviewed in Stage214.

This is a specific exception to Stage214's recommendation against another
generic-cache campaign: the earlier audit did not implement or screen this
remaining Domain path. It does **not** replace the proposed continuous block
producer/consumer port, nor demonstrate that port's economics. If this candidate
fails the early gates, stop this tactic rather than following it with another
observer or successive tiny variations.

The older Stage196 detailed evidence puts block Domain boundary preparation and
copies at75.226ms per ordinary step (8.217ms in WT plus67.009ms outside WT), with
132.250ms posting/waiting. That is a different historical window and includes
work this candidate retains. Stage207's reconciled RK Domain-boundary envelope is
165.551ms, versus304.082ms in legacy. These are **overlapping scopes from different
campaigns, never additive**, and rank-mean observations, not critical-path bounds.
A100ms net win would require substantial elapsed-time benefit beyond local route
preparation. Stage193 showed that preparation reduction can also reduce waiting;
that is a reason to screen, not a prediction. There is no evidence yet that this
candidate meets100ms. The path is shared with legacy; it is not an explanation of
block-only excess, and might also help legacy if later applied there.

## Work removed and retained

The candidate changes only `src/comm_mpi.f90`, delivered as a reproducible patch.
It compiles signed integer send/receive lists once for each requested level
interval, per actual communication topology. Plans contain nonempty Domain-pair
groups in the original peer/source/destination order. Fields and layers use those
lists in the original order. Cached counts size a retained send buffer and allow
direct indexed packing without per-value append. Receive installation uses the
same signed lists. Production receive-buffer zeroing is removed because posted
receives cover exactly its consumed prefix.

All routes are discarded when `comm_communication_mpi` completes an initialization
or refinement route rebuild, and when the communication module is initialized.
No field values, masks, owners, pointers into field arrays or numerical results
are cached. MPI calls, tags, wire values/counts, receive-before-send posting,
remote snapshots before local copies, negative-edge signs, field validity flags,
local-copy order and even all-level local copies for a single-level request remain
unchanged. This saves preparation/copy overhead; it does not remove MPI messages
or migrate the primitive producer to block owners.

The checked build re-compiles from the original lists on **every** plan use and
compares counts, group ordering and signed addresses. It fills receive storage
with NaNs and rejects any unfilled active element after completion. These checks
run only in validation, with bounds/FPE checks. Production timings have no added
observer or audit.

## Gates and stop rules

1. At16 tasks: baseline, checked candidate, candidate with existing numerical
   oracles, then candidate without oracles. Exact checkpoints and adaptive
   workload must agree before any screen.
2. One uninstrumented A/B pair at16. Advance only for at least3% gain. Otherwise
   report no clear gain or regression and stop; no repetitions to rescue it.
3. One A/B pair at83, with the established fixed42/41-core placement and local
   memory policy. Require at least3% and at least90ms gain to justify confirmation.
   Smaller gains stop as below target. This is a screening rule, not an estimate
   of statistical significance.
4. Fresh A–B–B–A at83. Keep the100ms target; both adjacent-pair gains must meet it.
   Repeat one batch only for disagreement or a result within10ms of the target.
   Two disagreeing batches are inconclusive. Screening pairs are not confirmation.

Each run contributes one ordinary-step median. Exclude the first timed step after
initial load and every restart, checkpoint writes, and logged field/data saving.
Keep ordinary vertical remaps. Timesteps within a run are not independent
repetitions. Use the existing frozen J5/J7 Z30 fixture and its restart/save
schedule for matching; no elapsed startup/output times enter the estimate.
`dt_write` is not changed because the solver clamps it and the existing fixture
also exercises reload correctness. Output exclusions are enforced independently.

Full source and climate callback identities, helper hashes, binary identities,
placement, memory/RSS and decisions are retained. Collection is atomic and verified;
an existing good archive can simply be retransferred. A purged Slurm job is
handled by querying the user's active-job list, not `squeue -j`.

See `parallel-blocks-stage215-local-results.json` for completed local numerical
validation and `parallel-blocks-stage215-preparation-results.json` for contract,
helper and packaging checks. Local J4 validation is not a J5 or performance claim.
The isolated idle-rank communication test does not establish complete solver
startup above160 tasks.

## Local reproduction

`candidate_boundary_215.py TREE` applies to a verified, uninstrumented Stage202
reference tree. Build with `STAGE202=on STAGE203=off FP_CONTRACT=off`; use
`PARAM=param_J4` for the local fixture. The checked build also uses `DEBUG=check`
and `-DSTAGE215_AUDIT`. Separate build/bin directories are required.

`run_boundary_test_215.py` accepts `--baseline`, `--candidate` and `--checked`
compiled object directories plus a fresh `--out`; it links the actual modules and
runs four MPI processes. `local_operations_208.py` runs each complete J4/J6 fixture
window, with `--oracles` only for the independent-oracle run. Final acceptance is
`validate_local_215.py --parent-kit ... --reference ... --legacy ... --optimized ...
--checked-off ... --oracles-on ... --halo-contract ... --out ...`. All reference
inventories and tested binaries must match; local test paths and hashes are
preserved in the validation JSON. No local timing is used for acceptance.
