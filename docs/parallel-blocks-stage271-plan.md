# Stage271: one fused migration proof and a production scaling baseline

Stage270 identified substantial Domain-bound shared producer/capture work and
scalar restriction, plus 64.341 ms of migration per ordinary J5 step. Inside
migration, three consecutive proof reductions cost 13.137 ms of mean-rank time,
checksum scans 7.403 ms, and payload MPI 8.163 ms. The candidate profiler passed
calibration; main's failed. These observations support a candidate-local change,
not an exact decomposition of the 209.311 ms legacy gap.

## One implementation change

In `exchange_block_migration_payloads`, replace three consecutive scalar
`MPI_Allreduce` calls with one three-element `MPI_INTEGER8`/`MPI_SUM` reduction:

1. Global outgoing byte checksum.
2. Global incoming byte checksum.
3. Global outgoing byte count.

These sums already use the same communicator and become available at the same
point after both checksum scans. They have no intervening numerical dependency.
The new code assigns the three results back to the same variables and retains
the global checksum equality failure, validated-state transition and diagnostics.
The reduction includes ranks with zero traffic. The existing representable
integer-sum semantics are preserved; no floating-point operation is regrouped.

Only `src/parallel_block_mpi.f90` differs from accepted Stage266 `ff740744`.
There is no selector, extra Domain omission, change to payload MPI, checksum
scan, migration ordering, physics, boundary refresh, or installation lifecycle.
The code removes exactly two collective invocations per payload transaction.

This is a bounded experiment, not a claim to solve the main scaling limitation.
The largest architectural target remains the Domain-owned shared producer and
its capture/restriction chain. The migration proof is small enough to isolate
without mixing numerical ownership changes into the same comparison. Its
13.137 ms observed cost includes waiting for late ranks; fusion does not remove
that waiting or checksum scans. No credible speedup percentage can be promised
before measuring the complete uninstrumented solver.

## Why the larger producer remains next

`block_tendency_compatibility_ml` still updates Domain boundaries and surface
pressure, computes diffusion/divergence inputs, calls `basic_operators`, executes
mass/velocity work on geometry owners, and captures products for the block
consumer. Native tendency destinations do not make all those producers
independent of Domain ownership. A replacement must preserve the layer/level
ordering, current boundary epochs, primitive/diffusion inputs and exact numerical
policy, then retire the corresponding Domain reads and transfers as one complete
producer/consumer unit. Stage267's reader audit remains applicable.

The scaling baseline below establishes whether useful block work at 256 ranks
is still constrained by these paths. Interpret the Stage270 component costs as
priorities; do not label the whole 194 ms shared-producer region avoidable.
After this one intervention is accepted or rejected, use the scaling evidence
to scope that substantial ownership change. Do not stack more small omissions
while postponing Domain-owner bottlenecks.

## Campaign and acceptance

One exclusive nk/bb13 job makes four serial builds with the same 23 physics
objects: accepted Stage266, legacy main `6a245db7` with the established checkpoint
I/O hardening, optimized candidate, and checked candidate. No profiler is added.
The kit carries sealed sources, fixtures and helpers; bbserv needs no old kit or
result directory. Physical placement uses one hardware thread per physical core,
with the established two-socket policy extended to 160 and 256 ranks.

Twenty solver runs use the unchanged J5 fixed-work trajectory:

| Purpose | Ranks | Versions/order |
|---|---:|---|
| Correctness | 16 | main, Stage266, candidate, checked candidate, candidate oracle |
| Correctness above Domain count | 256 | candidate oracle |
| Timing | 83 | main, Stage266, candidate, candidate, Stage266 |
| Timing | 160 | main, Stage266, candidate, candidate, Stage266 |
| Timing above Domain count | 256 | Stage266, candidate, candidate, Stage266 |

All runs retain science, states, checkpoint, 62 VTK members, restart/remap,
storage and placement checks. Candidate results must match Stage266 exactly.
Each oracle retains 48 exact inverse-return proofs and 84 installations. The
256-rank oracle includes 96 ranks without geometric Domains; its timing is
excluded from production scaling. Legacy main is run only at supported counts.

Two independent candidate/Stage266 samples at each count allow paired A-B-B-A
comparisons and expose mixed signs for a potentially small change. Main uses one
anchor per supported count; those main gaps are screens rather than replicated
precision estimates. Keep the nine ordinary steps (exclude 1, 8 and 9). Report
both paired gaps, run spreads, complete-step means of run medians, speedup and
relative parallel efficiency against 83- and 160-rank anchors. Steps within one
trajectory are not independent repetitions. No automatic repeat or promotion.

Promote only after the correctness gates pass and complete-solver evidence
supports benefit without a repeatable regression at another tested rank count.
Mixed or noise-sized effects remain inconclusive; a performance regression must
be resolved before any further change is stacked. This experiment does not
establish efficient scaling merely by running at 256 ranks.

The strategic goal is **at least 640 ranks for J5's 160 geometric Domains**.
That needs a separate multi-node campaign with suitable placement and enough
independently balanced blocks. The current four-subtree split caps capacity at
640 and permits only one block per rank at that endpoint. Further subdivision
and removing Domain-owner work remain part of the goal. This single-node screen
is an intermediate measurement, not completion of the refactoring.

## Validation and runtime records

Local validation uses serial builds and three sequential four-rank J4 replays:
optimized, checked and independent oracle. Reference binaries are reused only
after their source, binary and shared physics hashes match the sealed controls.
An additional harness compiles the actual baseline/candidate payload routines
and uses four real MPI ranks at O0/check and O2. Uneven signed-byte traffic with
an idle rank, all-zero traffic and deliberate receive corruption verify the
unchanged payloads, all three sums, three-versus-one collective count, and
retained checksum rejection. The harness does not measure performance.

Workflow tests cover all four builds and twenty runs, independent review,
collection after early failures, duplicate submission protection, purged job
handling, exact source scope, scaling arithmetic, accounting, and 256 unique
physical-core mappings. Fresh extraction is checked before delivery.

Stage270 took 15:38 for five builds/fifteen runs. Budget about 20–25 minutes
after allocation for Stage271's four builds/twenty runs; the longer 256-rank
oracle and filesystem load can change this. Queue wait is additional. The
three-hour Slurm limit is a ceiling. Runner UTC start/end and monotonic elapsed
are saved for complete or failed runs; collection captures raw Slurm accounting
and parsed Start/End/Elapsed, explicitly recording unavailable accounting.

All local validation above has now passed: both new serial builds, three exact
four-rank solver replays (62 VTK members each), all 48 oracle inverse proofs and
84 retained oracle installations, twelve actual-routine MPI cases at O0/O2,
and nine workflow tests. Fresh extraction verified all 45 sealed files, all
source members, unchanged controls, the candidate against the checkout, shell
syntax, helper imports, disabled audit selectors and no macOS xattr headers;
the nine delivered workflow tests passed again.

The first local MPI launch was blocked by the sandbox's socket restriction
before the solver started. The completed builds were reused for the successful
sequential replays with local socket access; no numerical failure occurred.
No compiler or second MPI replay ran alongside an active replay.

See the [no-edit run commands](../test/parallel_block_profile/collectives_271/README.md)
and [sealed preparation record](parallel-blocks-stage271-preparation.json) for
the archive identity and completed evidence. Cluster performance and final
acceptance remain pending.
