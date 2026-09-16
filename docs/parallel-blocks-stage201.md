# Stage 201: instructions versus memory stalls in complete dynamics paths

Objective: identify source-level mechanisms capable of removing at least about
0.10 s/step from the remaining 0.511 s block/legacy gap. This is a diagnostic
campaign, not another small optimization. Stage 197 candidate and legacy
numerics remain unchanged. Stage 200 supplies the verified 83-task, two-socket,
32-cache-group placement on exclusive bb13. No one-socket or rank-count sweep.

## Measurements and decisions

Compare complete, disjoint dynamics wavelets and remaining-dynamics phases;
retain adaptation and other phases for conservation. The same common phase
boundaries are used on both sides. Observer mode disables fine operation/MPI
clocks, work counters and access census. Count snapshots occur only at common
phase/transaction boundaries. Counters run in the solver main thread and report
user-space execution, including MPI progress/spinning. They exclude kernel and
other threads; they are not a solver-only arithmetic instruction count.

Two four-event groups, collected in separate repeated passes:

| Group | Events | Question |
|---|---|---|
| core | cycles, instructions, `ex_no_retire.not_complete`, `ex_no_retire.load_not_complete` | Is excess execution mainly more instructions, higher CPI, or substantial load-related retirement waiting? |
| memory | cycles, instructions, L2 data demand misses, near-or-far DRAM/MMIO demand fills | Is the slower execution accompanied by more demand misses/fills per step and per instruction? |

The AMD raw configurations are `0x02d6`, `0xa2d6`, `0x0864`, and `0x4843`.
The last combines event 0x43 masks 0x08 and 0x40. These are specific to the
validated family 25/model 160 CPU; the helper rejects other models. The Linux
perf CPU map selects the Zen 4 table for this family/model. Event provenance:
[Linux CPU map](https://raw.githubusercontent.com/torvalds/linux/master/tools/perf/pmu-events/arch/x86/mapfile.csv),
[core events](https://raw.githubusercontent.com/torvalds/linux/master/tools/perf/pmu-events/arch/x86/amdzen4/core.json),
[cache events](https://raw.githubusercontent.com/torvalds/linux/master/tools/perf/pmu-events/arch/x86/amdzen4/cache.json).

Load-retirement waiting is not equivalent to DRAM stalls. Demand fill counts
are not bandwidth, do not count all traffic, and include the documented MMIO
case. No conversion to physical DRAM bytes or assertion of bandwidth saturation.
Raw event acceptance is not a semantic hardware self-test; preserve CPU and
kernel identity, counts, scheduling and repeat agreement for review.

Each counter read reports group time enabled/running. Keep raw counts; never
silently scale heavy multiplexing. Invalid/zero major-phase counts or less than
95% scheduling stop the campaign for review. Per-rank/per-step intervals include
absolute CLOCK_MONOTONIC timestamps, thread CPU time, and event counts. They
partition the common phase ledger. Exclude restart steps when aggregating.
[Linux group-read semantics](https://man7.org/linux/man-pages/man2/perf_event_open.2.html).

Report absolute counts, instructions ratio, CPI/IPC and load-wait/miss evidence
for both implementations. A symmetric algebraic decomposition splits cycle
excess into instruction-count and CPI components; it is explicitly **not** a
causal decomposition of wall time or recoverable savings. MPI spinning can
increase both instructions and cycles. Do not normalize solely by the old
temperature-node proxy, which omitted other representation work.

## Source attribution

Sample cycles and instructions in separate event streams using fixed periods
of 20 million events and 8 KiB DWARF call stacks. No frame-pointer rebuild or
change to production optimization flags. `perf record --clockid mono` matches
solver timestamps. Export symbols/stacks on bb13 while executables and shared
libraries exist; retain raw perf data too. Sampling is limited to the main
thread (`--no-inherit`) consistently with the PMU counts.
[perf recording options](https://man7.org/linux/man-pages/man1/perf-record.1.html).

Use the same seven ranks for both implementations: **19, 46, 52, 63, 65, 79, 82**.
These are the union of low/median/high outside-MPI dynamics costs on each
socket, for each implementation, in Stage 200. The kit freezes the selection
and source-results hash. This is a diagnostic cohort, not a representative
node average. All 83 ranks still execute and all receive phase counters in PMU
passes. Samples are timestamp-matched to ordinary-step phases. Output retains
per-rank leaf symbols, DSOs and inclusive ancestry, with sparse/unresolved/lost
sample flags. Inclusive call stacks overlap and must not be summed.

Interpretation requires identifying MPI/UCX progress separately from solver
work; then locate instruction or cycle concentration in scalar record expansion,
lookup/address resolution, Domain/block publication/import, packing/copy/zero,
allocation/traversal, or numerical kernels. Preserve unknown symbols. Fixed
sample periods allow within-event sample shares, but short scopes and sampled
ranks do not justify exact per-function wall seconds. Do not claim a small
routine explains a whole phase because it appears on an inclusive stack.

## Bounded campaign and gates

Five fresh builds: bare candidate, checked candidate, candidate observer, bare
legacy, legacy observer. Optimized sources/flags match Stage 200: J5, RK4,
`-O2 -march=native -ffp-contract=off` built on bb13; common physics objects.
The observers differ only in the shared clock/counter module. Production
numerics remain frozen; no commit or source optimization is part of this stage.

Before any solver run, test both PMU groups on all 83 prescribed cores at once,
and verify perf sampling/unwinding with both events on a short standalone
program. Save capabilities and failures. If perf is missing or access is denied,
stop early, collect the same review archive, and resolve access with cluster
support. Do not use sudo, change kernel policy, or silently substitute software
counters. No network access is required on the cluster.

20 solver runs:

- Two checked candidate runs, oracle off/on.
- One warmup per implementation and two bare repeats per implementation.
- One counter-disabled timeline control per implementation.
- Two core and two memory passes per implementation, reversing side order.
- One sampled pass per implementation with counters disabled.

Use the unchanged checkpoint 3 fixture, time_end=0.315, dt_write=0.01,
checkpoint 4 save/reload and post-restart remap. All runs retain numerical,
workload and placement gates. Same-implementation comparisons are exact;
legacy comparison uses the existing fixed screen, not bitwise equality.

Compare every observer's printed timing against bare and its parent phases
against the disabled-counter control. Flag >3% printed overhead, or a parent
change exceeding both 2% and 1 ms. Flags require review, not a global correction.
Counters and sampling execute separately. Sampled-rank resource usage describes
the perf wrapper and descendants and is not compared to bare solver-only RSS.
Do not use perf-export or full-process elapsed time as timestep performance.

The outcome should rank complete paths by measured excess and explain whether
extra source work, load-related stalls, MPI spinning, or a mixture dominates.
Only then propose a source intervention with a credible substantial recovery
estimate. If evidence is ambiguous or sampling is biased, repeat only the
specific diagnostic needed; do not automatically repeat the entire campaign.
