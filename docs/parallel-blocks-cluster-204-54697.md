# Stage204 cluster review — job 54697

The campaign passes correctness and closes the measured operation ledger, but
**does not pass all observer-calibration gates**. It narrows the investigation;
it does not yet establish a calibrated causal ranking or justify a new
optimization with a 0.10 s/step savings claim. Stage202-on remains the accepted
baseline; Stage203 remains experimental. No numerical change or commit is made
by this review.

## Validation and reproducibility

Archive SHA256: `f3648577ea66581cd1179149744c7341ddd5c2de656bf5f5abb593f7909eb2c8`.
All 6,364 inventoried files match size/hash records, with no missing or extra files.
Seven recorded source inventories and build-log hashes match the frozen kit;
recorded executable identities match all run records. The operation-site map
matches its recorded SHA256. Builds use J5, contraction off, the same 23 physics
objects, and the prescribed Stage202/203 settings.

All 20 runs complete the same 12-step J5/J7 workload, restart from checkpoint 3,
write/reload checkpoint 4, and remap after restart. All 43 archived gates pass,
including checked oracle-off/on, exact block/observer comparisons, workload
checks and the unchanged legacy tolerance screen. Legacy agreement is not bitwise.
Checkpoint/binary payloads are intentionally absent: numerical comparator reports
were re-evaluated, but actual payload comparisons and executable hashing cannot
be rerun from this archive.

Independent reconstruction reproduces raw-log analyses, all 83-rank placements
(42/41 across the two sockets), nine observer datasets, exclusive operation and
MPI-owner conservation, primary phase summaries and saved summary calculations.
There is a secondary-path qualification: maximum full-step clocks tie in legacy
operations repeat 1. Filesystem iteration order then selects different ranks;
one operation's slow-rank-path average differs by 0.028496 s between valid choices.
The archived path rows are inside independently reconstructed tied-rank bounds.
**Do not use that secondary path ranking for causality.** Rank-mean ledgers and
bare timings are unaffected. Future analysis should specify deterministic tie
handling and report ambiguity, rather than imply a unique critical rank.

No bare run swaps. Available memory remains above 1,490 GiB. Sums of per-rank
maximum RSS are approximately 21.8 GiB for blocks and 11.1 GiB for legacy (these
are not simultaneous aggregate peaks). This excludes capacity/swap pressure in
these runs, not cache or memory-bandwidth inefficiency. Whole-process major faults
include startup/I/O and are not evidence of step-local memory stalls.

## Bare performance

| Implementation | Reversed samples, s/ordinary step | Median | Excess over legacy |
|---|---|---:|---:|
| Stage202 accepted baseline | 0.925, 0.929 | 0.9270 | 0.3835 |
| Stage203 experimental | 0.890, 0.884 | 0.8870 | 0.3435 |
| Legacy Domain | 0.544, 0.543 | 0.5435 | — |

Stage203 saves 0.040 s/step (4.32%); paired savings are 0.035 and 0.045 s.
Both again miss the agreed 0.10 s/step criterion. It remains a repeatable modest
improvement, not the substantial reduction required for promotion. Ratios are
1.706x legacy for Stage202 and 1.632x for Stage203; the historical 2.25x is not
the current ratio. Two repetitions do not establish precise confidence intervals.

## Calibration limits

- Stage203 detailed repeat 2 changes dynamics-wavelet rank-mean time by
  +0.007758 s (about 3.65%) relative to its coarse control, exceeding both 1 ms
  and 2% limits.
- Legacy detailed repeat 1 has a printed step median 3.22% above bare,
  exceeding the 3% limit.
- All Stage202 controls pass; the other individual controls pass. There is only
  one coarse control per implementation, so drift and instrumentation effects
  cannot be separated conclusively.

The correct campaign decision is `observer_effects_require_review`, despite
`status: passed` for execution/correctness. Do not relax thresholds retrospectively.
The following tables describe measured locations, not recoverable savings.

## Where the instrumented difference is located

These are signed, exclusive, rank-mean wall seconds per ordinary step, averaged
across the two operation runs. They are not the printed bare medians above.
Negative differences are retained, and inclusive parent timers are never added.

| Phase | Stage202 − legacy | Stage203 − legacy | Stage203 excess inside MPI | Stage203 excess outside MPI |
|---|---:|---:|---:|---:|
| Dynamics wavelets | 0.172617 | 0.145334 | 0.086368 | 0.058966 |
| Other dynamics | 0.104655 | 0.089511 | 0.022270 | 0.067241 |
| Adaptation/reconstruction | 0.088384 | 0.088962 | 0.052664 | 0.036298 |
| Other phases | −0.003254 | −0.003308 | −0.004477 | 0.001169 |
| Total | 0.362402 | 0.320499 | 0.156825 | 0.163673 |

Unassigned absolute time is only 1.93 ms for Stage202, 1.91 ms for Stage203 and
1.04 ms for legacy. Thus nearly all instrumented elapsed time has a location.
This is accounting completeness, not a complete causal explanation of the bare
gap. MPI time includes progress, payload movement and waiting for producers;
removing a producer cost and its resulting wait must not be counted twice.

For Stage203, the net MPI excess comprises 132.7 ms in collective/blocking
operations, 21.1 ms in completion and 3.0 ms in posting. This does **not** justify
another generic message-volume optimization. The source shows these operations
also implement compatibility transfers, ownership reconstruction and validation.

## Concrete source findings and exclusions

1. **Adaptation reconstructs and validates a distributed representation.**
   `migrate_blocks` calls manifest checking and size/payload exchanges even when
   `full_validation` is false. The manifest checker performs four unconditional
   Allreduces, size exchange performs six, and payload exchange performs three;
   the latter also scans outgoing and incoming byte buffers for checksums.
   See `src/parallel_block_mpi.f90` around lines 5635, 36452, 36619 and 36835.
   Production `full_validation=false` skips a larger diagnostic suite but does
   not remove these checks. This is a concrete maintenance cost, not evidence
   that all checks are unnecessary or that their time is individually known.
   The complete migration envelope is about 58.9 ms/step for Stage203, including
   nested routing, installation, MPI and plan work; do not add its children.
   Net adaptation excess is about 89 ms, so adaptation alone cannot supply a
   defensible 100 ms saving on this measurement.

2. **Domain/block compatibility remains in dynamics.**
   Domain-to-block import and block-to-Domain publication use Alltoallv,
   packing, staging, validation and installation. Forward vector propagation
   also exchanges values and produced-coverage flags and performs a convergence
   Allreduce (`parallel_block_mpi.f90` around 10360, 9912, 13618 and 13683).
   Some synchronizations therefore establish required dependencies, not merely
   diagnostics. Removing or coalescing them needs an explicit producer/consumer
   proof and an exact numerical gate.

3. **Several scalar/RK adapter costs are visible, but not independently removable.**
   Stage203 exclusive rows include 25.8 ms scalar capture transport/install,
   22.3 ms temperature boundary contract, 21.6 ms scalar ghost local work,
   21.5 ms RK capture/staging and 7.1 ms scalar capture preparation in tendency.
   Legacy represents related work under other labels, including negative
   differences. These positive rows must not replace the signed net comparison:
   the whole tendency context has only 65.6 ms net excess (55.5 ms outside MPI).
   The nested RK capture envelope of 108.7 ms includes other listed work and
   cannot be added or advertised as 108.7 ms of avoidable cost.

4. **Final inverse installation is not a large enough standalone target.**
   The entire Stage203 final-install envelope is only 22.8 ms/step. Even its
   complete elimination is far below 100 ms. The earlier large logical value
   count was not a measurement of DRAM traffic or proof of dominant cost.
   Production logs confirm the dynamics oracle is disabled; this run does not
   support blaming a complete duplicate oracle execution for the gap.

## Recommended next decision

Do not select another isolated kernel optimization from these data. First make
one bounded, lighter confirmation of **representation lifecycle work**:

- Retain Stage202 and legacy as the primary comparison; preserve Stage203 as a
  separate experimental control, not the new accepted baseline.
- Replace simultaneous 65-label instrumentation with a few batch-level scopes
  around migration/catalog/plan rebuild, scalar/RK capture, and Domain/block
  publication/import. Separate checksum/local validation, actual payload
  exchange and dependency/coverage collectives by call site. Record generation
  changes, reuse, actual calls and values consumed at batch boundaries.
- Pair observer-on/off controls closely in reversed order. Require the existing
  calibration limits, exact block/observer gates, fixed legacy screen and
  restart/remap checks. Use deterministic tie reporting; no per-value clocks
  or barrier insertion. Repeat the affected calibration with this reduced
  observer, not another full broad profiling campaign.
- Produce a producer/consumer and lifetime map: which fields/geometry/routes are
  already valid, which are rebuilt, which consumers require materialization,
  and which validation reductions affect control flow. This is necessary before
  proposing persistent state, fewer representation transitions, or check batching.
- Only then choose a coherent prototype whose **non-overlapping, calibrated**
  measured costs support the existing ≥0.10 s/step target. Confirm savings with
  bare runs; do not equate time inside MPI or a parent envelope with removable
  time. If the evidence does not support that target, report this rather than
  lower the threshold or combine unrelated small wins.

Machine-readable evidence is in `parallel-blocks-cluster-204-54697.json`.
