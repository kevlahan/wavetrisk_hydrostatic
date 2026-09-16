# Stage 198: account for the complete legacy/block time-step difference

Stage 198 is a diagnostic experiment. Its numerical candidate is the validated
Stage 197 code, compared directly with the pinned legacy Domain implementation.
Stage 196 remains an untouched historical control. There is no new numerical
optimization, ownership change, arithmetic flag change or relaxed correctness gate.

The objective is the [performance-gap diagnosis](parallel-blocks-performance-gap-diagnosis.md):
explain the main sources of the roughly 0.69 s/step excess before choosing another
optimization. The previous 83% figure located two expensive regions; it did not
establish a removable cause. This campaign measures those regions and the remaining
adaptation/reconstruction and outer-phase costs together.

## What the experiment records

The observer assigns each interval to exactly one operation and context. It covers
the full step: setup, dynamics wavelets, the remaining tendency/dynamics work,
physics, adaptation wavelets, remaining adaptation/reconstruction, remap and output.
There are 100 named operation IDs plus an explicit unclassified ID. Nested scopes
contribute only their exclusive time to the ledger; a parent is not added to its
children. Restart steps are recorded and validated but excluded from ordinary-step
comparisons.

| Matched task or question | Measurements and limits |
|---|---|
| Complete dynamics wavelet transaction | Transform/compression, scalar and vector inverse kernels, local preparation, packing, installation, ownership transfers and Domain publication; both complete transaction envelopes are retained. |
| Generic and specialized boundary work | Shared Domain scanning, routing, allocation/zeroing, local copies, posting, waiting and installation, alongside block inverse, scalar, mass and ghost paths. Include the generic-path savings, not only extra block costs. |
| Complete tendency/RK work | Shared producers and basic operators, native velocity/mass preparation, scalar capture/address/geometry setup, restriction, divergence, RK retention and compatibility closure. Legacy fused producer envelopes remain fused. |
| Adaptation/reconstruction | Masks, compression, Domain inverse reconstruction, catalog and source-block construction, migration/install, retained reconstruction and lifecycle completion. These intervals are no longer outside the detailed clock. |
| Storage representation | Separate sampled binary measures scalar range reads, value reads, shared-address resolution and writes. Logical slots, calls and selected work counts accompany the timings. Inclusive samples are never added to exclusive operation times. |
| Rank imbalance and completion | Per-rank records, rank-mean distribution and the operation ledger from the actual slowest full-step rank, selected once per step. Independent operation maxima are never summed. |

`operation-instrumentation.json` supplies exact IDs, source-hook inventories and
counter definitions. New call-site envelopes include conditional branch checks;
their visit counts are not necessarily counts of subroutine execution when a
condition is false. Existing fused kernel/traversal labels are not arithmetic-only
measurements. A shared name does not establish identical mathematical work.

Forty-two work counters distinguish selected evaluations, initialization, scanning,
packing, copies and logical storage requests. They are not a complete census of
every equation or hardware load. Hot legacy leaf counters use guarded integer
accumulation with reporting outside the leaf loops. Inverse dependency records
retain rank-local sequence, field/level/kind, requests, values and completion
intervals. They cannot establish globally synchronized arrival times or prove that
MPI wait is independently removable.

## Calibration instead of an expensive speedup campaign

The default is 27 short executions on the same exclusive 83-rank nk allocation:

| Executions | Purpose |
|---:|---|
| 2 | Checked candidate, production and oracle modes |
| 3 | Stage 196, Stage 197 and legacy warmups/reference checkpoints |
| 6 | Two bare runs of each numerical implementation |
| 4 | Two minimal common-clock runs of Stage 197 and legacy |
| 4 | Two exclusive operation-timing runs of Stage 197 and legacy |
| 4 | Two timing-plus-work-counter runs of Stage 197 and legacy |
| 2 | Separate storage-probe binary with probes disabled |
| 2 | Same binary with systematic samples enabled |

The second round reverses implementation/control order. Storage probes use stride
257 and different offsets, 0 and 127. Their estimates are systematic samples, not
random samples with statistical confidence intervals. The sampler is compiled out
of the ordinary observer binary. Its disabled control reveals binary-layout/hook
effects; enabling samples adds counting and timing costs as well. Neither effect
is subtracted as a uniform correction.

Common parent envelopes are compared across minimal, timers and work passes.
Changes exceeding 1 ms or 2% of the reference parent (whichever is larger) are
flagged for review; whole-step observer changes above 2% are also flagged. These
are diagnostic thresholds, not statistical error bounds. If instrumentation
distorts a region by as much as the effect being explained, that attribution is
insufficient and needs a less intrusive follow-up measurement.

## Output and acceptance

`ranking.csv` contains all 909 signed context/operation rows, largest block excess
first, including negative differences and unclassified intervals. Each row gives
legacy/block seconds, excess, call-site visits and the slow-rank-path difference.
`ranking.md` shows leading rows, matched parent totals and calibration. `results.json`
retains every run, gate, work count, rank distribution and storage sample summary.
The review archive retains the underlying per-rank common, operation, dependency
and sample records as well as logs and build/input/placement provenance.

The sum of signed operation differences must equal the difference of the same
instrumented rank means. The bare median gap, additive rank-mean gap and actual
slow-rank-path gap are reported separately. The report preserves negative offsets
and explicitly labels unclassified time. It does not rescale the ledger to the
bare median gap or treat a low timer residual as proof of complete causal knowledge.

Every run retains the existing checkpoint/restart, workload, remap and numerical
gates: exact equality against its own reference and the fixed cross-implementation
legacy screen. The latter is not bitwise legacy equivalence. Missing ranks,
operations, full-step coverage, dependency records, sample coverage, nonfinite
values and conservation failures are rejected. The schedule exercises 65 numerical
and workload comparison gates in addition to per-run validation.

Seven fresh J5 builds use identical physics objects and contraction off. Build
identities, final observer sources and executable hashes are verified before
submission. Provenance travels in a file, avoiding the prior command-line size
failure. Source instrumentation is confined to disposable trees; stripping the
new hooks restores the preceding observer source bytes. The frozen numerical
candidate patch is byte-identical to Stage 197.

The report deliberately leaves `optimization_ready` false. Review must still map
the significant measured intervals to matched producers/consumers, distinguish
required from avoidable work, and estimate recoverable completion time without
double-counting another rank's wait. Cache, bandwidth and NUMA causes are not
established by logical-byte counts or memory footprint alone. Significant fused
or observer-distorted regions remain explicitly unresolved and require further
diagnosis rather than a speculative optimization.

## Local clock correction

An initial local work-counter run stopped on a backward clock reading. Inspection
of the installed GNU Fortran 15.2 runtime showed its 64-bit SYSTEM_CLOCK calling
`clock_gettime` with clock ID 0, which the local macOS SDK defines as CLOCK_REALTIME.
The final observer instead links a small identical C helper using CLOCK_MONOTONIC
on both platforms, relative to a process-local origin. It fails on clock errors;
negative durations are not clamped or silently discarded. A million-read stress
check and fresh MPI runs validate the replacement. No synchronization between
rank clocks is implied. The C helper is compiled through the MPI C wrapper so it
uses the same configured compiler family and header search paths as the MPI
Fortran build. Initial runs using the previous clock are superseded.

## Running and collecting

Use the [cluster instructions](parallel-blocks-cluster-198.md). They build on bb13,
save submission state, monitor the correct log, and collect one review tarball for
transfer to `murier.math.mcmaster.ca:~/Desktop/`.

Local test results and the final frozen kit identity are recorded separately in
`parallel-blocks-stage198-local-results.json` and `parallel-blocks-stage198-kit.json`.
Local runs validate correctness and measurement coverage; they do not estimate
the 83-rank cluster performance gap.
