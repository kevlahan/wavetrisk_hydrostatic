# Stage205: confirm representation lifetime and conversion costs

## Decision to make

Determine whether maintaining Domain and block representations causes a coherent
set of avoidable operations large enough to support a **0.10 s/ordinary-step**
prototype target. Stage204 located time but did not prove redundancy, and two
calibration checks failed. Stage205 is measurement only: numerical statements,
MPI ordering, production validation and algorithm selectors remain unchanged.
Stage202 is the accepted baseline; Stage203 is a bare experimental control.

## Focused evidence

The observer has 25 batch categories plus three selected MPI categories and one
reserved slot. It does not time individual nodes, values, field accessors or
message posts. Selected completion timers and collective timers operate only
inside the chosen enclosing batches; the remainder deliberately includes other
computation and other MPI. Thus it is not a replacement full-MPI census.

The selected paths cover catalog/source construction, migration, manifest and
size/payload validation, payload checksum scans, plan construction, scalar
geometry/capture, Domain/block transactions, RK capture/assembly, forward-vector
and inverse preparation. Legacy retains the common workflow contexts and its
Domain RK/basic-operator batch envelopes for comparison.

Every original numerical statement is byte-recoverable from observer insertion.
Conditional call hooks execute inside the original predicate, so false calls are
not counted. The preparation tool rejects missing/ambiguous anchors. The kit
records source inventories and the original call-site map for both observers.

Runtime data include:

- Existing plan-generation predicates, rebuild versus reuse requests, destruction
  of a ready versus empty writeback plan, and installed writeback generations.
- Full versus compact scalar capture, retained-cache readiness and storage rebuild
  flags; geometry seeds and existing producer-counter deltas at batch boundaries.
- Ordered field import/publication records with payload family, component or
  vector-only selector, generation and RK stage. The RK interval marker precedes
  scalar production; subsequent tendency capture and candidate updates are also
  recorded. Stage zero means no RK tendency has been entered in this step; later contexts retain the most recent stage marker and must be interpreted with the context ID.
- Scalar operator consumer entries and catalog outcomes (total blocks and
  source-to-final owner changes). Equal counts are not equal topology, and an
  ownership change here is not a comparison with the previous step's owners.
- Actual selected collective call sites, send/receive element totals when
  available, per-call durations and rank-local step-relative start times.
  Allreduce counts are element counts, not network byte traffic. Alltoallv counts
  retain the site's MPI datatype; uncounted APIs report -1 rather than zero.
- Payload bytes scanned for checksums, with a separate scan timer. Validation
  reductions can be distinguished from data transfer or vector convergence
  using the recorded call site and enclosing scope.

Events are buffered, with a fail-closed 100,000-event/rank limit, and written only
at finish. There are no added MPI barriers, reductions or timed-step log writes.
No PMU or sampling support is required. Observer-off mode writes empty lifecycle
tapes and zero detailed timers; the common coarse clocks remain enabled.

## What constitutes actionable evidence

`consumer-map-205.json` names the producers, consumers, lifetimes and falsification
conditions for four hypotheses: plan reconstruction, geometry/residual capture,
field materialization and migration validation. Runtime events are interpreted
against these frozen source contracts, not automatically labelled redundant.

Expected review output is a short table of specific operations with (1) actual
frequency and validity state, (2) local preparation, checksum and collective
costs, (3) necessary consumers and state changes, (4) calibrated non-overlapping
cost eligible for an alternative, and (5) an explicit proposed replacement or a
reason to reject that target. Repeated generation/stage alone does not prove
unchanged live values. This campaign does not add per-value comparisons or claim
an exhaustive dynamic def-use trace. A source audit of writes and consumers is
required before deleting any operation.

Falsification is useful: predominantly reused plans, one required geometry seed
per generation, or small validation cost should rule out those targets. Do not
combine unrelated small improvements to evade the 0.10 s target. If maintenance
cannot support it, investigate execution efficiency within necessary computation.

## Controlled cluster design

Six fresh J5 builds: Stage202 bare and observer, Stage203 bare and checked,
legacy bare and observer. All retain contraction off and common physics objects.
Build and run on bb13/nk; 83 ranks, pinned 42/41 across sockets, local allocation,
exclusive node. Fixture remains J5/J7, resume 3, time_end 0.315, dt_write 0.01;
checkpoint 4 reload and post-restart remap remain required.

There are **19 runs**:

1. Two checked Stage203 runs, oracle off/on.
2. Three implementation warmups.
3. Two reversed bare triples (Stage202, Stage203, legacy, then reverse).
4. Four adjacent observer pairs: Stage202 off/on, legacy off/on, legacy on/off,
   Stage202 on/off. Each implementation therefore has two paired coarse controls
   and two focused measurements. Stage203 has no new detailed observer.

All runs retain exact block/checked/observer gates, unchanged tolerance-based
legacy comparison, workload, restart/remap and placement checks. Observer outputs
must conserve common phase and whole-step clocks, MPI owner totals and collective
call-site durations. The block observer must cover all lifecycle event kinds,
four RK stages on every ordinary step and rank, and nonzero physical production.

## Calibration and decision gate

Each focused run is compared with its adjacent coarse run. Flag any phase change
exceeding both 1 ms and 2%, and any coarse/focused printed median differing from
the implementation's bare median by more than 3%. Keep negative differences and
unselected time in the signed comparison. Inclusive envelopes overlap and cannot
be added. Deterministic lowest-rank tie handling plus all tied ranks is recorded
for the secondary slow-rank description; it is not a proven critical path.

`status: passed` means the execution/correctness campaign completed.
`measurement_status` is separately `observer_effects_require_review` or
`ready_for_consumer_lifetime_audit`. Neither authorizes an optimization automatically.
The confirmation must establish calibrated cost plus dependency/consumer evidence
before a prototype is selected; actual savings then require bare before/after
runs. MPI time includes producer waiting and is not independently removable.

## Artifacts

- `results.json`: complete run, numerical, placement, clock and trace evidence.
- `lifecycle-summary.json`: paired calibration, signed selected locations, plan
  reuse/rebuild, transition counts and collective-site timings/volumes.
- `REVIEW.md`: concise decision and measured locations.
- `operation-sites-205.json`, `consumer-map-205.json`: source linkage and contracts.
- `operations205-rank-*.txt`, `lifecycle205-rank-*.txt`: raw per-step evidence.

Use the companion cluster protocol and return one `stage205-review-JOB.tar.gz`.
Local J4/J6 checks validate observer correctness and restart equivalence; no local
speedup or cluster calibration claim is made.

## Local validation

All four final four-rank J4/J6 observer on/off runs pass exact checkpoint
equality to the unmodified Stage202 and legacy references, the fixed legacy
tolerance screen, 12-step restart/remap coverage and all clock/trace accounting.
The enabled block trace covers all required events and each of four RK stages
on every ordinary step/rank. Both disabled tapes are empty. The actual Fortran
observer also passes O0/fcheck=all tests, and the six-build source/hash/submit
protocol passes with a compiler stand-in. See `parallel-blocks-stage205-local-results.json`.
These are local correctness checks. Cluster job 54699 subsequently passed all
41 correctness/workload gates. Both block calibration pairs pass; one legacy
wavelet phase exceeds its calibration threshold slightly, so the campaign keeps
`observer_effects_require_review`. The lifetime evidence does not establish a
justified 0.10 s/step maintenance optimization. See
`parallel-blocks-cluster-205-54699.md` and its machine-readable companion.
