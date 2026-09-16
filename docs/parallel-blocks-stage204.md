# Stage204: matched operation accounting, without numerical changes

Stage203 passed correctness but saved only 0.047–0.062 s/ordinary step, below the
predeclared 0.10 s target. Its fresh coarse residual versus legacy is about 0.144 s
in wavelets, 0.092 s in other dynamics and 0.088 s in adaptation. These are phase
budgets, not measured removable work. Stage204 measures the surviving operations
before choosing another optimization. Stage202-on remains the accepted baseline;
Stage203 is a diagnostic comparator, not a promoted implementation.

## Frozen sources and experiment

The kit reuses the exact Stage203-job-54695 numerical source inventories and
patches: Stage202-on baseline, Stage203-on experiment and pinned legacy with the
common I/O fix. Only observer patches change. Every instrumentation transformation
checks that removing its inserted calls restores the original numerical bytes.
The kit records source hashes, patch hashes and a complete call-site inventory.
Seven fresh builds retain the existing checked candidate and identical-physics
protocol. Build and run on bb13/nk, contraction off, identical compiler options,
83 ranks with 42/41 socket placement, the same J5/J7 checkpoint window.

20 solver runs, sequential on the exclusive node:

- Two checked candidate modes, oracle off/on.
- Three warmups, one per implementation.
- Two reversed bare triples: Stage202/Stage203/legacy, then legacy/Stage203/Stage202.
- Three coarse controls: the same instrumented binaries with operation timing off.
- Two reversed operation triples, with all three implementations instrumented.

This replaces the four block-only PMU runs with six matched operation runs. No
scaling survey, sampling, PMU permissions/probes or per-value accessor clocks.
Repeated runs check stability, not a precise confidence interval. The bare runs
are controls; this diagnosis has no new optimization to accept.

## Measurements and interpretation

65 named batch/call scopes partition every common context, with operation zero
holding all unassigned time. A stack attributes each clock interval exclusively
to its innermost scope. Context transitions flush that same clock ledger. Each
rank/step/context must conserve its full duration. Checkpoint steps remain in raw
files and are excluded from ordinary-step summaries.

The measured responsibilities include:

- Inverse entry/seed; scalar, outer and inner kernels; gather/restoration;
  scatter/cache updates; complete final block installs; consumer alias snapshots;
  Domain publication. Packing, ordered local aliases and receive installation
  have separate batch boundaries.
- Forward scalar/vector wavelets, restriction, compression and ghost refreshes.
- Legacy boundary packing, route construction, buffer management, local copies
  and receive installation, with actual MPI calls measured separately.
- Scalar capture preparation and transport/install; flux, restriction and
  divergence batches; scalar ghost exchange; mass transport and velocity work.
- Adaptation masks/topology updates, catalog/ownership, source construction,
  migration routes/storage/install, reconstruction retention and lifecycle work.
- MPI posting, completion and collective/blocking calls. MPI time also retains
  its immediate enclosing scope, for example inverse transfer versus Domain
  boundary installation. This is a breakdown of existing time, not extra time.

Inclusive routine envelopes are recorded alongside exclusive time. They answer
whether a complete final install or preparation remains costly after its children
have been split out. Inclusive values are nested unions, never additive. Repeated
nesting of the same label is counted once per interval in that label's inclusive
value. Exclusive signed differences alone are used to reconcile total excess.

Raw per-rank operation times and call counts remain available. The report also
uses the actual slowest full-step rank, chosen once per step, for an additive
comparison; it never sums independently chosen phase/operation maxima. Call counts
are visits to wrapped sites, including conditionally skipped calls; they are not
weighted counts of arithmetic. Function expressions and unwrapped work remain in
their enclosing scope or the explicit unassigned row. No per-node or per-slot
clocks are introduced. Arithmetic/access batches still include their internal
field access: these labels do not prove purely mathematical cost.

## Evidence gates before recommending code changes

All exact blocks/checked/observer checkpoint gates, the fixed legacy tolerance
screen, workload sequence, restart/remap and rank-placement checks remain.
Legacy numerical agreement is tolerance-based, not bitwise.

Each detailed observer is compared with its implementation's coarse control:
flag a phase change exceeding both 1 ms and 2%, or a printed step median differing
from bare by more than 3%. A failed calibration changes the diagnosis to
`observer_effects_require_review`; the resulting ranking must not be used as a
causal cost estimate. A single coarse control is a limited calibration, not proof
that per-operation perturbations cancel nowhere. The block production ledger
printing introduced in Stage203 is removed from these observer patches.

The report gives three signed rankings: Stage202 minus legacy, Stage203 minus
legacy, and Stage202 minus Stage203. Every ranking must sum to the observed
instrumented whole-step rank-mean difference. All negative rows and unassigned
rows are retained in JSON. Source inventory links labels to actual wrapped calls;
shared labels express related responsibilities, not identical mathematical work.

Do not treat a large MPI wait or copy count as removable time. First identify the
producer/dependency or the actual required consumer. Do not add producer savings
to reduced waiting. Do not infer memory stalls from clocks. No new optimization
is justified until calibrated, stable operation evidence and consumer analysis
support recovering a substantial part of the remaining gap (at least the existing
0.10 s/step design target). If the remainder or mixed kernel/access envelopes are
still large, report that ambiguity rather than claiming complete causality.

Local checks exercise the actual observer with runtime checking, nested scopes,
context changes and corrupted data rejection, the real frozen patch/build/submit
sequence, and short J4/J6 restart runs of all three observers. Local timings are
not used to forecast bb13 speed or calibrate cluster observers.

See the companion cluster protocol for submission and the single review archive.

Local validation passed for all three final observer builds: 12 steps, exact
checkpoint equality to each unmodified implementation, post-restart remapping,
and exclusive/inclusive/MPI-owner accounting. The unchanged legacy tolerance
screen also passes. See `parallel-blocks-stage204-local-results.json`. No local
performance comparison is claimed.

## Cluster result

Job 54697 completed all 20 runs and passed all 43 correctness/workload gates.
The operation ledger closes, but two observer-calibration checks fail; causal
ranking is not accepted. See [the cluster review](parallel-blocks-cluster-204-54697.md).
Stage202 remains the accepted baseline and Stage203 remains experimental.
