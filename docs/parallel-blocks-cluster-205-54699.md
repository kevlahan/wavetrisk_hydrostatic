# Stage205 cluster review — job 54699

Stage205 passes execution and correctness. It provides useful negative evidence:
repeated scalar-plan construction, repeated geometry seeding, direct Domain/block
transactions, and migration checks do **not establish a justified 0.10 s/step
optimization**. It does not finish the causal explanation of the remaining gap.
Stage202 remains the accepted baseline; Stage203 remains experimental.

## Validation

Archive SHA256: `273b4e25ce387b34e662f8781d07f866b35ddd05d1345013b3c97d96142dcc62`.
All 6,528 inventoried files match their sizes/hashes, without missing or extra
files. Six recorded source inventories and build-log hashes match the frozen kit;
every run's recorded executable identity matches its build. Source-site maps and
consumer contracts match the kit. Builds use J5, contraction off and identical
23 physics objects. Placement is 83 ranks, 42/41 across bb13's sockets.

All 19 runs complete the same 12-step workload, restart from checkpoint 3,
write/reload checkpoint 4, and remap after restart. All 41 archived gates pass:
checked oracle off/on, exact block/observer comparisons, unchanged legacy
tolerance checks, and matching workloads. Legacy agreement is not bitwise.
Checkpoint and executable payloads are absent from the archive: their recorded
identities and numerical comparator reports were checked, rather than rerunning
payload comparisons.

Independent reconstruction reproduces raw log analyses, placement, resource
records, hardware-clock summaries, operation ledgers, lifecycle traces and the
saved summary. Exclusive clocks and MPI owners conserve enclosing clocks.
Stage205's deterministic slow-rank/tie results reproduce. The older secondary
scaling analyzer can select a different equally slow rank; these cases were
reconciled against exactly equal raw whole-step clocks, retaining the archived
selected rank's phase values. No unique causal critical path is inferred.

Bare runs show no swap and at least 1,489.85 GiB available memory. Sums of per-rank
maximum RSS are about 21.8 GiB for blocks and 11.1 GiB for legacy; these are not
simultaneous aggregate peaks. Capacity pressure is not the explanation here.
Cache/bandwidth effects remain possible.

## Bare performance

| Implementation | Reversed samples, s/ordinary step | Median | Excess over legacy |
|---|---|---:|---:|
| Stage202 accepted baseline | 0.925, 0.930 | 0.9275 | 0.3845 |
| Stage203 experimental | 0.878, 0.881 | 0.8795 | 0.3365 |
| Legacy Domain | 0.542, 0.544 | 0.5430 | — |

Stage203 saves 0.048 s/step (5.18%); paired differences are 0.047 and 0.049 s.
This is repeatable but below the agreed 0.10 s criterion. Current ratios are
1.708x and 1.620x legacy, respectively. Two repetitions are not a precise
statistical confidence estimate.

## Calibration qualification

Both baseline observer/control pairs pass every existing calibration threshold.
The first legacy pair also passes. In the second legacy pair, dynamics-wavelet
time rises from 0.069930915 to 0.071332838 s/rank-step: +1.401923 ms, or 2.004726%.
The threshold is 1.398618 ms, so the measurement exceeds it by 3.305 microseconds.
All other phase and printed-step checks pass.

Retain `observer_effects_require_review`; do not retrospectively loosen the
threshold or report full calibration. This small flag does not erase the stable
structural counts or passing block-cost controls. It also does not justify
another complete campaign simply to turn the status green. Any later numerical
wavelet experiment should include its own adjacent legacy control. The phase
table below describes measured locations, not certified removable costs.

## What the lifetime confirmation establishes

Costs below are Stage202 rank-mean elapsed seconds per ordinary step, averaged
over two operation runs. Plan/import/publication/capture rows combine exclusive
local time and MPI attributed to those same owners. They do not sum nested
inclusive parents. These are absolute block costs, not like-for-like net legacy
excess or promised savings.

| Target | Observed lifetime/frequency | Measured cost | Decision |
|---|---|---:|---|
| Scalar capture, restriction and temperature plans | Each has four requests per rank-step: one rebuild, three reuses | All selected plan work, including ghost/writeback plans: 0.03262 | No evidence of rebuilding these plans at every RK stage. Persisting across changed topology would need a separate validity proof. |
| Scalar geometry preparation | One seed per rank-step; no repeated seed within a generation | 0.00031 | Not a substantial target. This narrow timer is not all geometry computation in the solver. |
| Scalar capture setup/finish | One full capture, then three compact captures per rank-step | 0.02577, excluding nested plan/geometry work | Later captures transport current physical residuals; they are not repeated full geometry captures. |
| Domain-to-block imports | Four native velocity-tendency imports, four vector-wavelet imports, one edge-length import per rank-step | 0.01105 | No repeated visit within the recorded stage/generation for these families. |
| Block-to-Domain publications | Four temperature publications, about 9.18 full-SOL publications and five WAV publications per rank-step | 0.02931 | Repeated stage/generation is not an unchanged-value proof. Total direct import plus publication cost is only 0.04036. |
| Migration checksum scans and 13 validation Allreduces | Once per migration/rank-step | 0.00728 scan + 0.01126 reductions = 0.01854 | Too small to support a 0.10 s check-removal target. Preserve correctness/error ordering. |

Each repeat contains 913 ordinary rank-steps (83 x 11). For each of the three
scalar plan types, this means 3,652 requests, 913 rebuilds and 2,739 reuses.
Capture flags show 913 full and 2,739 compact captures. Stage-one storage is
rebuilt in 484 cases and retained in 429; the later stages retain it in all cases.
The same structural counts reproduce in both repeats.

The 13 validation reductions are separated from actual payload movement. Almost
all their time is in the first checksum Allreduce: 10.876 and 11.297 ms in the
two repeats. Across rank-step samples, its duration correlates negatively with
the preceding checksum-scan duration (-0.792 and -0.799). This is consistent with
waiting for uneven scan completion, not proof of intrinsically expensive small
messages. Local clocks are not globally synchronized. Scan plus all validation
reductions never exceeds 33.66 ms in any recorded rank-step; neither the mean
nor that observed maximum supplies a 100 ms savings case. Do not count producer
time and the resulting waits as independent step savings.

The whole migration envelope is 58.14 ms and includes nested checks, exchange,
installation and plan work; it must not be added to these rows. The broader
validation/size/payload scopes include a real payload exchange and must not be
mislabelled as disposable checking.

## Source interpretation

`prepare_block_multistage_tendency` in `src/time_integr.f90` calls the retained
compatibility producer, scalar capture, native velocity preparation and native
tendency evaluation. The full Domain tendency call is conditional on the oracle.
`block_tendency_compatibility_ml` in `src/multi_level.f90` likewise gates the old
mass/velocity comparisons on validation. The production evidence does not show
a complete duplicate legacy tendency execution.

The operation named **RK tendency capture** has 77.33 ms exclusive time and a
122.10 ms inclusive envelope. Its body calls
`evaluate_candidate_block_velocity_recomposition`, scalar restriction, and
`evaluate_candidate_block_exner_shadow`, which applies the physical tendency
kernel and retains the resulting tendency. Therefore this timer is not a
77 ms copy/adapter cost, and the 122 ms parent is not an available savings budget.
The retained-remainder flag also prevents a second preparation of the velocity
remainder in this call path.

Publications serve different source contracts: native temperature is published
after RK assembly for Domain-shaped wavelet consumers; wavelet coefficients are
published after transformation/compression; SOL is published after inverse
reconstruction; grid-change synchronization requests SOL and WAV before the
representation changes. These are in `retain_block_native_multistage_candidate`,
wavelet activation/inverse routines and `synchronize_parallel_block_checkpoint`
in `src/parallel_block_mpi.f90`. The final RK marker persists into later contexts,
so equal stage numbers can include both dynamics and later synchronization.
The trace is not a pointwise def-use/version proof: this review neither declares
every publication indispensable nor labels all repeated visits redundant.

## Remaining gap and next decision

The signed exclusive rank-mean operation ledger closes to 0.371101 s. Its phase
partition is:

| Complete phase | Stage202 minus legacy, s/rank-step |
|---|---:|
| Dynamics wavelets | 0.173283 |
| Other dynamics | 0.111863 |
| Adaptation/reconstruction | 0.088999 |
| Other phases | -0.003045 |
| Total | 0.371101 |

This instrumented rank-mean quantity differs from the 0.3845 s difference of
bare printed medians. Negative rows are retained. Unselected time includes
necessary computation and unobserved MPI; zero time under a legacy label does
not mean legacy performs no equivalent work.

The confirmation has **not established a coherent, avoidable 0.10 s cost**.
Do not implement a persistent-plan, geometry-cache, check-removal or generic
message-volume optimization on these results. Do not combine unrelated totals
into a fictional savings budget.

The next investigation should be a source-led, like-for-like audit of the
**complete forward/inverse dynamics-wavelet numerical paths**, the largest
remaining location, using existing Stage204/205 records first. Match legacy and
block mathematical outputs, count actual traversals and temporary preparation,
and distinguish repeated calculation, indirect scalar access, local data copies
and required boundary dependencies. Then examine the physical computation inside
the misleadingly named RK capture envelope. Scope a new measurement only to a
specific unresolved mechanism; another broad lifecycle campaign is not warranted.
A prototype still requires a dependency-preserving replacement with evidence
supporting the 0.10 s target, followed by bare timing and correctness tests.

Machine-readable audited results and derived budgets are in
`parallel-blocks-cluster-205-54699.json`. This review changes documentation only;
it makes no numerical changes and no commit.
