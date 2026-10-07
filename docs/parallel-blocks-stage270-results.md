# Stage270 review: useful block cost decomposition; main calibration still fails

Job 55467 completed all fifteen runs. Independent local review verified all
2,540 archive files and reproduced the science, placement, storage, timing and
raw observer checks. The archived identity, source manifest and kit manifest
match the delivered Stage270 kit byte for byte. Five serial builds shared 23
physics objects. Checked execution and both 16/167-rank independent oracles
passed, retaining 48 exact inverse-return proofs and 84 installations each.
All checkpoint, 62-member VTK, restart/remap and zero-Domain-rank checks passed.
Production remains Stage266 `ff740744`; this campaign tested instrumentation,
not an optimization or a further Domain removal.

## Total campaign runtime

The user supplied this Slurm accounting result:

| Job | State | Start | End | Elapsed |
|---|---|---|---|---|
| 55467 | COMPLETED | 2026-10-07 17:40:43 | 2026-10-07 17:56:21 | **00:15:38** |

The total is **938 seconds**, including the five builds and fifteen runs with
their validation work. It excludes queue wait and later collection/transfer.
The archive itself contains no Slurm accounting record; this elapsed time is
user-provided evidence, not inferred from file timestamps. Preserve it in the
[campaign runtime ledger](parallel-blocks-campaign-runtimes.json).

The earlier 30-minute estimate (20–45-minute range) was conservative. Budget
about **20–25 minutes after allocation** for a comparable campaign using the
same fixture and schedule. This is a planning allowance, not a prediction for
larger or multi-node tests. Future kits should collect Slurm Start/End/Elapsed
and explicit campaign timestamps so reviews do not require a separate query.

## Bare controls and calibration

| Measurement at 83 ranks | Seconds per ordinary step |
|---|---:|
| Legacy main, mean of two run medians | 0.540342276 |
| Accepted Stage266, mean of two run medians | 0.749653039 |
| Candidate minus main | **0.209310763 (+38.737%)** |

This remains a substantial gap. The source is unchanged from the previous
campaign; the small difference from Stage269's 212.308 ms gap does not establish
a performance improvement. Steps 1, 8 and 9 remain excluded by the fixed policy.

| Observer | Mean shift | Relative shift | Fixed calibration decision |
|---|---:|---:|---|
| Main | +8.715 ms | +1.613% | Failed |
| Stage266 | +2.758 ms | +0.368% | Passed |

Main's mean shift exceeds its 5.403 ms limit; the second paired shift is
12.348 ms, exceeding the 10.807 ms pair limit. The first pair is +5.083 ms.
Both spreads pass. Candidate shifts are +1.426 and +4.091 ms and pass all
prespecified limits. Main executes only 23 observer scopes per ordinary step;
the block observer executes about 467.794. Reduced instrumentation frequency
has improved compatibility, but does not prove the remaining main shift is
direct clock-call overhead; code generation/layout and run variability have
not been isolated.

Cross-version phase differences remain **unqualified**. Do not subtract a flat
main correction, relax the limits, or present signed region differences as a
measured decomposition of the 209.311 ms gap. Candidate calibration supports
using its substantial component costs to choose investigations. It does not
establish unbiased small-scope costs or recoverable savings.

## What the new decomposition establishes

The following values average ranks and ordinary steps across the two optimized
block observer runs. Self time excludes separately measured children; inclusive
time includes them. Nested inclusive costs must not be added together.

| Candidate operation | Cost per step | Interpretation |
|---|---:|---|
| Shared compatibility producer | **193.987 ms self; 195.821 ms inclusive** | Largest measured individual region; still contains Domain-owned shared producers and necessary science |
| Scalar capture preparation and publication | **43.906 ms combined inclusive** | Two separate sequential transactions; includes 2.244 ms capture MPI |
| Scalar restriction transaction | **87.574 ms inclusive** | Most of the 109.260 ms recomposition transaction |
| Final physical kernel batch | **19.623 ms self** | Much smaller than the old undivided recomposition label |
| Migration | **64.341 ms inclusive** | Transport, proof work, packing and preparation, not just network transfer |
| Inverse MPI wait/traffic | **51.746 ms self** | Includes waiting for other ranks; not a standalone removable budget |
| Native scalar forward | **51.267 ms self; 55.486 ms inclusive** | Significant transform work remains |

The scalar restriction transaction contains 36.017 ms local-kernel self time,
10.787 ms packing/local copies, 16.101 ms MPI posting/waiting, 5.181 ms received
installation, and a 19.170 ms temperature-edge transaction (including 13.274 ms
waiting). Some local-kernel inclusive time includes nested edge work; use self
times for additive accounting. These measurements redirect attention from
assuming all recomposition cost is Exner arithmetic to its restriction and
transport chain.

Migration's components include:

| Component | Mean self time per step |
|---|---:|
| Proof reductions | 13.137 ms |
| Manifest construction/check | 9.949 ms |
| Writeback plan preparation | 9.135 ms |
| Payload MPI | 8.163 ms |
| Byte-checksum scans | 7.403 ms |
| Packing batch | 6.473 ms |
| Ghost plan preparation | 4.965 ms |
| Unpacking batch | 3.351 ms |

Checksum scans and proof reductions together account for 20.540 ms of mean-rank
self time. This confirms that these production checks deserve a bounded source
audit, but **does not promise a 20.540 ms speedup**. Per-rank mean reduction time
ranges from 1.934 to 31.143 ms, consistent with different arrival times as well
as collective cost. Removing upstream imbalance can move time between scopes.
Combining compatible reductions while retaining the integer conservation check,
or integrating checksum accumulation with required traversals, are hypotheses
to qualify separately; do not simply disable safeguards.

Actual work counters show one rebuild and one reuse per ordinary step for each
of the ghost and writeback plans, and one rebuild plus three reuses for edge
metrics. This does not prove the rebuilt geometry or ownership was unchanged.
Any reuse change must establish valid topology, geometry and payload epochs,
including partially changed blocks, while preserving the post-physics refresh.

The component costs are consistent across both candidate runs: shared-producer
self time is 193.855/194.119 ms, restriction inclusive time 87.668/87.480 ms,
and migration inclusive time 64.741/63.942 ms. Per-rank local-kernel work and
inverse waits vary substantially. The saved analysis also retains each step's
longest-local-duration rank composition; it is not a synchronized global
critical path. Process CPU includes MPI spinning.

## Recommendation in light of the fourfold scaling goal

Use the measured **shared Domain producer → capture → block restriction** chain
as the primary architectural investigation. Its importance is both absolute
cost and the restriction of some useful work to Domain owners. Much of its
arithmetic is necessary, so migration must supply a complete block-owned
replacement with the established numerical ordering and consumer contracts;
the region's full cost is not overhead that can simply disappear. Retain the
Stage267 publication/reader constraints until those consumers are replaced.

The migration proof/packing/preparation chain supplies a smaller, concrete
optimization opportunity that can be tested independently. Select one bounded
change after auditing its dependencies; do not combine it with a broad producer
rewrite or accumulate multiple changes before measurement. Cache construction
(2.109 ms) and surface pressure (1.547 ms) remain lower priorities for closing
this gap or removing the principal scaling barriers.

The next performance campaign should include an uninstrumented J5 baseline at
**83, 160 and 256 ranks**, with matched legacy main where supported and any
candidate at the same counts. Rank-specific diagnostic comparisons require
their own calibration; the 83-rank result does not certify overhead at 256.
The 167-rank Stage270 run is an oracle and cannot establish production scaling.
There is no need to repeat this unchanged broad audit merely to seek a passing
main calibration before acting on well-supported candidate-only leads.

This is an intermediate screen: the strategic acceptance goal remains efficient
execution at **640 or more ranks for J5's 160 geometric Domains**, on multiple
nodes with appropriate placement and adequate block subdivision. Neither the
current 83-rank performance nor the successful zero-Domain-rank oracle qualifies
that goal. No optimization, speedup forecast or Domain removal is promoted by
this review.

The [machine-readable results](parallel-blocks-stage270-results.json) retain the
verified comparison, candidate component/work records, rank ranges and runtime.
Archive SHA256:
`953f1e57e2137e9baa678245618e6ba0c03e6bcb5248c6ac7662099099e268cb`.
The review archive omits binary/checkpoint/VTK payloads: their cluster-generated
hashes and numerical comparison reports were verified, not the omitted payloads
rerun locally. Eight-ulp tolerance applies only to derived floating reports;
raw hashes, counts and scientific gates remain exact.
