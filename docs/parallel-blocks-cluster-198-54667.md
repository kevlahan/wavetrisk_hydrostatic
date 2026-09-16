# Stage 198 cluster diagnosis: job 54667

Job 54667 completed the diagnosis-only Stage 198 campaign on the exclusive
`nk` node `bb13` with 83 ranks, J5/J7, RK4 and contraction disabled. The
numerical candidate is the validated Stage 197 implementation; Stage 196 is the
untouched control and the legacy executable is the pinned Domain implementation
with the common checkpoint-I/O fix.

The review archive SHA-256 is
`b77255e0b10867b1868c2fc38fff2a0ad8e0ab3b0e5ace6c556d8761d21b50f1`.
All 9,916 archived files match `collection.json`. All 27 executions and all 65
numerical/workload gates passed. The seven executables use 23 identical external
physics objects. All ranks ran on distinct physical cores of `bb13`; node swap
input/output and NUMA miss/foreign counters remained zero. The block executables
used about 21.7 GiB summed rank-peak RSS versus 11.0 GiB for legacy. This confirms
the footprint difference but does not by itself establish a bandwidth or NUMA
cause.

## Reconciled full-step result

The bare medians are 1.260 s for Stage 196, 1.245 s for Stage 197 and 0.556 s
for legacy. Stage 197 is therefore 2.239 times legacy, with a 0.689 s/step bare
excess. It retains the previously measured 1.19% improvement over Stage 196.

The exclusive operation ledger uses rank means from the timing-observer runs.
It measures 0.637382 s excess; its signed rows sum to the same value within
roundoff. The actual slowest-full-step-rank path gives 0.638032 s excess. These
are not rescaled to the 0.689 s bare-median gap because they are different
statistics and the observer affects legacy more strongly.

| Full-step task | Legacy, ms | Block, ms | Excess, ms | Share of instrumented excess |
|---|---:|---:|---:|---:|
| Dynamics wavelet transaction | 78.602 | 341.853 | 263.251 | 41.3% |
| Dynamics outside wavelets | 405.003 | 666.288 | 261.284 | 41.0% |
| Adaptation/reconstruction | 64.736 | 179.834 | 115.098 | 18.1% |
| Physics | 12.185 | 12.448 | 0.263 | 0.04% |
| Setup, remap and output, signed | 9.952 | 7.439 | -2.514 | -0.39% |
| Total | 570.479 | 1,207.862 | 637.382 | 100% |

The two rounds are stable. Their respective excesses are 262.231/264.271 ms
for wavelets, 262.201/260.368 ms outside wavelets, and 115.796/114.399 ms for
adaptation/reconstruction. This reproduces the Stage 196 location ledger after
allowing for the approximately 15--17 ms Stage 197 wavelet improvement.

## Complete signed classification

The following classification partitions all 101 operation IDs. It is a timing
classification, not a claim that every positive interval is removable.

| Operation family | Positive rows, ms | Block-faster offsets, ms | Signed excess, ms |
|---|---:|---:|---:|
| MPI completion/posting and waits | 330.572 | -96.811 | 233.761 |
| Boundary/transfer preparation, packing, copies and installation | 179.613 | -62.728 | 116.886 |
| Wavelet and numerical computation/traversal | 206.184 | -40.331 | 165.853 |
| Scalar storage, geometry and capture | 33.331 | 0 | 33.331 |
| RK/state/compatibility work | 31.579 | -5.371 | 26.208 |
| Adaptation construction and migration CPU work | 63.759 | 0 | 63.759 |
| Other/unclassified | 3.452 | -5.866 | -2.414 |
| Total | 848.489 | -211.107 | 637.382 |

Thus about 403.6 ms of the measured net gap is non-MPI CPU work and about
233.8 ms is recorded in MPI/post/wait intervals. The latter must not be read as
233.8 ms of avoidable network cost.

## Why the large wait rows are not a message-volume conclusion

The four largest wait/MPI rows are inverse wait (90.090 ms), native-mass wait
(75.189 ms), scalar transport MPI (46.070 ms), and migration MPI (36.231 ms).
Across the 83 ranks on every ordinary step, each interval is strongly
anticorrelated with other computation in its enclosing context:

| Interval | Median correlation with other context computation | Range over the two runs and ordinary steps |
|---|---:|---:|
| Inverse wait | -0.989 | -0.993 to -0.980 |
| Native-mass wait | -0.982 | -0.989 to -0.953 |
| Scalar transport MPI | -0.972 | -0.993 to -0.954 |
| Migration MPI | -0.942 | -0.962 to -0.904 |

The sum of non-MPI work and communication is consequently much more balanced
than either component. For example, candidate tendency computation averages
374.3 ms/rank but the most heavily computing rank is about 148.1 ms above that
mean. Communication averages 252.6 ms and ranges from about 101.7 to 458.8 ms;
their sum has only a 29.4 ms rank range. The analogous compute/communication
figures are 189.4/152.5 ms for dynamics wavelets and 102.7/50.2 ms for adaptation
remainder.

This is consistent with synchronization absorbing rank-local work imbalance.
It does not show that the transport is free, nor separate arrival delay from
network latency, but it rules out assigning the entire wait interval to message
volume. Work-volume correlations identify concrete sources of the imbalance:
wavelet non-MPI time correlates about 0.955 with inverse values sent; tendency
non-MPI time correlates about 0.850 with boundary candidates and 0.833 with
native-mass send values; adaptation CPU time correlates about 0.756 with migrated
bytes. These are correlations, not independently recoverable-time estimates.

Absolute spread alone must not be called worse load balance. Normalize the
maximum-minus-mean CPU time by the corresponding mean CPU time:

| Context | Legacy mean CPU, ms | Block mean CPU, ms | Legacy max-minus-mean, ms | Block max-minus-mean, ms | Relative excess, legacy / block |
|---|---:|---:|---:|---:|---:|
| Dynamics wavelets | 51.1 | 189.4 | 38.8 | 92.5 | 76% / 49% |
| Tendency tag | 199.7 | 374.3 | 131.4 | 148.1 | 66% / 40% |

On this aggregate context measure, blocks are relatively more balanced, while
performing much more CPU work. The larger absolute wavelet spread therefore does
not demonstrate a degraded balancing algorithm. Neither these aggregate spreads
nor negative compute/wait correlations resolve individual dependency arrival
causes. Extra representation work remains the primary demonstrated difference;
load-model inadequacy is a hypothesis requiring producer-level evidence.

The native inverse performs about 177.27 dependency completions per rank/step.
Its 90 ms wait is distributed over gather, scatter and boundary dependencies:
about 18.0, 31.8 and 38.8 ms/rank/step respectively. No single message dominates.
This is consistent with earlier volume-only changes producing little end-to-end gain.

## What makes the three parent regions slower

The 263.251 ms dynamics-wavelet excess decomposes exactly into:

| Wavelet mechanism | Signed excess, ms |
|---|---:|
| Native inverse transfer plan/pack/alias/post/wait/install | 149.268 |
| Block ghost packing/install and MPI | 34.179 |
| Domain compatibility import/writeback/publication | 46.731 |
| Wavelet/compression/inverse numerical work | 75.755 |
| Generic Domain boundary path saving | -42.957 |
| Other | 0.275 |

The 227.604 ms tendency-tag excess (the remaining 33.680 ms is dynamics/RK
outside that tag) decomposes into:

| Tendency mechanism | Signed excess, ms |
|---|---:|
| Block scalar capture/restriction/divergence pipeline | 156.948 |
| Native mass pipeline | 86.116 |
| Shared/native producers and preparatory operators | 38.465 |
| RK/state compatibility | 22.281 |
| Block ghost path | 31.010 |
| Domain compatibility imports | 7.140 |
| Generic Domain boundary path saving | -114.363 |
| Other | 0.009 |

The 115.098 ms adaptation/reconstruction excess is dominated by block lifecycle
work: catalog/source construction, packing, migration, installation and residual
management total 99.715 ms. Of this, migration MPI is 36.231 ms; the remaining
63.484 ms is CPU work. Domain reconstruction/import/ghost work contributes
12.604 ms, mask/compression arithmetic 1.778 ms, and other work 1.171 ms.

The individual top rows are stable between timing rounds. The first ten positive
rows total 421.219 ms, but they cannot be treated as 421 ms recoverable because
the full ledger also contains 211.107 ms of block-faster offsets and waits that
overlap the causal effect of rank imbalance.

## Numerical duplication versus representation overhead

The measurements do not support a wholesale reproduction of the legacy numerical
timestep. Selected work counts are lower in blocks: restricted scalar-flux edges
are 7,588 versus 45,527 per rank/step, divergence evaluations 35,559 versus
83,383, and temperature direct-flux evaluations 60,679 versus 79,842. Direct
scalar node-field samples are close (165,855 versus 159,684). Generic Domain
boundary sends, receives and candidate scans are approximately half their legacy
counts, yielding a 154.705 ms signed saving across the full step.

The systemic cost instead comes from maintaining the block representation and
its compatibility interfaces. The scalar pipeline alone contributes 156.948 ms:
scalar setup/storage 26.163 ms, direct-flux traversal excess 30.505 ms, transport
and installation 74.943 ms, replay excess 24.415 ms, and smaller closure/restriction
terms. Native mass contributes another 86.116 ms. Calls to `basic_operators` are
the same (349.1/rank/step), but the block routine also captures scalar physics,
captures divergence inputs, executes native mass restriction and initiates native
mass exchange; its exclusive excess is 20.970 ms. The surrounding shared producer
adds 17.402 ms and equal-count diffusion preparation adds 15.373 ms. Equal call
counts therefore do not mean equal work or equal memory state.

The storage census confirms a highly fragmented representation: approximately
842,859 range requests, 10.95 million logical requested slots, 10.58 million
`read_value` visits, 7.52 million shared-address resolutions and 3.12 million
`write_value` visits per rank/step. The sampled inclusive estimates must not be
added to the ledger. Merely compiling the disabled probe hooks changes the step
median by 4.78%, and enabling samples adds another 5.70%; two clock reads are a
large fraction of these tiny functions. The samples establish frequency and
nested access structure, but their absolute 138/436/264/121 ms estimates are too
observer-sensitive to claim as production costs.

## Calibration and remaining causal work

The minimal common clock changes printed medians by about 0.8% for both sides.
The full timer observer changes the candidate by 2.01% but legacy by 6.38%.
Legacy parent envelopes flag wavelets, outside-wavelet dynamics and adaptation;
the timing ledger therefore understates the bare gap and cannot be corrected by
one global factor. Work counters add little further parent distortion after the
timer pass.

Stage 198 fully accounts for where the instrumented time goes and identifies the
three systemic mechanisms: block-native/compatibility pipelines, rank-local work
imbalance exposed at frequent synchronization points, and per-step block lifecycle
construction/migration. It still does not establish how much of each mechanism is
avoidable while preserving the block design. In particular, wait time cannot be
added to producer savings, and storage sample time is not a recoverable estimate.

Before optimizing, the next diagnosis should split the large CPU envelopes at
batch granularity without per-element clocks:

1. Partition `basic_operators` into the legacy-equivalent pressure/`step1` work,
   scalar capture, divergence capture, native-mass restriction/exchange, and
   Bernoulli/Exner restriction.
2. For every `(layer, level)` scalar and mass exchange, record preceding producer
   duration and work volume beside its wait. This will identify which producer
   makes the lagging rank late without introducing barriers.
3. Partition the 63.5 ms adaptation CPU cost into catalog traversal, allocation,
   route construction, serialization, local installation and retained-state work.
4. Replace per-access storage clocks with outer-batch timers and counter-only
   access censuses. Any isolated replay must use the same records and preserve
   the measured producer/consumer sequence.

No production optimization should be selected from job 54667 alone. A candidate
needs a measured avoidable CPU interval and a conservative completion-time saving
that accounts for the corresponding reduction in MPI wait exactly once.
