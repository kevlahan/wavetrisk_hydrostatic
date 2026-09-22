# Stage229 return review — job 54737

The combined Stage203 inverse workspace and scalar batching option meets the
current acceptance target on the accepted Stage226 source. Fresh uninstrumented
ABBA confirmation at 83 tasks reduces ordinary-step time from 0.806890 to
0.761174 s: **45.715 ms, or 5.666%**. The paired gains are 46.652 and 44.779 ms;
both exceed the predeclared 25 ms target and are outside its 5 ms near-threshold
band. No extra confirmation batch is required by the protocol.

This is a successful requalification of existing code, not a new algorithm.
At return review, production remained Stage226 (`1d7bc050`), with `STAGE203=off`;
the attachment review itself did not promote the candidate. The subsequent
user-authorized integration below changes the defaults. Main remains `b82647b3`
and was not modified.

## Performance evidence

| Comparison | Baseline (s/step) | Candidate (s/step) | Reduction |
|---|---:|---:|---:|
| 16-task exploratory pair | 2.570368 | 2.327113 | 243.254 ms / 9.464% |
| 83-task early transfer pair | 0.805181 | 0.761527 | 43.653 ms / 5.422% |
| 83-task independent ABBA confirmation | 0.806890 | 0.761174 | 45.715 ms / 5.666% |

Confirmation uses the mean of two separate run medians per version, excluding
screening and diagnostic runs. Ordinary steps are 2–7 and 10–12. Startup, the
saving/reload step 8 and the first full post-reload step 9 are excluded; normal
remapping at steps 5 and 10 remains included. The nine timesteps within each run
are correlated observations, not independent repetitions. This short restarting
fixture is not a continuous no-output or long-duration production benchmark.

The old Stage203/204/205 gains of 40–55 ms justified this retest, and the new
45.7 ms result is consistent with them. Current performance is established by
the new runs, rather than by subtracting historical campaign results.

## What was verified

The independent return reviewer reproduced the complete 29-run schedule and
all decisions from raw timing logs. It verified 4,373 archive member hashes,
ten source/build records and compiler logs, 23 shared physics-object identities,
all three placement probes and 1,879 rank placements. Builds and runs used the
same host, bb13, with fixed physical-core placement and FP contraction disabled.
The candidate and baseline source inventories are identical; both Stage203
compile definitions and runtime markers differ as prescribed.

All 54 numerical comparison reports pass the unchanged numerical policy; 50
are exact. The four nonexact comparisons are legacy versus blocks, with maximum
observed absolute differences of 3.724e-12 (velocity), 2.502e-12 (mass) and
7.859e-10 (temperature), below their fixed limits. Soil, topology and thresholds
agree exactly. Candidate/checked/observer block runs agree exactly with block
references, and legacy observer/repeat runs agree exactly with their own legacy
reference. Recorded final checkpoint hashes agree within each same-count family.

The J4 checked and full-oracle candidate runs at 41 tasks/40 Domains pass, as do
the J5 16-task checked/full-oracle and 83-task checked runs. Every run completes
12 steps, both checkpoint loads and both remaps with the same printed adaptive
trajectory for its fixture. This validates zero-Domain compatibility at 41/40;
it does not independently requalify Stage203 over the full former 256/160 envelope.

All 498 raw observer rank files, 5,976 rank-steps and 89,640 events pass coverage,
interval and phase-conservation checks. Four disjoint WT/tendency transactions
per ordinary step and rank are required. Cross-Python summation differences of
about 1e-16 s are accepted only in recomputed floating summaries (tolerance
1e-14); numerical acceptance limits and exact comparisons are unchanged.

The archive excludes executable, physics-object and checkpoint payloads. This
review checks their recorded identities, raw logs and comparison reports; it
does not rerun checkpoint arithmetic or regenerate binaries from the return.

## Recovered work and remaining gap

The block-only observer comparison passes the prescribed whole-step overhead
compatibility limits for both versions. Its signed differences locate the gain:

| Location | Baseline minus candidate (ms/ordinary step, rank means) |
|---|---:|
| Dynamics wavelets (`WT_after_step`) | +34.095 |
| Tendency preparation context | +13.353 |
| Adaptation/reconstruction | −1.091 |
| Other locations combined | +0.106 |
| Additive observed reduction | +46.464 |

Both reversed observer comparisons show the same main pattern: 33.773/34.417 ms
less wavelet work and 11.034/15.673 ms less tendency-context work. This supports
removing repeated representation transfers and batched access as a useful
tactic. It does not isolate the inverse and field-batching features: no ablation
was performed. Whole-step compatibility also does not prove every small phase
timer unbiased. The uninstrumented ABBA result remains the gain estimate.

The separate fresh uninstrumented diagnostic runs give:

| Version | Mean of two ordinary-step run medians (s) |
|---|---:|
| Stage226 baseline | 0.813473 |
| Stage203 candidate | 0.758268 |
| Legacy | 0.542883 |

Within this lane the candidate–legacy gap is **215.386 ms**, compared with
270.590 ms for baseline–legacy (about one fifth removed). The lane's 55.205 ms
A–B difference is descriptive and must not replace or be pooled into the
45.715 ms independent confirmation. The older 0.385 s gap is not the current
gap; accepted changes and fresh measurements must be distinguished.

**Legacy's observer fails calibration.** Its mean shift is +11.223 ms (2.067%),
above the 5.429 ms mean limit; one paired shift is +14.773 ms, above the 10.858 ms
pair limit. Baseline and candidate pass their own limits. Therefore the overall
A/B/legacy attribution is `observer_costs_unqualified`. Retain the raw signed
locations below as provisional source-search guidance, not validated exact
decomposition or recoverable savings budgets:

| Location | Candidate minus legacy (ms, provisional) |
|---|---:|
| Dynamics wavelets | +136.817 |
| Adaptation/reconstruction | +77.664 |
| Other dynamics | +19.765 |
| Tendency context | −27.882 |
| Setup, physics, remap and ordinary output combined | −4.161 |
| Additive observed difference | +202.202 |

The two observer comparisons put the wavelet excess at 135.249/138.385 ms and
adaptation/reconstruction at 77.688/77.639 ms. Their consistency supports the
search order, but cannot remove the calibration qualification. Tendency contexts
have implementation-specific boundaries; the negative row does not prove that
identical numerical kernels are faster in blocks. MPI waiting is included and
cannot be treated as independently removable communication.

The bookkeeping closes: the bare run-median gap of 215.386 ms becomes a bare
ordinary-step-mean gap of 210.969 ms (−4.417 ms), then an observer rank/step-mean
gap of 202.202 ms (−8.766 ms). These are descriptive statistic/observer bridges,
not extra routines or corrections to phase values. Positive rows alone would
overstate the gap by omitting the negative offsets.

## Recommended next action

The user subsequently authorized production integration. This commit promotes
`STAGE202=on STAGE203=on` as the Makefile defaults on the development branch.
The numerical source, checked/oracle/reference paths, zero-Domain guards and
option compatibility constraints are unchanged. Fresh build and binary
directories are required when switching settings. `STAGE202=on STAGE203=off`
retains the tested Stage226 baseline; both off retains the original path.
See `parallel-blocks-stage229-integration.json` for build-selection verification.
No new cluster performance campaign is needed for changing defaults to the exact
compile-time selection already tested by Stage229. The larger-ratio coverage
limitation above remains explicit.

Then make a bounded source comparison of the **remaining complete wavelet
representation transitions with Stage203 enabled**. Start at `WT_after_step`
(`src/adapt.f90:380`) and follow forward production, compression, inverse seed /
alias preparation and final publication into their next consumers. Compare
legacy's direct field lifetime with the block/Domain/workspace lifetimes. Check
whether a proposed reusable result also serves adaptation/reconstruction; the
two coarse excesses must not be assumed to be the same operation.

For each possible removal, establish which consumer still needs the data,
whether its topology/value epoch has changed, and the cost of the replacement.
The inverse still has real seed, alias, consumer-epoch and final-install
obligations (`src/parallel_block_mpi.f90:14682`). Their presence alone is not
evidence of redundant work. Avoid reopening Stage228's sub-ms publication copy,
Stage227's flat exchange, local-halo caches or rejected RK closure variants.

This result is enough to justify the bounded source audit without another broad
instrumentation campaign. A new performance prototype still needs a concrete
whole-pass removal and a credible net saving proportional to its complexity.
If that requires an unresolved timing fact, measure only that fact. Do not
launch an observer rerun merely to make the accounting table pass calibration.

Machine-readable evidence: `parallel-blocks-stage229-review.json`.
Reproducer: `test/parallel_block_profile/review_cluster_229.py`.

Archive SHA256:
`d842d14d530b69f1ece20d6a05748b88145b5c34d80a4b8e68a09257d345b3a5`.
