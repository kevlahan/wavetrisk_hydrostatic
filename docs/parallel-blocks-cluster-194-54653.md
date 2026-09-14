# Stage 194 cluster review: job 54653

Accept Stage 194 for the tested configuration and commit the implementation,
split profiling, cluster tooling and evidence. The optimization produces a small,
consistent timestep gain and a large reduction in the intended inner-vector
kernel. This review changes no numerical source or runtime defaults.

## Validation and provenance

The supplied archive `stage194-review-54653.tar.gz` has SHA-256
`8f015e1cfe2d58ba7f69ecdfd29eb5cc12aecfe334827244a39192e4aa43bbeb`.
All 4,853 inventory entries passed byte-count and hash checks. Baseline and
candidate source identities match the frozen kit and current numerical source.
The baseline is validated Stage 193, reconstructed from the committed Stage 191
bootstrap; it is not the unmodified bootstrap executable.

- Baseline patch: `e00d6dbb6d010c34c1f0fcc56b0d127d80551cda7d4174e0078cc7527d0cfc33`.
- Candidate delta: `ab37d1f6b458aa233aa6cfeafdeda9b430bffdad696881dda478bdff11944796`.
- Five build records and compiler logs match; optimized commands use O2,
  `-march=native`, and contraction off; checked commands use O0 and runtime checks.
- J5/J7, RK4, 83 ranks on bb13, with identical 23 external physics object hashes
  recorded by the builder. External physics is outside the solver runtime-check scope.
- All ranks retained the same 83 distinct single-core bindings across runs.
  No swap-in/out was recorded.

All 28 runs completed the same 12-step workload, read checkpoint 3, wrote and
reloaded checkpoint 4, and continued through remapping. All 64 gates passed:
30 exact checkpoint reports, 26 workload comparisons and eight fixed legacy
screens. The audit reparses the logs and reevaluates the recorded comparison
reports. Dynamics, adaptation and remap oracles were enabled in checked-on,
with Stage 193 scalar columns and Stage 194 geometry caching enabled.

Maximum legacy differences remained velocity 3.72324e-12, mass 2.50111e-12,
and temperature 7.85803e-10, within the established limits 1e-10, 1e-10 and
1e-7 respectively. Legacy screening is not bitwise equivalence. Checkpoints and
executables are absent from the review archive; this audit does not independently
rerun those payload comparisons. Saved checkpoints cover the restart boundary,
with later workload/remap coverage rather than a final-state checkpoint.

## Repeated timing

Six rotating triples, after warmups; controls and profiles are excluded:

| Measure | Stage 193 | Stage 194 | Legacy |
|---|---:|---:|---:|
| Median ordinary-step time | 1.280 s | 1.260 s | 0.558 s |
| Median elapsed time | 31.0901 s | 30.8695 s | 22.5768 s |
| Ordinary-step ratio to legacy | 2.2939x | 2.2581x | 1x |

Stage 194 reduces ordinary-step time by approximately **1.56%**, and elapsed
time by **0.71%**. All six paired ordinary-step comparisons improve by
0.01–0.03 s, and all six elapsed comparisons improve by 0.192–0.363 s.
Elapsed ranges do not overlap: baseline 31.029–31.261 s, candidate
30.764–30.950 s. Printed step times have finite precision; do not overinterpret
the percentage digits. The single cache-disabled candidate control returns to
1.28 s per ordinary step (30.9538 s elapsed).

This removes about 2.77% of the previous excess timestep time above legacy.
The remaining 2.26x gap is substantial. Memory is essentially unchanged:
median sums of rank peak RSS are 21.6734 GiB baseline, 21.6829 GiB candidate,
and 11.0116 GiB legacy. These are sums of separate peaks, not simultaneous RSS.

## The split identifies the improvement

Same Stage 194 binary and split instrumentation, cache disabled versus enabled.
Whole-window rank-mean kernel self times include the checkpoint step:

| Scope | Cache off | Cache on |
|---|---:|---:|
| Scalar inverse | 0.0460749 s | 0.0461548 s |
| Inner-vector inverse | 0.2412416 s | 0.1008489 s |
| Outer-vector inverse | 0.0809681 s | 0.0810857 s |
| Combined parent, inclusive | 0.3692180 s | 0.2290200 s |

Inner-vector time falls **58.20%**, with essentially unchanged scalar and
outer-vector work. There are 7,636 inner, 7,636 outer, and 15,272 scalar calls
in both profiles. Children plus parent overhead conserve the combined parent.
The original Stage 193 parent measures 0.361552 s; the disabled split control
measures 0.369218 s. This small difference includes instrumentation and run
variation; one profile cannot isolate their separate contributions.

Physical cell updates remain 11,592,480. Geometry reads fall from 370,959,360
to 12,365,312 (**96.67%**), preparing 386,416 spatial cells once across 30
layers. Actual disabled-control reads match the candidate reference count.

Matched complete wavelet time falls from 4.88730 to 4.65359 s (**4.78%**);
legacy measures 1.28361 s, leaving **3.63x** for this matched transaction.
Inverse gather/alias/scatter retain exactly 42,263,016,000 remote bytes and
496,004 messages. Inverse MPI-wait self time falls from 1.143164 to 1.074055 s
despite unchanged traffic, consistent with less local work reducing waiting.
This is diagnostic evidence, not a bandwidth measurement.

## Implications for the next computational stage

In the candidate's ordinary-step profile, scalar, inner and outer inverse
kernels plus their parent overhead total only **0.01881 s**, about **1.49%**
of the 1.25894 s rank-mean step. Inner-vector alone is 0.00828 s, about 0.66%.
Further tuning these kernels cannot by itself close the 2.26x timestep gap.

Larger computational scopes include scalar replay self time (0.09640 s),
native wavelet work (0.06091 s), and block ghost pack/install (0.05556 s).
Scalar input/restriction as a broader exclusive group totals 0.20043 s,
including its separately measured transport. Legacy component boundaries differ;
subtracting its scalar group from this one would not measure duplicate work.

The next target should be the direct-flux/restriction production chain. Split
direct-flux preparation, flux restriction, divergence and temperature-boundary
computation before attributing all scalar replay to one kernel. Source inspection
already provides two specific candidates in `compute_block_scalar_direct_flux`:

1. A 50-value record is read before checking the production mass path, which
   uses only its physics entry; it is also read before rejecting nonphysical
   levels for the other scalar fields. Narrow the data reads while preserving
   required mass/soil handling and validation coverage.
2. Horizontal neighbour selection repeats inside the scalar/layer loops.
   Evaluate preparing addresses once per spatial stencil and consuming the
   live field values across the physical column, preserving expression order.

Neither candidate has a measured individual cost yet. Production restriction
and divergence already have column paths; do not assume that work is still the
older cell-at-a-time implementation. Oracle comparison branches are disabled in
timing runs. The word "replay" denotes native flux/restriction production and
does not establish duplicate legacy execution. Compatibility boundary inputs
remain, but their avoidable share requires a producer/consumer audit.

Keep direct Stage 194 and legacy comparisons, unchanged numerical gates, and
disabled controls for the next stage. Retain the split inverse timers; they have
successfully resolved the previous ambiguity. This campaign alone supports no
claim for other rank counts, multiple nodes or longer integrations.

## Reproduce the review

`parallel-blocks-cluster-194-54653-campaign.json` preserves the supplied campaign.
`parallel-blocks-cluster-194-54653-analysis.json` records the independent audit.
With the extracted archive, frozen kit and tested source available:

```sh
python3 test/parallel_block_profile/analyze_stage194_campaign.py \
  /private/tmp/stage194-review-54653 \
  --kit /private/tmp/wavetrisk-stage194-cluster-kit
```

The audit intentionally checks current numerical source hashes; after further
source changes, use the Stage 194 source snapshot with `--repo`.
