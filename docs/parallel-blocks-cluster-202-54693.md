# Stage202 cluster review — job 54693

**Stage202 meets the substantial-gain target and is ready to commit with this
validation record.** The combined prototype removes 0.130 s per ordinary step,
a 12.3% time reduction and 25.3% of the measured excess over legacy. Both reversed
pairs exceed the predeclared 0.10 s/step threshold. This review does not change
the source default or make a commit.

## Integrity, correctness and provenance

Archive `stage202-review-54693.tar.gz`, SHA256
`43723d1e31a3bba68cb02573cea9a54927f1f08fd2adae65c8c1b1a28506d8ad`.
All 4,616 inventoried files pass size/hash checks; no missing or extra files.
All 17 runs and 37 recorded gates pass, including the combined checked-plus-oracle
run, exact baseline/candidate/checked/observer checkpoint comparisons, identical
printed adaptive workload, checkpoint reload and post-restart remapping. Legacy
passes the existing fixed numerical screen; it is not bitwise equal to blocks.

Independent reconstruction from raw logs reproduces the saved analyses and timing
medians. Re-reading all six observer runs reproduces every counter and phase
summary, including per-rank/per-step coverage and the common clock ledger.
All 17 sets of 83 rank placements match the independently reconstructed map:
exclusive bb13/nk, 42/41 ranks across sockets, using the prescribed cache-spread
mapping. All four PMU runs have 100% running/enabled time, no quality flags and no
observer-calibration flags. Printed observer medians differ from bare medians by
less than 1%; major dynamics envelopes remain within the calibration thresholds.

All six recorded source inventories and build-log hashes match the frozen kit.
Every run's executable identity matches its build record; build/test hosts match,
contraction is off, and the recorded 23 common physics objects pass the build
protocol's identity gate. The archive intentionally excludes checkpoint and
executable payloads: numerical equality relies on the archived comparator reports,
not a fresh payload comparison during this review.

## Bare timing and remaining gap

| Version | Repetition 1 | Repetition 2 | Median seconds/ordinary step |
|---|---:|---:|---:|
| Stage201 baseline | 1.050 | 1.060 | 1.0550 |
| Stage202 prototype | 0.923 | 0.927 | 0.9250 |
| Legacy Domain | 0.540 | 0.543 | 0.5415 |

Each run completes the same 12-step window; checkpoint steps are excluded from
ordinary-step medians. Paired recoveries are **0.127 and 0.133 s/step**. Even the
smallest baseline median minus the largest candidate median is 0.123 s/step.
Warmups support the same separation. Two repeats establish this prototype screen,
not precise statistical confidence or a claim about other rank counts/problems.

The blocks/legacy ratio falls from **1.948 to 1.708**. The absolute excess falls
from 0.5135 to **0.3835 s/step**. Thus 25.3% of the current measured gap is removed;
71% extra time relative to legacy remains. These figures use the same-job baseline
and legacy, not the historical 2.25x ratio under different measurements.

## Where the improvement occurs

These are ordinary-step **rank-mean coarse phase times**, not the printed median
statistic above. Phase savings include changes in MPI progress and waiting.

| Phase | Baseline seconds | Candidate seconds | Saving | Retired-instruction reduction, two pairs |
|---|---:|---:|---:|---:|
| Dynamics wavelets | 0.283765 | 0.244529 | 0.039236 | 14.11–14.86% |
| Remaining dynamics | 0.574779 | 0.485254 | 0.089525 | 17.55–18.03% |
| Adaptation/reconstruction | 0.149223 | 0.149659 | -0.000435 | -1.92 to -0.63% |
| Whole ordinary step | 1.026284 | 0.897711 | 0.128573 | 13.45–14.09% |

The coarse total recovery, 0.1286 s, closely matches the 0.130 s bare recovery.
The two dynamics phases account for effectively all of it: approximately 30% in
wavelets and 70% in other dynamics. Adaptation is unchanged at this resolution.

Whole-step cycles fall 12.66–13.23%. Wavelet CPI is essentially unchanged
(candidate/baseline 0.995–0.998); remaining-dynamics CPI is about 2.4% higher.
The benefit therefore accompanies substantially fewer retired instructions,
consistent with removing repeated record resolution and inverse staging work.
It does not depend on an improvement in aggregate CPI. User-space MPI execution
and spinning are included in these counters: the instruction reductions cannot
all be labelled deleted solver instructions, and reduced waiting cannot be added
to producer savings as a separate benefit. There is no feature ablation here,
so exact savings cannot be assigned separately to record batching, direct local
copies and alias decoding.

No swap-in/out occurs during any bare run. Available node memory stays above
1,490 GiB. Summed per-rank peak RSS stays about 21.75–21.78 GiB for blocks and
11.07–11.08 GiB for legacy. These peaks are not synchronized live usage. Stage202
has improved execution time without materially reducing the memory footprint.

## Decision and next priority

Commit the tested combined prototype and this evidence. The predeclared threshold
is met with clean numerical and observer gates; another broad timing/scaling
campaign is not needed to accept this particular change. Enabling the combined
path by default is a separate source change; the current prototype remains opt-in.

The residual gap still warrants substantial work. Comparing Stage202 candidate
coarse phases with the previous Stage201 legacy coarse controls gives provisional
remaining excesses of about 0.175 s in wavelets, 0.110 s in other dynamics and
0.090 s in adaptation. Legacy bare timing is stable between campaigns (0.5425
versus 0.5415 s), but these are cross-campaign phase comparisons, not a new matched
legacy profile. Use them to prioritize the surviving inverse staging/publication
and representation-maintenance paths, not as guaranteed recoverable amounts.
Stage202 retains the inverse workspace, consumer alias epochs, scaffold state
and several full-record interfaces. Any further removal must identify the exact
remaining copies and ownership dependencies and retain the substantial-gain gate.

Machine-readable review: `parallel-blocks-cluster-202-54693.json`.
