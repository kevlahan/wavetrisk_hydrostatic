# Stage203 cluster review — job 54695

Stage203 passes the correctness and measurement gates, but **fails the predeclared
performance acceptance criterion**. Preserve it as an experimental result; do not
promote it as the accepted next performance baseline or lower the threshold after
seeing the results. Stage202-on remains the accepted baseline. No commit or
production-default change is made by this review.

## Validation

Archive `stage203-review-54695.tar.gz`, SHA256
`0992d1216242ae4fe4da1d23b5190fbf0835026fd3529df5cbe101eb67c7532c`.
All 5,034 inventoried files pass size/hash checks, with no missing/extra files.
Seven recorded source inventories and build-log hashes match the frozen kit;
recorded executable hashes match every run. Baseline builds use Stage202-on and
Stage203-off; candidate builds use both on. Contraction is off, J5 parameters
and 23 common physics objects are recorded, and compilation occurred on bb13.

All 18 runs complete the same 12-step workload, reload checkpoint 4 after seed 3,
and remap after reload. Independent reconstruction reproduces the raw-log analyses,
all workload gates, 83-rank placement for every run (42/41 across sockets), all
seven observer datasets, phase-balance summaries, inverse ledgers and final
summary. All 39 archived gates pass. This includes runtime-checked oracle-off/on,
exact semantic checkpoint equality across blocks/checked/observers, and the fixed
legacy tolerance screen. Legacy is not bitwise equal: maximum atmospheric
velocity, mass and temperature differences are respectively 3.72e-12, 2.50e-12
and 7.86e-10, inside the unchanged limits.

Numerical comparator reports were re-evaluated against the policy; checkpoint and
binary payloads are intentionally absent from this archive, so this review does
not rerun payload comparisons or independently hash the actual executables.
All four PMU controls have running/enabled fraction 1.0 and no quality flags.
Observer calibration passes: printed block medians differ from bare controls by
at most 0.65%; legacy coarse differs by 1.67%. There is one legacy coarse control,
not repeated independent calibration of each legacy phase.

## Performance decision

| Implementation | First bare run | Second bare run | Median s/ordinary step |
|---|---:|---:|---:|
| Stage202-on baseline | 0.940 | 0.924 | 0.9320 |
| Stage203 combined | 0.878 | 0.877 | 0.8775 |
| Legacy Domain | 0.541 | 0.539 | 0.5400 |

Stage203 saves **0.0545 s/step (5.85%)** on these medians. The two reversed paired
savings are **0.062 and 0.047 s**, both below the required **0.10 s**. Even the most
favourable pairing is only 0.063 s. Warmup/observer results support a real benefit;
this is not evidence of zero improvement, but it is insufficient for the agreed
design goal. Two repetitions are a prototype screen, not a precise confidence
interval or a guarantee for other problems/rank counts.

The same-job gap falls from 0.3920 to **0.3375 s/step**, removing **13.9%** of it;
the ratio falls from 1.726 to **1.625 times legacy**. The prior 0.3835 s gap used
Stage202 job 54693 timings; using that older baseline to score this trial would
mix campaigns. The historical 2.25x ratio is not the current measured ratio.

## Where time was saved and where the excess remains

These are same-job coarse **rank-mean phase seconds per ordinary step**, distinct
from the printed bare medians above. Phase intervals include MPI progress/waiting.

| Phase | Stage202 | Stage203 | Legacy | Saved | Remaining excess |
|---|---:|---:|---:|---:|---:|
| Dynamics wavelets | 0.244274 | 0.214171 | 0.070201 | 0.030103 | 0.143969 |
| Other dynamics | 0.483282 | 0.467907 | 0.375764 | 0.015375 | 0.092143 |
| Adaptation/reconstruction | 0.148797 | 0.148776 | 0.060332 | 0.000021 | 0.088444 |
| Whole step | 0.894845 | 0.849300 | 0.527751 | 0.045545 | 0.321549 |

Other phases contribute about -0.0030 s net to the residual. Thus the three main
phases account for the coarse excess, but this is a phase budget, not proof that
all of that time is redundant/removable work. Wavelets remain about 3.05x legacy;
adaptation remains about 2.47x. Approximately two thirds of this trial's coarse
saving is in wavelets and one third in other dynamics. Adaptation did not improve.

Wavelet retired instructions fall 12.67–13.15%, cycles 12.51–12.93%, with CPI almost
unchanged. Other-dynamics instructions fall 4.25–5.12%, cycles 3.03–3.78%, while CPI
worsens about 1.3–1.4%. Whole-step instructions fall 5.73–6.20%. This supports
reduced execution work, rather than improved aggregate CPI, as the mechanism.
Counters include user-space MPI progress/spinning; they cannot assign every saved
instruction to solver computation. No feature ablation separates the two changes,
and no per-operation timing identifies the cost of surviving representation work.

The coarse saving (0.0455 s) and bare-median saving (0.0545 s) use different
statistics/runs; the difference is not another separately attributable saving.

## Why the larger expectation was not met

The implementation eliminated the intended intermediate operations, but the
opportunity envelope was larger than their measured net benefit. It was not a
measurement of removable work. The combined test cannot say precisely how much
the deleted passes cost before replacement or how much the new execution owner,
compiled kernels and final installs offset their removal.

The ledger supplies useful constraints, identical across all three candidate
observers (whole 12-step run, including checkpoint steps):

- 48 inverse transactions per rank; 488 replaced gather and 488 replaced scatter
  calls per rank. The intended resident path is active.
- 92,739,840 scalar outputs and 11,592,480 inner cells across ranks.
- **534,760,800 logical interior scalar values installed at exit**, excluding
  boundary installs, entry imports and Domain publication. This is about 4.28 GB
  of double payload over the whole run, not measured DRAM traffic, unique data,
  network volume or recoverable seconds. Full installs remain, including the
  second install for the provisional final-boundary contract.
- Zero sparse canonical-restoration values on this fixture. The nonempty restore
  operation is therefore not a plausible large cost here; checking/indexing and
  fixed scaffold restoration can still cost time. Local overlap unit coverage is
  separate from this cluster fixture.
- Scalar/inner work has a maximum/mean ratio of 2.64, with one zero-work rank;
  installed interior values have a ratio of 2.25. These are unequal work counts,
  not measured kernel-time imbalance or comparisons with the old owner map.

The wavelet phase's mean maximum-minus-mean rank time is 6.20 ms before and
5.99 ms after; other dynamics is 8.37 versus 8.68 ms. There is no large new
phase-envelope imbalance explaining the shortfall. However, collectives can
hide producer imbalance inside waiting, so these observations do not exonerate
ownership imbalance. Phase maxima must not be summed as a critical path.

No swap occurs in any bare run; available node memory stays above 1,491 GiB.
Summed per-rank peak RSS is about 21.75–21.79 GiB for blocks versus 11.07–11.08 GiB
for legacy. These unsynchronized peaks show no material footprint reduction;
they do not establish cache behaviour. Whole-run major faults include startup/I/O
and cannot be interpreted as ordinary-step stalls.

## Recommended next action

Preserve this prototype and evidence without declaring the stage accepted. Before
another optimization, close the remaining operation-level attribution gap:

1. Measure mutually exclusive inverse entry/import, scalar/inner/outer kernels,
   alias local preparation/installation, MPI completion, final block installs and
   Domain publication. Use matching legacy kernel/boundary scopes. Include the
   forward transform separately so the entire wavelet excess is reconciled.
2. Measure surviving scalar adapter/ghost/transport work versus restriction and
   divergence arithmetic, and split adaptation into topology construction,
   reconstruction, representation rebuild and publication with matched legacy
   scopes. Adaptation's 0.088 s excess is currently untouched.
3. Calibrate focused observers against the existing coarse envelopes, report
   unassigned time, and retain rank-local work/time pairing. Attribute waiting to
   dependencies without adding it to producer savings. Reuse the fixed placement
   and short workload; do not repeat the broad scaling campaign or chase precise
   timing statistics before the causal accounting is complete.

The entry/exit and duplicate-publication paths are candidates for investigation,
not proven next optimizations. Do not promise another 0.10 s from their value
counts. A design should identify enough measured, removable work to justify that
threshold before another implementation campaign. The independent fields/inverse
selectors are available for a small ablation if the exclusive measurements cannot
separate their effects; no ablation is present in this archive.

Machine-readable audit: `parallel-blocks-cluster-203-54695.json`.
