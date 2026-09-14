# Stage 190 cluster campaign 54553

Accept Stage 190 for the tested J5/J7 RK4 configuration on nk. All 36 recorded
checks pass, and all three candidate elapsed times are below all three baseline
times. This is a repeatable improvement within this short campaign, not a
universal speedup or a multi-node result.

## Correctness and provenance

All 16 runs complete 12 timesteps, reload checkpoints 3 and 4, and exercise a
vertical remap following checkpoint 4. Input settings are time_end=0.315,
dt_write=0.01, max_level=7, 30 atmospheric layers and 10 soil layers; compilation
uses param_J5 and RK4.

The 36 checks comprise 18 exact checkpoint reports (including three reference
self-checks), 14 workload comparisons, and four fixed legacy numerical screens.
All exact reports show zero differences in semantic fields, thresholds, headers
and topology. These include checked oracle-on/off, checked/optimized,
Stage 189/190, repeated runs and profiled runs. Workload comparisons require the
same printed time, dt, Jmax and degree-of-freedom sequence.

Legacy Domain with the common checkpoint I/O fix passes the existing fixed
numerical screen, but is not bitwise equal to the block code. Maximum absolute
field differences are 3.72324e-12 for velocity, 2.50111e-12 for mass and
7.85803e-10 for the stored temperature field, below respective limits 1e-10,
1e-10 and 1e-7. Soil fields, headers, topology and thresholds match exactly.
Raw checkpoints were not supplied; this review verifies the supplied reports
rather than repeating their checkpoint comparisons.

The embedded build record equals the separately supplied build record. All
candidate source hashes match the local candidate, including symlink identities.
The baseline changed translation unit matches commit a0c75adb4dc1741f56ecb572ab1809e8a70a32e4.
Every run binary hash matches its build role. All three supplied build-log hashes
match the manifest. The optimized baseline and candidate compilation lines use
-O2 -march=native -ffp-contract=off; checked uses -O0, runtime checks and
-ffp-contract=off. Builds ran on bb13, AMD EPYC 9754. The build record confirms
23 shared physics object identities and four passing builds. The legacy build
log was not supplied, so its compilation provenance relies on the manifest and
build harness rather than a separate inspection of expanded compile lines.

## Timings

Warmups, correctness runs and detailed profiles are excluded. Three timing
triples rotate baseline/candidate/legacy, candidate/legacy/baseline, and
legacy/baseline/candidate in one allocation.

| Metric | Stage 189 | Stage 190 | Legacy Domain + I/O fix |
| --- | ---: | ---: | ---: |
| Median elapsed | 33.64821 s | 32.92663 s | 22.53499 s |
| Median of ordinary-step medians | 1.49 s | 1.43 s | 0.558 s |

Stage 190 reduces elapsed time by 2.14449% and ordinary-step time by 4.02685%.
Candidate elapsed times range from 32.77928 to 32.94249 seconds; baseline ranges
from 33.64119 to 33.75145 seconds. Stage 190 remains 1.46113 times legacy elapsed
time and 2.56272 times legacy ordinary-step time. This campaign establishes the
remaining gap directly on the same workload and allocation.

## Profiles and memory

The scalar replay region's exclusive rank-mean wall time, summed across the two
reporting windows, decreases from 1.519247 to 1.169265 seconds (23.04%). Its
corresponding CPU time decreases from 1.507280 to 1.158794 seconds. Scalar
transport MPI time also decreases, but this is not proof of reduced network
work: waiting can reflect changed arrival times and load imbalance. No MPI
protocol was changed by Stage 190. Profile sums are not a decomposition of the
solver critical path.

The new closure reports 252344 edges, 4037504 column locations, 1009376 cached
reuses and 18240 fallback layers. Both windows have 16 locations and four reuses
per edge, confirming use of the intended path.

Every run has 83 child-resource records collected with os.wait4. Core assignments
are identical across runs and distinct per rank: 42 cores on socket 0 and 41 on
socket 1, all on bb13. Both NUMA nodes are allowed for memory. Node snapshots
show zero swap-in/out increments for every run and abundant available memory;
there is no evidence of the local machine's memory-pressure problem here.

Timing-run maximum rank RSS is approximately 397--399 MiB for both block
versions and 166--168 MiB for legacy. Summed per-rank peaks are approximately
23.0 GiB for block and 11.0 GiB for legacy. These are sums of separate peaks,
not simultaneous node peaks; they exclude the waiting Python parents. Stage 190
therefore improves computation without materially reducing the measured memory
footprint. Major page faults are recorded, so this is not a zero-page-fault claim.

## Decision and next investigation

Recommend committing Stage 190 implementation, tests, runner repair and these
validation records. No additional cluster run is required to accept this bounded
change. Keep contraction off and compile native binaries on their target node.

The remaining performance investigation should prioritize the boundary and
inverse communication paths, scalar replay, and repeated storage/copy work.
The candidate profile's Domain boundary MPI wait and inverse MPI wait sums are
1.879965 and 1.421739 seconds respectively; these require rank-level imbalance
analysis before attributing them to communication itself. The legacy memory gap
also remains. Checkpoint/bootstrap overhead is substantial in this short test;
ordinary-step timings are the better indicator of its ongoing compute gap.

## Evidence

- parallel-blocks-cluster-190-54553-campaign.json: byte-preserved supplied results.
- parallel-blocks-cluster-190-54553-analysis.json: identities, checks, timings,
  profiles, placement and per-run memory summaries.
- Raw supplied files preserved in /private/tmp/wavetrisk-stage190-cluster-54553/.

No solver source or compiler settings changed during this review. Stage 190 is
not committed by this analysis.
