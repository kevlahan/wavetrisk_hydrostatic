# Stage 192 cluster review — job 54581

Stage 192 passes this campaign's correctness requirements, but does not
establish a useful performance improvement. Keep the accepted Stage 191
configuration as the production default. Preserve the experiment, validation
records and three-way profiling infrastructure; do not label the sparse route
as a demonstrated speedup. No production selector was changed during this review.

## Evidence verified

The supplied `stage192-review-54581.tar.gz` contains 4,683 inventoried files;
all sizes and SHA-256 hashes match. All 27 runs completed and all 62 reported
gates passed: 29 exact-checkpoint reports (including reference self-checks),
25 workload gates and eight unchanged legacy screens. Fresh checked-off and
checked-on runs passed. Every run completes 12 timesteps, loads checkpoint 3,
saves/reloads checkpoint 4 and remaps afterward.

Baseline is `30c831149653205210531427f81b184663a3acbf`. The candidate patch is
`f23cf5222e285981ffadad1440ec53b23e7d9a74d5eae2436ceed1914f2b6b4c`.
All five build source identities match the issued kit. Compiler-log hashes and
expanded commands confirm J5 O2/native/contraction-off optimized builds and the
O0 checked build. Build records contain the common 23 external physics objects.
The legacy numerical sources remain the pinned revision plus common I/O fix;
its instrumented copy passes exact comparison against uninstrumented legacy.

All runs use the same 83 nonoverlapping core bindings on bb13. Reparsed rank
resource files match the embedded records. There is no node swap activity.
Maximum rank RSS is essentially unchanged: Stage 191 379.04–382.27 MiB versus
Stage 192 379.65–382.12 MiB. Median sums of rank peaks are 21.670 versus 21.698 GiB
(legacy 11.006 GiB); these sums are not simultaneous node memory measurements.

The archive omits checkpoint payloads and executables. This review validates
supplied gate reports, raw logs, resource files, source/build identities and
profile calculations; it does not independently rerun the checkpoint comparison
against absent payloads.

## Timing

Six rotating timing triples, excluding warmups, validation, controls and profiles:

| Metric | Stage 191 | Stage 192 | Legacy |
|---|---:|---:|---:|
| Median ordinary-step time | 1.345 s | 1.340 s | 0.556 s |
| Median elapsed time | 31.832 s | 31.752 s | 22.532 s |
| Ordinary-step ratio to legacy | 2.419x | 2.410x | 1.000x |

Observed reductions are only 0.37% for ordinary steps and 0.25% elapsed. The
ordinary-step inputs are rounded to 0.01 s at this scale; baseline/candidate
values repeatedly overlap at 1.34–1.35 s. Candidate elapsed is lower in four of
six triples, with paired baseline-minus-candidate differences of -147, +266,
+46, +69, +100 and -31 ms. These results do not demonstrate a convincing,
repeatable gain. The disabled control records 1.35 s ordinary steps.

## What the profile resolves

The complete matched WT transaction still makes 60 calls per rank across 12
steps. Whole-window rank-mean inclusive wall time is 5.59833 s for Stage 191,
5.60518 s for Stage 192, and 1.28261 s for legacy. Stage 192 is 4.37x legacy in
this scope. The profile includes the checkpoint step and is not a headline
timing measurement.

Remote scatter payload falls from 25,348,300,800 to 23,565,148,800 bytes: 7.03%
less scatter and 4.22% less total inverse payload. Gather and alias bytes are
unchanged. Message counts remain exactly 496,004 overall, including 81,044
scatter messages; helper call counts are unchanged too.

Despite fewer bytes, scatter helper wall time rises from 0.77549 to 0.79219 s
(+2.15%). Across all inverse phases, packing self time rises from 0.3100948 to
0.332953 s (+7.37%), while explicit MPI wait falls only from 1.414772 to
1.409602 s. This is consistent with extra edge-record packing/dispatch work
absorbing the small transport saving. One profile per variant is insufficient
to prove a general hardware bottleneck or a repeatable component regression.

The experiment shows that this amount of payload pruning is insufficient. It
does not show that bandwidth is irrelevant. Prioritize keeping values on final
owners across dependent inverse phases and eliminating complete gather/alias/
scatter rounds, preserving alias epochs and outer-edge write order. Measure the
same complete WT transaction and ordinary-step medians against Stage 191 and
legacy. Avoid attributing the remaining gap by subtracting unmatched legacy and
block component labels.

## Records and reproduction

- `parallel-blocks-cluster-192-54581-campaign.json`: supplied results verbatim.
- `parallel-blocks-cluster-192-54581-analysis.json`: verified identities,
  correctness reports, resources, profiles and derived measurements.
- `test/parallel_block_profile/analyze_stage192_campaign.py`: reproducible audit.

After extracting the archive, run:

```sh
python3 test/parallel_block_profile/analyze_stage192_campaign.py /path/to/extracted-review
```

Commit recommendation: preserve the numerical experiment separately from the
cluster infrastructure and result records. Keep the sparse route opt-in if it
is retained in the main branch. Changing that default is a subsequent source
change, not something represented by job 54581's source hashes.
