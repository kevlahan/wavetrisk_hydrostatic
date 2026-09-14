# Stage 193 cluster review: job 54646

Accept the Stage 193 scalar inverse column kernel for this tested configuration
and commit the implementation and evidence. It produces a consistent timestep
gain with unchanged inverse transfer volume. Keep the Stage 192 sparse scatter
experiment disabled. The numerical source and defaults have not been changed by
this review; Stage 193 remains opt-in until its default is explicitly promoted.

## Provenance and validation

The supplied `stage193-review-54646.tar.gz` has SHA-256
`692139246cc7f11676a4b9fa554a6c06ebd7e8c0465c734668847e77963c8586`.
All 4,683 inventory entries passed byte-count and hash checks. Build provenance
matches the frozen Stage 193 kit for all five variants, including the original
legacy sources, common I/O patch and separately instrumented legacy copy.

- Accepted baseline: Stage 191 `30c831149653205210531427f81b184663a3acbf`.
- Candidate patch: `e00d6dbb6d010c34c1f0fcc56b0d127d80551cda7d4174e0078cc7527d0cfc33`.
- J5/J7, RK4, 83 ranks on bb13, GNU Fortran 13.2.0 / MVAPICH environment.
- All optimized compiler commands use O2, `-march=native`, contraction off.
  Checked commands use O0, runtime checks and contraction off.
- The builder records 23 identical external physics objects across variants;
  solver runtime checks do not imply those external objects were checked.
- All ranks retained the same distinct single-core bindings across all runs.
  No swap-in/out was recorded. Maximum rank RSS was approximately 379–382 MiB
  for blocks and 166–167 MiB for legacy.

All 27 cases completed their 12-step window. Each read checkpoint 3, wrote and
reloaded checkpoint 4, and continued through subsequent vertical remapping.
All 62 gates passed: 29 exact checkpoint reports (including reference
self-checks), 25 workload comparisons and eight fixed legacy numerical screens.
The raw logs confirm the checked-on dynamics, adaptation and remap oracles were
enabled; the new scalar kernel also remained enabled during that run. All printed
step/time/dt/Jmax/dof sequences matched the baseline. Checked-off took 56.45 s;
checked-on took 185.41 s, both excluded from performance medians.

Maximum reported differences against legacy were:

| Field | Maximum absolute difference | Existing limit |
|---|---:|---:|
| Velocity | 3.72324e-12 | 1e-10 |
| Mass | 2.50111e-12 | 1e-10 |
| Temperature | 7.85803e-10 | 1e-7 |

Legacy equality is the existing fixed numerical screen, not bitwise equality.
Checkpoint payloads and binaries are excluded from this review archive. The
review verifies the recorded comparisons, hashes, logs and provenance; it does
not independently rerun absent checkpoints or executables.

## Repeated uninstrumented timings

Six rotating baseline/candidate/legacy triples, after separate warmups:

| Measure | Stage 191 baseline | Stage 193 | Legacy |
|---|---:|---:|---:|
| Median ordinary-step time | 1.340 s | 1.280 s | 0.556 s |
| Median elapsed run time | 31.7442 s | 31.0899 s | 22.4788 s |
| Ordinary-step ratio to legacy | 2.4101× | 2.3022× | 1× |

Stage 193 reduced ordinary-step time by **4.4776%** and elapsed run time by
**2.0612%**. It won all six paired elapsed comparisons by 0.387–0.742 s and
all six ordinary-step comparisons by 0.05–0.06 s. The ranges did not overlap:
baseline ordinary medians were 1.34–1.35 s; candidate medians were 1.28–1.29 s.
The elapsed ranges were 31.624–31.839 s and 31.024–31.237 s respectively.
Ordinary CPU times are printed at finite precision, so percentage precision
beyond approximately 4.5% should not be overinterpreted.

The same candidate binary with the new kernel disabled returned to 1.35 s per
ordinary step and 31.7677 s elapsed. That single control is diagnostic and is not
included in either timing median. It supports attributing the repeated gain to
the new kernel.

This removes approximately 7.65% of the excess ordinary-step time above legacy
in this campaign. The remaining factor of 2.30 is substantial. No claim is made
for other CPU architectures, rank counts, longer integrations or multiple nodes.

## Computational and communication evidence

Matched whole-window rank-mean profile times include the checkpoint step. Each
implementation executed 60 complete WT_after_step transactions per rank. These
are separate instrumented runs, not the repeated timing runs above.

| Matched scope | Baseline | Candidate | Reduction |
|---|---:|---:|---:|
| Combined scalar, inner-vector and outer-vector inverse kernels | 0.759201 s | 0.361888 s | 52.33% |
| Complete native inverse parent | 3.586475 s | 2.892992 s | 19.34% |
| Complete wavelet transaction | 5.588820 s | 4.903210 s | 12.27% |
| Inverse MPI wait self time | 1.415956 s | 1.149575 s | 18.81% |

These nested scopes must not be added together. Inspection during Stage 194
established that the historical "scalar/inner" timer also encloses the
outer-vector kernel (detail region 21). The measured reduction is unchanged,
but the remaining 0.362 s cannot be assigned to scalar or inner-vector work
alone. Stage 194 adds separate scalar and inner-vector children and accounts
for the existing outer-vector child and parent overhead. The preserved raw
campaign records retain their original labels. Only one profile per
implementation was collected.

The complete legacy wavelet transaction took 1.28804 s, leaving the candidate
at **3.8067× legacy** for this matched scope. Comparing bare inverse timers to
legacy would mix differing call counts and boundary scopes.

The scalar inverse produced 92,739,840 field samples using 1,545,664 compiled
stencils. Field-address resolutions were 10,433,232 instead of a
reference-equivalent 625,993,920. Geometry reads were 10,819,648 instead of
649,178,880. Both preparation counts fell 98.33%; arithmetic samples were not
removed. Reference-equivalent counts follow the retained scalar kernel and are
verified against actual accessor calls by the focused Fortran test.

Inverse gather, alias and scatter counts were exactly unchanged: a combined
42,263,016,000 remote bytes and 496,004 remote messages across the profiled
window. Self copies are excluded. The reduction in MPI wait despite unchanged
traffic is consistent with less rank-local preparation reducing waiting at
communication boundaries; it does not demonstrate reduced network bandwidth
cost. Packing itself did not improve (0.30866 to 0.32189 s in these single
profiles). This evidence supports continuing the computational priority.

Memory usage did not materially improve. Median sums of rank peak RSS were
21.6659 GiB baseline, 21.6739 GiB candidate and 11.0153 GiB legacy. Those sums
are not simultaneous node memory measurements.

## Recommended follow-up

1. Commit the validated Stage 193 kernel and its test/review records, preserving
   a disabled control. Promote it to the production default with the selector
   tests updated for the new default. Keep Stage 192 sparse scatter opt-in.
2. Continue removing repeated computational preparation. Inner-vector inverse
   geometry is a concrete candidate, but first separate scalar and inner-vector
   kernel times so the available gain can be assessed. Reuse the same direct
   baseline/legacy comparison and fixed validation gates.
3. Retain the complete wavelet transaction as a matched legacy boundary. The
   remaining 3.81× wavelet gap and roughly doubled block memory footprint remain
   unresolved; the existing differently nested component timers do not justify
   attributing all of that excess to one routine or to message volume.

Machine-readable evidence is in `parallel-blocks-cluster-193-54646-analysis.json`
and the preserved `parallel-blocks-cluster-193-54646-campaign.json`. Reproduce the
main audit with:

```sh
python3 test/parallel_block_profile/analyze_stage193_campaign.py /private/tmp/stage193-review-54646
```

The additional all-variant identity check in the saved analysis compares the
archived build records to the frozen local Stage 193 kit manifest.
