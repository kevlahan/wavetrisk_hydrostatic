# Stage 200 cluster review — job 54689

Stage 200 passes all validation gates and identifies a substantial placement
sensitivity specific to blocks. It weakens explanations based primarily on
high task counts or worse relative computational balance. It does **not** yet
identify a source change with a defensible large recovery estimate. No new
production optimization or commit is part of this review.

## Evidence and validity

Input: `stage200-review-54689.tar.gz`, SHA256
`ac19ec48be8f1c1660e528198cfdadbebf94a98f2bce03877fd43dae9a5fd120`.
All 6,848 inventoried files pass size/hash checks. All **42 runs and 95 gates
pass**, including restart/remap coverage. Same-configuration repeats and
observer comparisons pass exact checkpoint gates; the cross-task-count and
cross-legacy comparisons pass the existing fixed numerical screen and workload
checks. The numerical screen is not a claim of bitwise legacy equivalence.
Same-rank one-/two-socket checkpoint comparisons are exact.

Independently reconstructed the placement plan from archived topology, verified
all 42 rank-resource sets, and reread all 16 raw observer reports using the
local analyzer. Integer counts, calls and volumes agree exactly; recomputed
aggregate times agree within 1e-10 seconds and individual rank times within
1e-12 seconds. Archived source provenance matches the prepared kit. Builds use
GCC 13.2.0, J5/J7, RK4, contraction off, common 23 physics objects, and optimized
`-O2 -march=native` on bb13. Each implementation uses the same bare binary in
all four cases. Numerical sources remain the frozen Stage 197 candidate and
pinned legacy with common checkpoint I/O fixes.

bb13 has two 128-core EPYC 9754 sockets, one hardware thread per core, two NUMA
nodes, and 32 last-level-cache groups. One-socket cases use nested core subsets
spread over all 16 cache groups on socket 0. The two-socket case uses 42/41
ranks spread across all 32 groups. Both implementations use identical explicit
rank affinity, verified after MPI initialization in observer runs, and local
first-touch allocation policy. The job owns the node exclusively.

## Bare production timings

These are medians of two ordinary-step medians, in seconds. Warmups, restart
steps and observers are excluded. Two repeats detect gross inconsistency;
they do not establish precise statistical confidence intervals.

| Tasks / sockets | Blocks repeats | Legacy repeats | Blocks median | Legacy median | Excess | Ratio |
|---|---:|---:|---:|---:|---:|---:|
| 21 / 1 | 2.77, 2.77 | 1.22, 1.20 | 2.770 | 1.210 | 1.560 | 2.29 |
| 42 / 1 | 1.75, 1.74 | 0.851, 0.852 | 1.745 | 0.8515 | 0.8935 | 2.05 |
| 83 / 1 | 1.19, 1.19 | 0.546, 0.546 | 1.190 | 0.546 | 0.644 | 2.18 |
| 83 / 2 | 1.06, 1.05 | 0.543, 0.545 | **1.055** | **0.544** | **0.511** | **1.94** |

Moving 83 ranks from one socket to two reduces block time by **0.135 s/step
(11.34%)**, versus just 0.002 s for legacy. The two block repeats show 0.130
and 0.140 s reductions. The corresponding excess shrinks by 0.133 s, or 20.7%.
This is a material configuration result, not a demonstrated code optimization.

Do not attribute the entire change from Stage 199's 1.25/0.560 seconds to this
experiment: its two-socket placement used only 12 cache groups, whereas Stage
200 spreads over 32. The within-Stage-200 one-/two-socket contrast is the
controlled comparison. Keep the new verified 83-task two-socket mapping for
subsequent code comparisons.

## What the experiment resolves

### Placement changes execution cost, not the counted work

For both implementations, all recorded per-rank work counts, MPI call counts
and MPI volumes are **identical** between the 83-task one-/two-socket probe
runs. Work ownership is unchanged. Yet coarse block phase times fall by about
37 ms in wavelets, 61 ms in remaining dynamics and 27 ms in adaptation.

In the probe, block whole-step rank-mean time falls by 133 ms: 101 ms outside
the selected MPI scopes and 32 ms inside them. Thus most of the observed
improvement occurs outside those MPI scopes. This favors an execution-throughput
or memory-hierarchy explanation over a reduction in messages or mathematical
work. Cache capacity/contention, memory resources, socket power/frequency and
MPI progress can all change with placement. No DRAM-bandwidth, cache-miss or
usable measured CPU-frequency data establish which mechanism dominates.
“Outside MPI” is elapsed time including unclassified work, not pure arithmetic
or hardware CPU execution time; unselected MPI can remain there.

No swap-in or swap-out counters increase. At 83 tasks, peak solver RSS is
roughly 389,000 KiB/rank for blocks versus 171,000–174,000 for legacy. Larger
RSS is consistent with a larger representation but does not measure active
working set, memory traffic or memory pressure.

### High task counts are not the primary explanation

The block/legacy ratio is already **2.29 at 21 tasks**. On one socket, 21-to-83
task speedups are 2.33 for blocks and 2.22 for legacy; neither scales ideally,
but blocks do not show worse aggregate strong scaling over this interval.
Communication/ownership overhead remains relevant; a slowdown caused mainly
by using too many ranks is not supported.

The mathematical-work counters are invariant across task counts within each
implementation. They do not indicate wholesale duplicate legacy numerical
steps in blocks. For example, blocks count about 5.04 million temperature
flux nodes versus 6.63 million in legacy, and 2.95 million divergence evaluations
versus 6.92 million. These counters have specific source scopes, not a complete
instruction census. They cannot rule out redundant representation maintenance,
record expansion, zeroing or compatibility copies. Boundary work grows with
task count in both implementations; lower logical byte/count totals do not
prove lower physical memory traffic.

### Blocks do not have worse relative balance in the main measured phases

At 83 tasks on two sockets, mean per-step rank-maximum/rank-mean time outside
selected MPI scopes is 1.46 for blocks versus 1.78 for legacy in wavelets, and
1.36 versus 1.62 in remaining dynamics. This agrees with the earlier distinction
between absolute spread and relative imbalance. Blocks can incur larger
absolute delays because their local execution takes longer; that does not
establish a uniquely poor ownership distribution.

The one-/two-socket block improvement occurs without improving these relative
balance ratios. Do not prioritize a load-balancer rewrite on this evidence.
Negative correlations between local time and MPI time are also not causal
producer/consumer traces: a nearly conserved phase duration itself induces
such a correlation.

## Low-overhead accounting of the remaining gap

Use coarse observers for the primary location accounting. Their printed-step
medians differ from bare timings by less than 1.7%. Probe overhead is asymmetric:
about 2.4–3.6% for blocks and 5.7–9.1% for legacy. At 83 tasks/two sockets it is
2.37% versus 6.62%. Probe phase envelopes flag changes in legacy's main phases;
there is no defensible universal overhead correction.

The coarse rank-mean signed phase differences are:

| Phase excess, seconds | 21 / 1 | 42 / 1 | 83 / 1 | 83 / 2 |
|---|---:|---:|---:|---:|
| Dynamics wavelets | 0.5457 | 0.3343 | 0.2571 | **0.2194** |
| Remaining dynamics | 0.8212 | 0.4189 | 0.2688 | **0.2060** |
| Adaptation/reconstruction | 0.1906 | 0.1322 | 0.1161 | **0.0897** |
| Other phases, net | 0.0034 | -0.0004 | -0.0029 | **-0.0033** |
| Sum | 1.5609 | 0.8849 | 0.6392 | **0.5118** |
| Bare median excess | 1.5600 | 0.8935 | 0.6440 | **0.5110** |

The two dynamics phases still account for about **83%** of the best-placement
coarse excess. The close sum establishes accounting by location, not that 83%
is avoidable. Rank means and bare printed-step medians are distinct statistics.
Do not sum independently selected rank maxima or overlapping parent/child scopes.

In the 83-task/two-socket probe, the whole-step difference is 0.493 s: 0.328 s
outside selected MPI scopes plus 0.165 s inside. Because the probe perturbs
legacy more, this split is diagnostic, not a precise decomposition of the
0.511 s bare gap or a sum of independent recoverable savings.

## Recommended next decision

The agreed causal goal is **not yet fully met**. This campaign gives a useful,
substantial configuration improvement and removes weak explanations, but it
does not justify promising a large production-code gain.

Retain the verified 83-task/two-socket/cache-spread placement. Focus the next
bounded comparison on the **complete wavelet and remaining-dynamics paths**:
extra instructions and representation/compatibility work versus cache/memory
stalls. Stage 199 already identified full scalar-record expansion and Domain
publication/import/writeback as concrete mechanisms; Stage 200 adds strong
placement sensitivity. Neither licenses assigning the entire deficit to one
getter, message exchange or bandwidth bottleneck.

A useful next diagnostic would collect matched, low-frequency call-stack and
hardware-counter evidence for instructions, cycles and available cache/stall
metrics on those paths, preserving rank identity and actual affinity. Establish
which counters the cluster exposes before depending on them, calibrate observer
impact on both implementations, and retain an uninstrumented timing reference.
No event name or availability is assumed here. A whole-path controlled
intervention may then test the dominant source-level mechanism. Keep necessary
numerical work and removable representation work separate, and avoid counting
both producer cost and consequent waiting as independent savings.

Reuse the established numerical/workload gates and focus on the best-placement
case; another full 42-run scaling sweep is unnecessary unless the next change
alters ownership or placement. Require a credible route to at least **0.10
s/step** recovery before selecting a production optimization. This is a proposal
screen, not a promised saving or statistical bound. No such code saving is
established by the current measurements.

Machine-readable evidence: `parallel-blocks-cluster-200-54689.json`.
