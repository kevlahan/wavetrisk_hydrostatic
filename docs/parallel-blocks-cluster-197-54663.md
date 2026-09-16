# Stage 197 cluster review: job 54663

Current direction: the user has required diagnosis of the main legacy/block cost causes before any further optimization. The [diagnosis objective and evidence requirements](parallel-blocks-performance-gap-diagnosis.md) supersede optimization recommendations below. Existing results are retained as historical evidence.

Stage 197 passes its correctness gates and removes the intended ownership work.
The observed whole-step gain is small: 1.265 to 1.250 seconds (1.19%), with the
candidate still 2.246 times the legacy step. This validates the limited change;
it does not establish the source of the complete remaining performance gap.

Evidence: `stage197-review-54663.tar.gz`, SHA-256
`9cf514afb580bfce5d2b8df08bffb4cbdc855cd61c72c0bb6dcdff16a8a931b0`.
The archive's 3,137 inventory entries match their sizes and hashes. Derived
values and verification status are in [the review record](parallel-blocks-cluster-197-54663.json).

## Acceptance and execution conditions

All 14 executions completed, all 34 gates passed, and every execution took the
same 12 steps, reloaded checkpoint 4 and exercised post-restart remapping.
Stored exact checkpoint reports satisfy the original acceptance function;
the fixed numerical legacy screen also passes. Bitwise block/legacy equivalence
is not claimed. Checkpoints are intentionally absent from review archives, so
this review audits their recorded comparison reports rather than rereading the
checkpoint payloads.

Build records match the frozen baseline/candidate and legacy source inventories,
patch identity and six recorded build-log hashes. Executable hashes in the run
records match their build records. Builds use J5, RK4, contraction off and the
same 23 external physics objects; the optimized block compile lines include
`-march=native`. All runs have the same rank-to-core mapping on bb13, 83 distinct
single-core affinities across both sockets (42/41), and zero solver exit status.
Node swap-in/out deltas are zero throughout.

The sum of rank peak RSS is about 21.66 GiB baseline, 21.70 GiB candidate and
11.02 GiB legacy. These sums are not simultaneous memory measurements. Stage 197
has not meaningfully reduced peak memory, but node swapping does not explain
this timing gap.

All three raw operation/workflow profiles reparse successfully, conserving their
timing envelopes and matching recorded results. Work counts match exactly;
floating timing aggregates agree to 1e-12 relative/absolute tolerance, allowing
rounding differences between Python summation implementations. This tolerance
applies only to recomputed profiling summaries, not numerical checkpoint gates.

## Bare timings

Warmups, checked runs and profiles are excluded. The entries below are each
run's printed ordinary-step median, followed by the median across the two runs.

| Implementation | Round 1 | Round 2 | Median step | Median elapsed |
|---|---:|---:|---:|---:|
| Stage 196 | 1.270 s | 1.260 s | 1.265 s | 30.998 s |
| Stage 197 | 1.250 s | 1.250 s | 1.250 s | 30.764 s |
| Legacy | 0.558 s | 0.555 s | 0.5565 s | 22.674 s |

The candidate is faster in both timing rounds, including the reversed second
round. Observed reductions are 1.19% in step time and 0.76% in total elapsed time.
Two samples and rounded printed times do not support a precise statistical
speedup claim. The bare candidate-minus-legacy gap is still 0.6935 seconds.
More repetitions could refine this small gain, but are not needed to decide
that this particular ownership change does not remove most of the gap.

## What the change actually removed

Per rank per ordinary step, the candidate records four avoided scalar imports
and four avoided uncompressed publications; the control records zero. Full-import
buffer initialization falls from 1,354,423 to 406,327 doubles (70%), and full-import
patch installation visits halve. Local coefficient initialization adds 414,638
doubles. These are distinct logical work counters, not a net DRAM-byte estimate.

Mean time across ranks in the complete WT context falls from 0.354893 to
0.338110 seconds, a 16.784 ms reduction (4.73%). The selected import/writeback/
publication groups, net of local initialization, account for about 15.438 ms of
that reduction. Outside-WT dynamics is essentially unchanged: 0.668097 versus
0.669342 seconds. The profile location of the saving is consistent with the
15 ms observed bare-step reduction.

The larger inverse exchange sequence is unchanged. Both block profiles contain
176,956 dependency records. Every ordinary dependency key has identical call,
request-bearing call, request and sent/received value counts between variants.
Mean ordinary completions per rank per step are about 56.55 gathers, 56.55
scatters and 64.18 alias/boundary refreshes, including empty completions.
Candidate completion intervals sum to about 17.71, 31.09 and 38.27 ms respectively.
Planning/packing/install still costs 55.59 ms in the inverse groups, and the
posting/self-copy/wait group still costs 90.11 ms. Rank-local wait durations
include arrival imbalance/progress; they are not network service times or a
promise of recoverable savings. Stage 197 removed an outer import/publication
round trip, not these repeated inverse phases.

## What the computational measurements establish

Counters 4, 5 and 12–15 now report actual work. They do not show an excess of
ordinary generic Domain boundary traversal in blocks: outside-WT boundary pack
candidates are 3.825 million versus legacy's 6.896 million, and buffer zeroing is
203.620 million versus 290.327 million doubles per rank per ordinary step.
Blocks also have separate native transport and ownership paths; lower counts
in the shared Domain routine do not imply lower total boundary work.

| Selected outside-WT count, per rank per ordinary step | Stage 197 | Legacy |
|---|---:|---:|
| Direct scalar node-field samples | 165,855 | 159,684 |
| Temperature direct-flux node evaluations | 60,679 | 79,842 |
| Restricted scalar flux edges | 7,588 | 45,527 |
| Divergence samples | 60,679 | 139,767 |
| Divergence evaluations | 35,559 | 83,383 |

These have explicitly limited source scopes, particularly legacy's regular
step1 branches and its fused producers. They are not an exhaustive count of
all special boundary work. They do not support attributing the whole slowdown
to repeating more scalar arithmetic or reproducing the complete legacy timestep.
The direct-flux/traversal scopes take 55.23 ms in blocks and 24.95 ms in legacy,
but the latter includes shared step1 work. This is not an arithmetic-only
per-evaluation ratio.

Source inspection identifies a concrete redundant preparation pattern:
`compute_block_scalar_direct_flux` requests a complete 50-slot record through
`scalar_read_range` before testing the mass shortcut or physical vertical level.
The production mass shortcut uses one of those values; the temperature flux
uses three edge lengths and three physics flux values. With compact storage,
this 50-slot mixed geometry/field range takes the generic reconstruction path,
which dispatches individual slot reads and resolves shared geometry addresses.
The measured sample count corresponds to about 8.29 million requested record
slots per rank per ordinary step at this call site alone. This is source-level
logical work, not measured instructions or memory traffic. Its actual time
contribution still needs isolation.

Relevant source: `src/parallel_block_mpi.f90`, `compute_block_scalar_direct_flux`;
`src/parallel_block_scalar_storage.f90`, `scalar_read_range`, `read_value` and
`shared_address`. There are nine `scalar_read_range` call sites in the MPI module;
those in production hot paths need a caller/field-use audit rather than an
assumption that they all have equal cost. The compressed storage saves capacity,
but recreating its old full-record interface can add execution overhead.

## Profiling limit and next decision

The profile step medians are 1.270 s baseline, 1.270 s candidate and 0.605 s legacy:
0.40%, 1.60% and 8.72% above their respective bare medians. The new per-leaf legacy
counter hooks can perturb the execution and waiting distribution. Consequently,
do not turn this campaign's profile differences into precise percentages of
the bare gap. For example, the raw profile WT, outside-WT and adaptation excesses
are about 256.8, 252.0 and 114.9 ms, but the total profile step excess is about
621.6 ms versus 693.5 ms bare. Those are different measurements.

Recommendation: commit the validated Stage 197 candidate, its measurements and
submission-guard repair, describing the small gain as observed and retaining
the reversible selector. Do not spend another long campaign precisely estimating
this one-percent result before investigating the larger mechanism.

For the next targeted stage, audit and replace full compatibility-record
reconstruction in production scalar kernels with the fields they actually use,
resolve horizontal geometry once per column, and consume contiguous field arrays.
Keep independent oracle records and checkpoint/workload gates. Measure record
reads, address resolution and copying across their production callers, not just
another flux arithmetic subexpression. This class of overhead is a stronger
computational hypothesis than excess temperature/divergence evaluation counts;
it is not yet a demonstrated explanation of the entire gap.

At the same time, batch legacy work-count reporting outside leaf loops, then
check profile overhead again. Retain the short direct previous-stage/candidate/legacy
comparison and exact gates; add timing rounds only if a precise small-gain claim
becomes necessary. No further production change or commit was made during this
review.
