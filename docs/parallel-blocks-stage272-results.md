# Stage272: confirmed work concentration; phase timing is unqualified

Job 55471 completed all fourteen runs. Independent review verifies **5,257**
files and reproduces the science, placement, storage, checkpoint, VTK, raw
observer and timing checks. Archived identity, source and kit manifests match
the delivered kit byte for byte. All four builds share 23 physics objects.
Both 16/256-rank oracles retain 48 exact inverse-return proofs and 84 installations.
There is no production implementation to promote: Stage272 is diagnostic-only,
and production source remains exact accepted Stage266 `ff740744`.

The useful new evidence is **actual work distribution**, not a precise breakdown
of the slowdown. Retained Domain producer work does not spread to the extra
ranks, and those extra ranks receive less measured block work per rank as well.
The next proposed isolated intervention is to correct the front-loaded block
assignment, then measure it uninstrumented. Moving the complete Domain producer
and its dependent consumers to block owners remains the architectural objective.

## Untouched solver performance

Each entry summarizes two independent run medians over ordinary steps 2–7 and
10–12. Timing order brackets the 256 group with 160-rank runs.

| Ranks | Run medians, ms | Mean, ms |
|---|---|---:|
| 160 | 547.282, 547.733 | 547.507 |
| 256 | 588.490, 600.595 | 594.543 |

Adding 60% more ranks makes this fixed J5 trajectory **8.591% slower**, a
47.035 ms increase. Relative speedup is 0.9209× and efficiency is 57.56% versus
160 ranks. Both 256 samples are slower than either 160 sample. Stage271 likewise
found worsening scaling (6.436%); this is another campaign showing that direction,
not a precise cross-campaign estimate of how much the regression changed.
Main appears only as a 16-rank correctness anchor in this campaign; there is no
new legacy performance comparison.

## Calibration: retain both failures

| Ranks | Mean observer shift | Reason attribution is unqualified |
|---|---:|---|
| 160 | +9.516 ms (+1.738%) | Exceeds the preset 5.475 ms mean-shift bound |
| 256 | −0.491 ms | Untouched-run spread 12.105 ms exceeds the preset 11.891 ms spread bound |

At 160, both observer pairs are slower (+8.981 and +10.051 ms); the overhead
problem is repeatable within this allocation. At 256 the observer/control pairs
have opposite signs (+3.294 and −4.275 ms). That does not establish negligible
observer overhead. The 256 spread failure is small but must not be waived after
seeing the results. Cross-count phase timing remains **unqualified**.

Observer scope entries average 3,156.835 per rank/ordinary step at 160 and
2,420.826 at 256. These include the added Domain-boundary and native-mass wait
scopes. Their instrumentation is a plausible source of overhead, but this
campaign does not isolate the overhead of each hook. Do not subtract a global
9.516 ms correction from individual phase measurements or call such corrected
values production costs.

Exploratory records locate a strong waiting pattern: the 96 Domain-free ranks
spend about 152.8 ms per ordinary step in scalar-capture MPI in both observer runs,
versus about 2.4 ms for Domain owners. Source order places this collective after
the Domain compatibility producer. This supports investigating late arrival at
that handoff; it is **not calibrated evidence of 150 ms of recoverable step time**.
Inclusive timings overlap, local durations are unsynchronized, and MPI time
combines waiting and transfer. Exact timing attribution is deliberately deferred.

## Work counts remain informative

Counts are independently recomputed from the raw per-rank files. They repeat
exactly between the two observer runs at each count. Global Domain/patch
inventories, producer patch-layer visits, restriction-parent/edge work and
divergence replay counts match exactly at each ordinary step across 160 and
256 ranks. Selected remote payload totals conserve sent/received values.

At 256 ranks, every recorded ordinary step has 160 Domain owners and 96 ranks
without Domains. Verified placement splits these groups 80/80 and 48/48 across
the two sockets, respectively; the groups are not segregated onto one socket.
Ownership is measured, not inferred from rank IDs.

| Group at 256 ranks | Rank share | Block share | Block patch / divergence-replay share | Producer patch-visit share |
|---|---:|---:|---:|---:|
| Domain owners | 62.5% | 83.18% | 77.21% | 100% |
| Without Domains | 37.5% | 16.82% | 22.79% | 0% |

A Domain owner averages 3.327 blocks and 12.810 block patches per RK-stage
snapshot; a Domain-free rank averages 1.122 blocks and 6.303 patches. The extra
ranks do perform useful block work, but their divergence-replay workload per
rank is roughly half that of Domain owners. Equal block counts alone would not
mean equal cost; the independent patch and executed-work counters are essential.

The Domain owners' mean producer patch-layer visit count remains 1,991.083 per
ordinary step at both rank counts; Domain-free ranks execute zero such visits.
This confirms a structural parallelism ceiling in that producer. It does not
establish the producer's exclusive time or prove it accounts for the entire
47.035 ms whole-step loss. Retained communication, preparation, hardware
bandwidth/frequency and run variation remain possible contributors.

## A specific partitioning weakness

`build_parallel_block_catalog` in `src/parallel_block_mpi.f90` (target assignment
around lines 5243–5302) uses the original fixed average load for each rank,
consumes contiguous blocks in increasing rank order, reserves one block per
remaining rank, and relaxes an overload tolerance until every block fits.
It does not adjust the next target from remaining weight/ranks or cumulative
ideal cuts. Early ranks can consume excess weight, leaving later ranks below
target. In this campaign the Domain-free ranks occupy the later rank IDs.

A small offline counterexample reproduces the algorithm with **640 equal-weight
blocks on 256 ranks**. It assigns three blocks to the first 192 ranks and one
to the last 64. A feasible assignment has 128 ranks with three blocks and 128
with two, with the same maximum load. This proves a weakness in the assignment
rule even when the weight proxy is perfect. It is a synthetic counterexample,
not a replay of the J5 catalogue: Stage272 did not archive individual catalogue
weights, so the exact share caused by that rule has not been isolated.

The weight itself is also provisional: `subtree_weight_Domain` in
`src/domain.f90` counts descendant leaf patches times `PATCH_SIZE**2`. It does
not explicitly include fixed per-block overhead, communication, or retained
Domain-owner work. These are distinct limitations; changing all weights and
ownership rules together would make the next result harder to interpret.

## Next proposed implementation and test

Correct the rank-order front-loading with deterministic balanced contiguous
cuts using cumulative/remaining workload targets. Initially retain catalogue
order and the existing weight definition, so this is one isolated partitioning
change. Include unequal indivisible-block cases and the equal-weight
counterexample in offline qualification; retain valid ownership for every block.

Use actual catalogue weights and lightweight ownership/work evidence to verify
what changed, then compare **uninstrumented whole-step performance** against
Stage266 at 83/160/256, with supported legacy anchors at 83/160. Preserve exact
checkpoint/VTK, checked/oracle, restart/remap and zero-Domain-rank gates. Phase
profiler calibration need not precede this experiment: it is supported by exact
work counts and an independently reproducible algorithmic weakness. Reject a
regression before stacking another change. No speedup percentage is justified yet.

This correction alone cannot distribute the retained Domain producer beyond
160 owners. The subsequent substantial ownership change must replace a complete
producer-to-consumer chain, including its Domain-bound communication and capture,
while preserving Stage267 input epochs and numerical oracles. Efficient execution
at **at least 640 ranks** remains the goal, requiring finer useful granularity and
multi-node evidence; neither correctness at 256 nor a partition fix completes it.

## Runtime and retained evidence

Slurm records **00:18:52 (1,132 seconds)**, from 2026-10-08 12:29:35 to 12:48:27.
The runner independently records 1,131.675 seconds. Four serial builds and
fourteen runs fit the 16–22-minute estimate. Queue and collection/transfer are
excluded. The runtime ledger now retains the accounting record.

The [machine-readable review](parallel-blocks-stage272-results.json) records
calibration failures, individual timing samples, exact work totals, group shares,
source provenance and the counterexample. Archive SHA256:
`7c103ed35cce3e39b9387d1dda53bd741e8b8b5c3ee397628c7f6661d3d2ece8`.
The archived numerical reports and hashes were verified; omitted binary,
checkpoint and VTK payloads were not rerun locally.
