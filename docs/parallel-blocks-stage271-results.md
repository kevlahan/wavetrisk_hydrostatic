# Stage271 result: no consistent gain; scaling deteriorates beyond the Domain count

Job 55468 completed all twenty runs. Independent review verifies all **5,314**
archive files, reproduces the scientific/placement/storage/timing checks, and
matches the delivered identity, source and kit manifests byte for byte. Four
serial builds share 23 physics objects. Optimized and checked candidates pass;
both 16/256-rank oracles retain 48 exact inverse-return proofs and 84
installations. Exact block checkpoint, same-count 62-member VTK, trajectory,
restart/remap and zero-Domain-rank checks pass.

The experiment is numerically valid but **not promoted for performance**.
Restore the exact accepted Stage266 solver source (`ff740744`) before further
implementation. The fused candidate remains recoverable in commit `ab9c8463`
and the immutable Stage271 kit. No user edits existed in the restored file;
both its prior candidate hash and restored baseline hash match the reviewed
source manifest. This restores already-tested code and introduces no new solver
implementation requiring another local replay.

## Uninstrumented performance

Times are milliseconds per ordinary step. Each block entry is the mean of two
independent run medians; main has one anchor at each supported count. Steps 1,
8 and 9 are excluded. Oracle and checked runs are excluded from these timings.

| MPI ranks | Legacy main | Stage266 | Fused candidate | Candidate change |
|---|---:|---:|---:|---:|
| 83 | 539.122 | 755.298 | 760.462 | +5.164 ms (+0.684%) |
| 160 | 356.142 | 549.270 | 546.914 | −2.357 ms (−0.429%) |
| 256 | — | 584.620 | 583.426 | −1.194 ms (−0.204%) |

At 83 ranks both candidate pairs are slower: +8.918 and +1.411 ms. At 160 both
are faster, by 3.350 and 1.364 ms, but that modest change does not transfer
consistently to 256: its pairs are +1.971 and −4.358 ms. Candidate run spans
are 4.961, 4.500 and 6.772 ms respectively. These results do not establish a
useful consistent gain from fusing the three proof reductions. They also do
not justify claiming a statistically established general regression from one
allocation and two runs per configuration. Do not repeat this micro-optimization
solely to resolve noise or stack further edits on it.

The accepted baseline remains substantially slower than main: +216.176 ms
(40.098%) at 83 ranks and +193.128 ms (54.228%) at 160. The main anchors are
single-run screens. Legacy main was not run beyond its supported Domain count.

## Scaling conclusion

For J5 there are 160 geometric Domains. Increasing from 160 to 256 ranks adds
**60% more ranks**, but Stage266 step time increases by **35.350 ms (6.436%)**.
The candidate likewise slows by 6.676%. Both 256-rank baseline samples are
slower than both 160-rank samples. The small collective change does not alter
the scaling conclusion.

| Block version | Speedup, 160 → 256 | Relative parallel efficiency |
|---|---:|---:|
| Accepted Stage266 | 0.9395× | 58.72% |
| Fused candidate | 0.9374× | 58.59% |

Efficiency is `T(160)*160 / (T(256)*256)`, relative to the 160-rank configuration,
not serial efficiency. From 83 to 160 ranks, Stage266 achieves 1.3751× speedup
and 71.33% relative efficiency. These are observations for this fixed trajectory
on bb13, not a claim about all resolutions or nodes. Rank groups ran in ascending
count order, so cross-count comparisons retain potential time/order/load effects;
within each count, block comparisons used A-B-B-A. Distinct physical cores and
the prescribed two-socket placement were verified for every rank.

The strategic goal of efficient execution at four or more ranks per Domain is
unmet. J5 needs at least 640 ranks for that goal, plus adequate independent block
granularity and eventual multi-node validation. Passing the 256-rank oracle
establishes correctness, not efficient scaling.

## Reasons being investigated, with evidence boundaries

The source supplies a concrete architectural explanation for a ceiling in parts
of the computation:

- `prepare_block_multistage_tendency` still calls the production
  `block_tendency_compatibility_ml(domain_stage)`, then finalizes scalar capture
  and prepares the native velocity remainder.
- `basic_operators` traverses `size(grid)` and each Domain's patch lists for
  pressure integration, fluxes, vorticity, Bernoulli and Exner products. That
  useful producer work remains assigned to geometric-Domain owners.
- `prepare_native_mass` sizes native workspaces from `size(grid)`, as does the
  preparation of native velocity programs. Native arrays and destinations do
  not by themselves distribute all computation over final block owners.
- Stage270 measured about 194 ms of shared producer self time, 43.906 ms of
  capture preparation/publication inclusive time, 87.574 ms of scalar restriction
  inclusive time, and 51.746 ms of inverse wait/traffic self time at 83 ranks.
  These identify meaningful places to examine; they are not additive independent
  savings or measurements of those regions at 256 ranks.

Extra block owners can therefore contribute to block-native phases while still
depending on a producer restricted to at most 160 Domain-owning ranks. Transfers
between ownership layouts, repeated preparation and collective waiting can
compound this ceiling. Fusing two redundant collective invocations has not
materially changed the whole-solver result. This does **not** rule out other
communication costs or prove what fraction of the 35.350 ms loss belongs to
each mechanism.

The archived per-rank CPU/RSS records were examined. They span the complete
process, including initialization, I/O/reloads and MPI spinning; they are not
ordinary-step useful-work measurements. They cannot safely distinguish idle
Domain-free ranks from useful block work or identify bandwidth saturation.
Placement verification rules out obvious oversubscription, but does not rule
out cache, NUMA, all-core frequency or memory-bandwidth effects. Per-rank Domain
and block ownership/work counts were not captured by this uninstrumented kit.

## Next substantial direction

Prioritize the **shared Domain producer → capture → block restriction** ownership
chain. Before choosing its first complete replacement, obtain targeted
Stage266-only phase and per-rank work evidence at **160 and 256 ranks**, with
observer overhead checked separately at each count. The 83-rank calibration is
not transferable automatically. No further main-profiler calibration campaign
is needed to answer this scaling question.

The diagnostic should distinguish Domain owners from ranks without Domains using
actual ownership records, count useful Domain/block/patch work and exchanged
payloads, and report distributions of computation and waiting for those groups.
Measure whole transactions and bounded batches; avoid clocks in node/layer loops
and added timed-step synchronization. If producer time stops decreasing while
Domain-free ranks wait, move that complete producer and its consumers to block
owners. If replicated preparation or communication grows instead, address that
specific design cost first. Do not interpret CPU spinning as useful computation.

Preserve Stage267's numerical-reader, boundary-epoch and oracle constraints.
Evaluate each resulting intervention uninstrumented at 83/160/256 before
acceptance. The aim is to remove the architectural bottleneck with a complete
replacement, not accumulate additional small omissions and optimize afterward.
The current evidence determines priorities but does not yet measure a recoverable
speedup or prove the dominant cause of the beyond-Domain slowdown.

## Runtime and evidence

Automatically collected Slurm accounting records **00:21:47 (1,307 seconds)**,
from 2026-10-08 10:15:09 to 10:36:56. The runner independently records
1,305.975 seconds from 14:15:10.383472 to 14:36:56.358387 UTC. Slurm elapsed
includes batch setup; both exclude queue wait and collection/transfer. This
matches the 20–25-minute planning estimate for four builds and twenty runs.
The [runtime ledger](parallel-blocks-campaign-runtimes.json) retains the record.

The [machine-readable review](parallel-blocks-stage271-results.json) preserves
both pairs at every count, scaling arithmetic, accounting, resource-record
limitations and the rejection decision. Archive SHA256:
`6a6215971760327a3797dadb26fa40a82a30630e0ac75faabe5e93550b6f6654`.
The review omits binary/checkpoint/VTK payloads; their cluster-generated hashes
and numerical comparison reports were verified, not the omitted payloads rerun.
