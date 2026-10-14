# Stage 278 structural preparation repair v7

Prepared and locally qualified against the exact attribution-only source measured
in job 55503. **Stage 266 remains accepted; cluster qualification is pending.**
The v6 canonical-fallback shortcut remains preserved separately and is excluded.

## Why this repair

Job 55503 measured RHS route construction at 151.61/92.75/84.20 ms at 83/160/256
ranks, plus input-plan construction at 34.36/26.98/31.67 ms. The remaining bare
cuts to reach 266 were 43.71/52.79/30.01 ms. Preparation is a material budget;
these profiled means are not directly subtractable bare-time predictions.

An exact cross-epoch route-cache experiment compared producer homes/order/slots,
reader homes/destinations/signs, canonical slots and buffer layout, with global
consensus. **No route was reused across 12 adaptive preparations:** 45 initial builds and
495 subsequent checks. Its science
replay passed, but the cache was rejected: no demonstrated reuse justified its
memory and consensus costs. This does not rule out finer local substructure
reuse; unchanged counts do not establish reusable route identity.

The repair instead removes a general final-writer/halo-plan compiler from
single-use physical geometry metadata transfers. Their unique producer rank and
address are already determined by the current geometry layout. One directory
per preparation supplies those addresses to direct indexed pulls. Repeated
requests are deduplicated by exact rank/slot, preserving the old compiler's
remote metadata value count. A second structural fix replaces the existing
rank/slot deduplication hash in both this pull and the generic RHS route compiler.
The old linear hash clusters contiguous address ranges; the mixed hash retains
exact key comparisons and first-appearance ordering. A targeted dense-address
benchmark needed 1,154,000,000 probes with the old hash versus 471,980 with the
new one for 20 passes over 16,000 keys. Other tested sizes include layouts where
the old hash is already fast; synthetic results are not solver speedups.
The campaign's cost files contain 32 rank/epoch rows with 993 program nodes
at 160 ranks, and 89 with 140 nodes at 256 ranks. Full vector-alias metadata
queries therefore create tables of 5,959 and 841 entries. The old rank coefficient
65599 reduces to strides of 50 and 1 in those tables. This ties the vulnerability
to observed campaign layouts; the archive lacks actual lookup probe counts,
so it does not quantify their contribution to the scaling loss.
Numerical-stage communication is unchanged.

## Dependencies, invalidation and bounded lifetime

| Dependency | Producer and readers | Validity rule |
|---|---|---|
| Physical owner, offset and node extent | Current geometry owners; alias, provider, retained-geometry and surface readers on block owners | Construct once during every RHS preparation; expire on return |
| Physical key | Dense `[Domain+1,node,component,1]` source namespace | Check every source key, current local layout and requested extent |
| Alias/mask metadata | Existing accepted boundary and geometry compilers | Recompute from current grid; no retained plan selected by counts alone |
| Numerical RHS/RK | Existing block owners and direct fused kernels | Placement, arithmetic, reader demands and state/halo stamps unchanged |
| Metadata epoch | Preparation caller | Distributed topology, phase and column-width agreement checked before transfer |

Legacy writeback invalidation remains conservative. No real adaptive invalidation
is bypassed. The directory, requests and replies have preparation/call lifetimes;
there is no added persistent cache, borrowed state pointer or MPI communicator.
No Domain numerical pathway is restored and no new ownership retirement is claimed.

For the three-level J4/J5 fixtures, 16 specialized pulls replace disposable halo
plans: six boundary aliases, five resident/RK provider lookups, two input
boundaries, two retention records and one surface field. Explicit construction
collectives fall **192 → 65 per preparation**, including the new layout consensus:
127 fewer. Sixteen communicator duplication/free pairs disappear. This counts
removed construction work; it does not identify network time as the bottleneck.
Catalogue/provider-key transfers outside this physical namespace retain the
generic compiler, including the transform-geometry path; its rank/slot
deduplication receives the same hash correction.

The directory's array payload is 3D default integers: **1,920 bytes at D=160**,
with a 4,480-byte array peak while constructing it. Indexed-pull temporary integer
arrays are conservatively bounded by `4*(5P+9Q+R+1)` bytes, including the caller's
home array; double arrays use `8*C*(R+U+Q)` bytes including returned values.
Q is local queries, U≤Q unique queries, R incoming unique requests, C metadata
width. Source arrays, small scalars/descriptors and MPI internal storage are
additional. All buffers expire after the call. This is a storage/lifetime bound,
not an observed whole-process RSS improvement.

## Final local evidence

Serial builds; at most four MPI ranks. All twelve full replays pass: matched coarse
control/candidate, bare ABBA+BAAB, candidate independent oracle and candidate checked.
They retain twelve exact accepted-family states, exact checkpoint fields/topology,
two reloads, two remaps and 62 VTK members. No main-relative or above-Domain full
solver qualification is inferred from these local tests.

| Ordinary-step run median, seconds | First | Second | Third | Fourth | Mean |
|---|---:|---:|---:|---:|---:|
| 55503 attribution-off control | 4.208613 | 4.324102 | 3.701969 | 4.176575 | 4.102815 |
| Final v7 repair | 4.325962 | 4.118669 | 3.971342 | 3.480904 | 3.974219 |

**The mean of ordinary-step run medians changes by -3.134%**. Control repeat spread is
15.164%; candidate spread is
21.263%. These are four run repetitions/family,
not independent timestep samples or a confidence interval. The median point
estimate is smaller than the observed run-to-run spread, so it does **not**
establish a reliable typical-step speedup. The additional BAAB
set was predeclared because the first final ABBA set's median gain was small
and variable while ordinary total time improved; all eight timings are included.
The mean ordinary-step time (equivalently total time over the same nine ordinary
steps, a separate statistic) is 4.742523 →
3.975427 s, a -16.175% change.
Every candidate ordinary-step run mean is below every control run mean:
control range 4.445587–5.503657 s,
candidate range 3.510378–4.262204 s.
This supports an ordinary-throughput improvement in this local sample, but the
aggregate is sensitive to the long fourth control and the fast fourth candidate.
Neither the aggregate nor the noisy median is a cluster speedup forecast.
The final replays
include metadata deduplication with a mixed rank/slot hash; the earlier prototype's results are archived
under `initial-proposal/` and are not the final source's timing claim. A first
deduplicated variant recovered only 1.63% ordinary-step time and 9.21% route
construction. Dense source addresses produced long linear-probe clusters.
That variant is retained under `poor-index-experiment/`; the final version
mixes the slot before table reduction and retains exact key comparisons.
The targeted synthetic dense-address benchmark and all intermediate results
are preserved; synthetic times are not production speedup evidence.

Remap/non-remap grouped run means are a separate statistic:
non-remap -19.104%, remap -4.044%.
Ordinary steps exclude the established startup/reload/save indices and retain
both ordinary remaps.

| Coarse self region on the same slowest whole-step rank, ms | Control | Repair | Change |
|---|---:|---:|---:|
| Route construction | 445.431 | 221.195 | -50.34% |
| Input-plan construction | 187.369 | 73.901 | -60.56% |
| Retention-plan construction | 58.775 | 5.702 | -90.30% |

The profiled whole-step mean is 4.895678 →
4.463474 s. The new indexed value pull falls within its
preparation parent; the old one-use halo exchange was nested in region 78.
Do not reinterpret region 78 as pure numerical transport or subtract these
regional changes directly from the bare medians.

The arithmetic region is 754.090 →
757.228 ms; the measured exchange region is
833.633 → 832.036 ms.
The observed reduction concentrates in preparation, consistent with the source
change rather than an assumed improvement in numerical transport.

Optimized/checked indexed-transfer tests compare exact values **and remote value
counts** with the independent generic compiler across 40 changing-layout cycles
at one/four ranks. They cover empty providers/readers, duplicate requests and
bad-address/distributed epoch/phase/width rejection. Real geometry-directory tests
cover 12 ownership/order/extent changes, repeated topology labels, one/three
components, zero-geometry consumers and stale-directory/layout/key rejection.
All pass. The existing optimized/checked exchange suite also passes its
negative epoch/identity/ordering checks and large rendezvous stress test.
All integer work/ownership/numerical-traffic/narrow-storage fields match
the prior qualified attribution-only oracle: four 48-row production files and
four 12-row cost files. This narrow ledger does not establish whole-process RSS.

## Dual-goal decision and cluster gate

The repair preserves connected block-owned RHS/RK and every existing paired
Domain retirement while reducing recurring preparation costs. It does not
complete whole-timestep independence or restore Domain-owned numerical work.

The 55503 gaps correspond to roughly 29%/57%/36% of the route-construction budget
at 83/160/256. The local route reduction is 50.3%, with
additional input/retention savings. This supports a material 266-parity campaign
candidate; the 160-rank case needs those additional savings. It is a screening
comparison across machines/statistics, **not proof that any cluster gap is closed**.
The local median point estimate is a 3.13% reduction but remains noisy;
ordinary-step total time improves in every candidate versus every control in this
sample. Performance acceptance stays open until the matched cluster measurement.

There is a concrete reduction in work that grows expensive with rank count:
127 fewer construction collectives and 16 fewer disposable contexts per preparation,
without increasing remote metadata value counts. Useful MPI strong scaling is
still unverified. Numerical RHS waits, inverse/writeback dependencies and producer
imbalance remain. Stage 266 also regressed at 256; 78–83% of the prior main gap
already existed in 266. Legacy parity, total RSS and efficient ≥640-rank execution
remain open; no unconditional scaling forecast or numeric efficiency threshold
is invented.

The next cluster test should use matched bb13 placement at 83/160/256, repeated
bare v5/repair comparisons and 266/main controls, preserving checked/oracle and
zero-Domain science gates. Report absolute ordinary-step times, remap groups,
160→256 efficiency and rank RSS separately. Coarse profiles should verify the
specific preparation saving. No broad fine-profile repeat is needed.

## Reviewable source and evidence

The six-file patch applies to the exact 55503 attribution candidate and reproduces
all qualified source hashes. The working implementation contains the final repair.
Earlier proposals, the rejected cache and the first prototype remain preserved.
Stage 266 and all sealed kits are unchanged. No remote upload or cluster job ran.

[Exact patch](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/repair-v7/delta-from-55503.patch), [source identity](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/repair-v7/IDENTITY.json),
[qualification](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/repair-v7/qualification.json), [additional timings](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/repair-v7/timing-repeat.json), [metrics](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/repair-v7/metrics.json),
[indexed tests](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/repair-v7/indexed-units/results.json), [geometry tests](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/repair-v7/geometry-units/results.json),
[exchange tests](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/repair-v7/exchange-units/results.json),
[numerical counters](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/repair-v7/traffic-equivalence.json), [campaign hash layouts](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/repair-v7/campaign-hash-layouts.json).
