# Stage 181: targeted live-allocation census

Stage 180 was committed as `6c3b3f8c`. This next stage prepares a diagnostic-only
archive; it does not change the working repository's numerical Fortran sources,
authoritative main, or the numerical acceptance limits.

## Why this measurement, rather than more broad timing profiles

The existing profiles already identify a concrete block-specific target: the
50-double scalar workspace repeats 33 horizontal-geometry/metadata slots across
scalar variables and vertical levels. We do not need another broad profiling
campaign to choose that target. We do need a live-allocation accounting to
distinguish those records from retained block state, Domain compatibility data,
and concurrently allocated transport buffers.

`prepare_memory_census.py` generates size inquiries inside an isolated copy of
four owning modules: `domain`, `comm_mpi`, `parallel_block`, and
`parallel_block_mpi`. It inventories 548 allocation paths from their declarations.
It counts allocated capacity (not only used lengths), including nested owned
allocatables and their compiler-reported inline storage. Pointer targets are
excluded so aliases of existing allocations are not counted again. A manifest
records every included allocation path and excluded pointer path. Unknown owned
types/declarations fail preparation rather than silently disappearing.

The generated code inspects allocation descriptors and sizes; it does not scan
field values, mutate numerical arrays, or introduce MPI collectives. It writes
per-rank records at scalar-storage preparation and timestep end. Its I/O and
descriptor traversal make this a **diagnostic executable, not a timing build**.
An explicit `WAVETRISK_ALLOCATION_CENSUS=1` enables the sampler; the ordinary local
protocol continues to clear inherited diagnostic switches.

## Scope and interpretation

- Within each rank, all categories in one sample refer to the same numerical
  execution point, with no solver work between their inquiries.
- Report each rank's sampled peak and the allocations present at that peak.
  Do not sum independently peaked categories or ranks and call that a live peak.
- If rank sample sequences align, also report the sum at matching logical
  samples. This is **not a synchronized wall-clock observation**; no barrier is
  added. If sequences differ, omit that sum.
- This is a four-module census, not total process memory. Other modules
  (including native inverse/velocity/mass workspaces and external RK storage),
  automatic/local temporary arrays, allocator-retained memory/metadata, MPI
  internals, libraries and OS memory remain outside the count. Transient peaks
  between sampling points can be missed. RSS/compression measurements remain
  complementary, not interchangeable with allocated capacity.

## Reproduction

Prepare from the same verified block archive used for Stage 180:

```sh
python3 test/parallel_block_profile/prepare_memory_census.py \
  --repo "$PWD" --baseline /private/tmp/wavetrisk-stage180-block-opt-01 \
  --out /private/tmp/wavetrisk-stage181-census-03
make -C /private/tmp/wavetrisk-stage181-census-03 \
  -j1 DEBUG=false PARAM=param_J5 TEST_CASE=climate
python3 test/parallel_block_profile/run_memory_census.py \
  --build /private/tmp/wavetrisk-stage181-census-03 \
  --reference /private/tmp/wavetrisk-stage180-j6-campaign-01/block \
  --fixture /private/tmp/wavetrisk-stage180-j6-fixture-01/seed \
  --out /private/tmp/wavetrisk-stage181-census-run-01
```

Use new output directories. The diagnostic runner verifies the baseline binary
identity, fixture/grid hashes, rank count and interval. It runs four MPI ranks
sequentially with all three oracles disabled and then requires **zero numerical
difference in the compared checkpoint fields** against the existing block
baseline at both checkpoints 5 and 6. This tests instrumentation transparency;
it does not waive the separate unresolved block-versus-legacy discrepancy.

## Next performance change to prepare

Replace the scalar execution workspace's repeated geometry with shared
horizontal-node storage while retaining the 17 per-field slots. Carry this
through interior, boundary and ghost storage, inverse/restriction consumers,
and oracle access. Preserve operation order, ownership and generation lifetimes.
In particular, preserve existing zero/mask semantics for soil/surface slots:
sharing geometry must not accidentally activate fields previously zeroed there.

Acceptance requires candidate-versus-baseline numerical equivalence, independent
RK/adaptation/remap oracles, restart coverage and an actual decrease in counted
capacity. The legacy comparison remains separately visible and unresolved.
Reliable speedup claims require a pressure-free timing case and later the
matched 83-rank cluster comparison. No runtime speedup follows directly from a
percentage reduction in allocated bytes.

## Local results

The isolated optimized build succeeded. All 42 helper tests pass, including a
compiled Fortran allocation-capacity regression with bounds checking, nested
allocatables, empty arrays and non-unit lower bounds. Parser tests reject
incomplete/duplicate/negative samples and explicitly distinguish matching-phase
sums from sums of independent rank peaks.

The four-rank J6 run completed both checkpoint/reload cycles. All compared
checkpoint fields, headers, thresholds and topology match the existing block
baseline exactly at checkpoints 5 and 6. No numerical acceptance tolerance was
changed. This is instrumentation transparency, not new block-versus-legacy
acceptance; the known small legacy discrepancy is unchanged.

All four ranks produced 55 complete samples with aligned phase sequences. Their
largest matching sample is number 42 (`scalar-storage-ready`):

| Allocation category at that logical sample | GiB, sum over four ranks |
| --- | ---: |
| All counted allocations in four modules | 10.487 |
| Scalar workspace including coverage/inline metadata | 6.517 |
| Its 50-double scalar records alone | 6.504 |
| Retained local block objects and their arrays | 0.985 |
| Staged coarse-vector boundary plan | 0.346 |

The scalar records comprise 1,324,070,400 interior bytes, 3,321,590,400 boundary
bytes, and 2,337,984,000 ghost bytes. Thus boundary/ghost records account for
81.0% of this allocation. Current source/received migration-block storage is
zero at this sample, and the scalar divergence transport-plan root totals only
1,687,588 bytes. The dominant persistent storage is the **expanded computational
workspace**, not that particular transport buffer.

For this two-scalar, 41-field-level layout, retaining 17 field slots and sharing
the 33 geometry slots once per stored horizontal node estimates a **4.240 GiB
reduction** in these records, leaving approximately 2.264 GiB. This is a layout
estimate from allocated record counts, not an implemented result, a whole-process
memory prediction, or a runtime speedup claim. Mask/soil semantics and additional
mapping storage still require audit.

No additional broad profiling is required before choosing this target. A
controlled layout experiment should use exact candidate-versus-block-baseline
checks, with the unresolved independent legacy comparison retained separately.
It must not be described as production acceptance while that gate remains open.

Raw evidence: `/private/tmp/wavetrisk-stage181-census-run-01/`, with preparation
coverage/source hashes in `/private/tmp/wavetrisk-stage181-census-03/`.
The sibling `parallel-blocks-allocation-census-181-results.json` retains a compact
copy of the numerical gate, identities, memory observations and allocation
summary. Stage 181 changes are left uncommitted.
