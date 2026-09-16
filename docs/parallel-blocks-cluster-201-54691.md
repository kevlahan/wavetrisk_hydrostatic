# Stage 201 cluster review — job 54691

Stage 201 provides stronger evidence for **extra representation/data-movement
work, amplified by MPI progress**, than for a general slowdown of every block
instruction. Memory costs remain important within wavelets. The aggregate IPC
alone would obscure both facts. The next target should be complete inverse
staging and scalar-storage execution paths, not another isolated message-volume
or arithmetic-kernel tweak. No source optimization or commit is made here.

## Validation and timing

Archive `stage201-review-54691.tar.gz`, SHA256
`e77dfc5fe1936486446250016b09014723804280b404cec8d54616712fab2664`.
All **6,931** inventoried files pass size/hash checks, with no missing or extra
files. All **20 runs and 42 gates pass**, including checked/optimized comparisons,
restart/remap, fixed workload and the existing legacy numerical screen. The
screen does not imply bitwise equivalence between implementations.

Independently reconstructed placement from topology, checked all 20 rank sets,
and reread the raw hardware records from all 12 observer runs and all 14 sampled
rank exports. Integer counter values and sample tallies reproduce the archived
results; timing arithmetic agrees within 1e-12 seconds. List order for tied
sample tallies and filesystem enumeration is not significant. Source/build
records match the frozen kit; bare binary hashes are consistent across repeats,
and the 23 common physics objects match. The job uses the prescribed 83 ranks,
42/41 across both sockets and all 32 cache groups on exclusive bb13.

| Version | Ordinary-step medians, seconds | Median |
|---|---:|---:|
| Blocks | 1.05, 1.07 | **1.0600** |
| Legacy | 0.547, 0.538 | **0.5425** |

The ratio is **1.954**, with **0.5175 s/step excess**, consistent with Stage 200's
best-placement result. These are two repeats, not precise confidence bounds.
The coarse rank-mean phase excesses are 0.2155 s wavelets, 0.2021 s remaining
dynamics, 0.0892 s adaptation and -0.00345 s elsewhere. Their 0.5033 s sum is a
different statistic from the bare printed-step median. The two dynamics phases
still contain about 83% of the coarse excess.

## What the hardware counters establish

The table gives ranges over the four core/memory measurement pairs. All counts
include user-space MPI progress/spinning as well as solver instructions.

| Phase | Block/legacy instructions | Block/legacy CPI |
|---|---:|---:|
| Dynamics wavelets | **3.68–3.88** | **1.013–1.029** |
| Remaining dynamics | **1.495–1.506** | **1.006–1.009** |

The excess aggregate cycles mostly follow additional instructions, rather than
a large aggregate CPI penalty. It would be incorrect to call all those extra
instructions redundant solver computation: MPI progress is a major contributor.
Likewise, the algebraic instruction/CPI cycle decomposition in the automatic
report is not a causal wall-time recovery estimate.

Memory evidence is also clear. Mean demand-fill counts per rank per ordinary
step, averaged over the two memory passes:

| Phase | L2 demand misses, blocks/legacy | DRAM/MMIO demand fills, blocks/legacy | Fill ratio |
|---|---:|---:|---:|
| Wavelets | 2.839M / 0.582M | 1.352M / 0.166M | **8.13** |
| Remaining dynamics | 6.163M / 5.011M | 1.217M / 0.312M | **3.90** |
| Adaptation | 1.180M / 0.268M | 0.859M / 0.070M | **12.30** |

These event counts are not total memory traffic or a bandwidth measurement.
The DRAM/MMIO event retains its documented MMIO qualification. Demand fills per
thousand instructions are 0.428 versus 0.196 in wavelets and 0.180 versus 0.069
in remaining dynamics: the increase is not explained only by executing more
instructions. Load-wait retirement cycles comprise about 27.8% versus 26.7% of
wavelet cycles, and 29.6% versus 33.5% outside wavelets. Those mixed-phase
fractions do not isolate solver load stalls from MPI polling, and waiting on a
load is not uniquely a DRAM bottleneck. No swapping occurred.

## Separating MPI progress from solver work

Classified a sample as an MPI/UCX stack if any frame resolves to libmpi,
libmpifort, libuct, libucp or libucs. Otherwise retain it outside that class.
This also assigns libc/TLS work called inside MPI to MPI, rather than treating
it as independent solver work. Other leaves are classified by their resolved
source symbols, preserving inline frames. This is a post-processing view of
the existing exports; no extra run or timing correction was needed.

The same seven selected ranks and eleven ordinary steps give:

| Phase / event | Blocks MPI / other samples | Legacy MPI / other samples |
|---|---:|---:|
| Wavelets, cycles | 1,453 / 1,688 | 254 / 571 |
| Wavelets, instructions | 7,661 / 4,484 | 1,356 / 1,831 |
| Other dynamics, cycles | 3,162 / 3,602 | 2,381 / 2,128 |
| Other dynamics, instructions | 15,647 / 10,885 | 11,374 / 6,347 |

Thus **63%** of block wavelet instruction samples are on MPI/UCX stacks. This
is mostly progress/wait execution, not evidence that transmitting fewer bytes
will save the same fraction. The leading block wavelet MPI caller is
`finish_exchange_values`, with 810 cycle and 4,308 instruction samples.
Writeback, vector-wavelet validation and ghost-exchange callers also appear.
MPI time can reflect local/remote producer delay as well as communication.

After excluding these MPI stacks, the sampled block wavelet instruction count
is still **2.45×** legacy, with an approximate sampled CPI ratio of **1.21**.
Remaining-dynamics instructions are **1.71×**, with sampled CPI ratio **0.99**.
The fixed sampling period is identical for both events and implementations.
These are diagnostic ratios for the selected cohort, not unbiased all-rank
estimates. Nevertheless they materially refine the aggregate result: wavelets
show both excess instructions and a local execution penalty; other dynamics
show principally excess instructions in this cohort.

## Concrete source families

### 1. Complete native inverse staging, aliasing and publication

Native inverse own leaves, including the compiled local route-copy helper,
account for **25.5% of block wavelet cycle samples**. Other block interfaces add
14.6%; these are disjoint sampled leaf groups, not overlapping parent timers.
The leading leaves include `transfer_nodes`, `alias_value`, `apply_aliases`,
`copy_values`, `prepare_native_inverse`, `set_alias`, and scalar-record resolution.

The source explains what this work does:

- `parallel_block_inverse.f90:708`: `transfer_nodes` packs/scatters fields and
  levels between inverse staging and routed block storage.
- `parallel_block_inverse.f90:988`: `apply_aliases` snapshots remote sources,
  performs ordered local copies through alias decoding, then installs receives.
- `parallel_block_inverse.f90:366`: preparation seeds compatibility/scaffold
  state, gathers scalar/vector families, applies aliases and scatters them;
  consumer-visible alias state is retained separately.
- `parallel_block.f90:958`: compiled routes still execute field copies between
  storage representations. Compiling routing metadata does not eliminate them.

These are representation and dependency-maintenance operations, rather than
another full inverse arithmetic evaluation. Some preserve necessary historical
alias epochs and scaffold semantics. They cannot simply be deleted. Their own
work and subsequent waiting identify a **whole-path redesign candidate**:
reduce staging/translation between authoritative owner storage and transform
consumers, retain prepared addressing, and publish only at required consumer
epochs. Removing only a few arithmetic operations or message bytes would leave
much of this structure intact.

### 2. Scalar storage adapters throughout the dynamics pipeline

Scalar storage/addressing leaves account for **16.5% of block remaining-dynamics
cycle samples**, plus 4.8% in other block interfaces. Examples are `read_value`,
inlined `shared_address`, `write_range_values`, `write_value` and
`scalar_read_range`. The wider tendency-capture routine is a prominent inclusive
caller, but its samples must not be added to these leaf totals.

`parallel_block_scalar_storage.f90:179` reconstructs mixed 50-slot records via
individual `read_value` calls. Geometry reads repeatedly decode node/field/zone
addresses at line 122. The analogous write fallback at line 334 visits slots
individually, including geometry writes that are then rejected because another
canonical field owns that geometry. The direct-flux caller at
`parallel_block_mpi.f90:24098` reads the full record before its mass shortcut and
nonphysical-level rejection, consistent with the Stage 199 access census.

This is concrete excess execution caused by presenting compact field/shared-
geometry storage through an old full-record interface. The target should span
read, write, capture, restriction and flux consumers: prepared field/geometry
views and bulk operations, without repeatedly materializing legacy-shaped
records. Preserve the useful memory sharing; merely expanding all storage back
to full records would trade away memory capacity without testing the intended
solution. Nor should the next stage be limited to the approximately 55 ms
single direct-flux batch already identified in Stage 199.

### Other costs and offsets

Generic Domain boundary work remains present in blocks, but legacy spends more
sampled own-leaf cycles there: 302 versus 47 in wavelets and 822 versus 513 in
other dynamics. Those savings offset some new block costs and must be retained
in a signed comparison. This is not evidence that every retained boundary call
is redundant. Memory allocation/copy/zero also appears, especially outside
wavelets, but is not a single independently removable amount.

## Measurement qualifications

All recorded PMU groups have **100% running/enabled time**; no multiplex scaling
is needed. No lost samples or shallow stacks are reported. Numerical PASS and
valid counters do not waive calibration review:

- Block PMU printed medians equal its 1.06 s bare median. One wavelet envelope
  changes by +7.5 ms (2.63%); the other block PMU phase changes are small.
- Legacy PMU medians are 1.38–2.30% above bare. Wavelet changes reach +3.68 ms
  (5.32%); other-dynamics changes are at most +2.74 ms (0.73%).
- Sampling changes printed medians by +0.94% for blocks and +4.33% for legacy.
  Other-dynamics envelopes rise by 12.4 and 14.1 ms respectively.

These effects do not explain the multiple-fold wavelet instruction excess or
50% other-dynamics excess. They do prevent claiming precise per-function wall
seconds or applying a universal overhead correction.

The unresolved-leaf flags are mostly missing **UCX transport symbols**, not
unidentified solver code: 6,375 of 6,790 unresolved block leaves and 3,601 of
3,808 unresolved legacy leaves are in known libuct. Their library and caller
stacks remain identifiable, allowing the broad MPI/solver separation above.
Exact transport-routine attribution remains unavailable. Small additional
unknown counts occur in libm and a few other libraries. Seven stratified ranks
are not a statistical sample of the full node.

## Decision

The campaign has answered the main discriminator sufficiently to select source
families: there is substantial extra software work in representation adapters
and inverse state movement, with increased memory demands and progress waiting.
It does not support a general claim that the numerical arithmetic is twice as
slow, or that network message volume alone explains the gap.

Do not repeat the full scaling or profiling campaign. The next useful step is
a bounded, correctness-gated **whole-path representation-removal prototype**,
prioritizing inverse staging/publication and the scalar storage adapter pipeline.
Specify which copies/record reconstructions disappear, which consumer epochs
remain, and which source counters should fall before implementing it. Compare
against this unchanged candidate and legacy at the same placement, with two bare
pairs and only the matching counter checks needed to test the mechanism.

A demonstrated **0.10 s/step** reduction remains the acceptance target. The
current samples identify sizable gross costs; they do not establish how much
can safely be removed, and MPI waiting cannot be added to producer savings as
an independent benefit. No narrow getter fix, single exchange deletion, or
production speedup is approved by this review alone.

Machine-readable evidence: `parallel-blocks-cluster-201-54691.json`.
