# Stage 199 cluster review — job 54687

Stage 199 passes its validation campaign and narrows several suspected causes,
but does not establish a large, safely recoverable fraction of the production
gap. No production optimization or commit is made as part of this review.

## Evidence and validity

Input: `stage199-review-54687.tar.gz`, SHA256
`85a388b04feb82a902ba3ae964910d69afd296aa889de04e5af19965991f5e12`.
All 11,261 inventoried files passed size and SHA256 checks. All 27 runs and 65
gates passed. Same-implementation checkpoint comparisons are exact; cross-legacy
comparison passes the existing fixed numerical screen, not bitwise equivalence.
The run includes checkpoint 3 input, checkpoint 4 save/reload, and post-restart
remapping. There are 12 steps, including one restart step excluded from ordinary
step analysis.

Seven builds use J5, RK4 and contraction off, with common physics objects. The
optimized candidate build log contains `-O2 -march=native -ffp-contract=off`.
Build and execution host is bb13, AMD EPYC 9754. All 27 runs use the same
83-rank CPU affinity, distributed across the two sockets. No node swap-in or
swap-out counter increased during any run. Peak solver RSS in bare candidate
runs is about 390,000 KiB/rank, versus 173,000 KiB for legacy; this is not
evidence of memory pressure on this allocation. Cache and bandwidth costs are
not measured by these snapshots.

The four audited production source files (`multi_level.f90`,
`parallel_block_mpi.f90`, `parallel_block_scalar_storage.f90`, and
`parallel_block_build.f90`) match the archived candidate source hashes.
Production remains the Stage 197 candidate; Stage 199 changes observers only.

## Production timing and calibrated accounting

| Version | Ordinary-step medians, s | Median, s |
|---|---:|---:|
| Stage 196 baseline | 1.26, 1.26 | 1.26 |
| Stage 197 numerical candidate | 1.25, 1.25 | 1.25 |
| Legacy Domain | 0.557, 0.563 | 0.560 |

Blocks/legacy = **2.232**, an excess of **0.690 s/step** (123.2% slower).
This reproduces Stage 198's 0.689 s gap. It is not a new optimization result.

The exclusive operation ledger sums to 0.6335525 s excess. Selecting the actual
slowest full-step rank separately for each step gives 0.6337727 s. Neither is
the same statistic as the bare printed-step median. Positive location
differences sum to 0.866582 s and negative offsets to -0.233029 s; ignoring the
offsets exaggerates the gap.

| Matched phase | Legacy, ms | Blocks, ms | Excess, ms |
|---|---:|---:|---:|
| Dynamics wavelets | 81.33 | 343.73 | 262.40 |
| Remaining dynamics | 412.37 | 671.26 | 258.88 |
| Adaptation/reconstruction | 66.03 | 180.84 | 114.81 |
| All other phases, net | 22.48 | 19.94 | -2.54 |

These are rank-mean ordinary-step times in the timer-only runs. The two dynamics
phases account for 82.3% of the **instrumented net excess**, not a demonstrated
82.3% removable fraction of the bare gap.

Detailed timers change printed candidate medians by +1.2% but legacy by +7.68%
relative to their bare runs. Legacy parent envelopes also change appreciably.
Work/event recording adds about 6 ms to outside-wavelet dynamics on each side.
There is no defensible global overhead correction. Storage census has the same
1.285 s printed median as its own disabled-census control; this does not prove
zero overhead inside each individual batch.

The disjoint ledger can also be split by the explicitly timed MPI operation IDs:

| Phase | Excess outside designated MPI scopes, ms | Excess in designated MPI scopes, ms |
|---|---:|---:|
| Dynamics wavelets | 138.16 | 124.24 |
| Remaining dynamics | 193.32 | 65.56 |
| Adaptation/reconstruction | 74.95 | 39.86 |
| Complete step | 402.29 | 231.27 |

“Outside MPI scopes” includes arithmetic, copies, lookup, allocation, validation,
observer overhead and unclassified work; it is not hardware CPU execution time
and does not exclude MPI calls hidden inside unsplit routines. “MPI” includes
waiting for late producers and MPI progress, not just data-transfer cost. These
columns cannot be interpreted as independent recoverable savings.

## What Stage 199 resolves

### 1. The retained basic numerical operators are not the large deficit

Inclusive basic-operator batches cost 116.97 ms in blocks versus 138.56 ms in
legacy. Their non-designated-MPI portions are 45.80 versus 57.88 ms. These
scopes contain different supporting work, so they are not pure matched kernel
benchmarks, but they do not support attributing the 259 ms remaining-dynamics
deficit to a generally slower `basic_operators` call.

The more tightly matched pressure/step1 patch scope costs 21.41 versus 19.94 ms
in timer-only runs, a 1.47 ms difference. Inclusive Bernoulli and Exner
restriction together cost 0.644 ms in blocks and 0.721 ms in legacy. Native
mass restriction adds only 0.326 ms. These are not substantial targets.

Block-specific scalar divergence capture costs 17.18 ms and scalar physics
capture about 1.59 ms. The wider block scalar pipeline, storage and interfaces
remain relevant; optimizing the shared primitive arithmetic alone would miss it.

### 2. The repeated mass-edge exchange is real but small

The raw event trees distinguish the exchange inside `basic_operators` from the
one immediately after its return. In each work pass, both have exactly
349.09 calls per rank per ordinary step, with no mismatched `(step, layer,
level)` counts. The second post-plus-finish batch costs 16.85 and 16.63 ms,
including 12.76 and 12.54 ms in the wait scope.

Source inspection shows the intervening optional rotational-Laplacian path
operates on vorticity and its Laplacian, not the native mass-flux array. This
supports treating the second exchange as a redundancy candidate. It is not
evidence that removing it saves its full measured duration: synchronization
can move to the next dependency. Even its whole measured batch is only 2.4%
of the bare 0.690 s gap. It should not become the headline next optimization.

### 3. Scalar storage expansion is concentrated in direct flux

The two integer-census passes agree exactly on the work counts. Per rank per
ordinary step, the direct-flux batch requests **165,855 complete 50-slot
records**, causing **8,292,760 scalar reads** and **5,473,221 shared-geometry
address resolutions**. This one batch accounts for 92.5% of full-record range
requests, 82.1% of scalar reads and 77.7% of shared-address resolutions in the
tendency context. Its outer batch costs **55.28 ms** in the work passes.

The source explains the excess access work: `scalar_read_range` reconstructs
the full mixed geometry/field record through scalar accesses, before checking
the production mass shortcut or rejecting nonphysical field levels. Production
mass uses one returned slot; physical temperature flux uses three geometry and
three physics slots. Only 60,679 of the 165,855 samples proceed through physical
temperature flux evaluation. Thus there is concrete unnecessary representation
work, not merely a correlation between a large array and a slow routine.

However, the batch also includes field access, neighbor selection, flux
arithmetic, checks and output. The census counts logical accesses, not DRAM
traffic or cycles. No measurement assigns the whole 55 ms to record expansion.
Even eliminating the whole batch would account for only 8% of the 0.690 s gap.
This evidence is a useful component of a broader scalar-path redesign, not
justification for claiming a large total-step gain from a narrow getter change.

### 4. Adaptation costs are distributed, not dominated by the new allocator scopes

Source extraction costs 9.60 ms inclusive: interior copy/validation 1.58 ms,
boundary construction/copy 4.78 ms, ghost construction/copy 3.12 ms, and small
remaining work. Boundary plus ghost construction is about 82% of that scope.
These children must not be added to their parent.

Separate measured costs include migration route construction 10.23 ms,
serialization 7.59 ms, migration buffer allocation 3.15 ms and retained-state
preparation 1.54 ms. Catalog ownership balancing is 0.012 ms, local catalog
traversal about 0.001 ms, and descriptor/store allocation negligible. Source
field/topology allocation is about 0.14 ms in the exclusive ledger.

The older catalog residual remains 9.44 ms and migration residual 16.47 ms;
they are not resolved into a single avoidable mechanism. Migration MPI costs
36.24 ms. Stage 199 rules out several allocation/traversal hypotheses but does
not explain away the full 115 ms adaptation excess.

## What remains unresolved and the next test

Wavelet overhead remains systemic: inverse planning/packing/local copies/
installation/scaffold staging alone total 63.06 ms outside designated MPI
scopes; inverse posting/wait adds 94.58 ms. Separate Domain compatibility
publication/import/writeback scopes total 46.71 ms, and block ghosts 33.97 ms.
These are measured block costs, not all additional to legacy: generic Domain
boundary savings and differently placed numerical work must remain in the
signed comparison. Reducing message bytes alone does not remove representation
maintenance or producer delay.

Ordered batch records prove invocation order and nesting, not which remote
producer delays a particular receive. The analyzer's nearest-preceding-producer
links are chronological associations. They do not match send/receive edges and
cannot establish a network-versus-load-balance causal split. Nor do these
results overturn the Stage 198 finding that blocks have larger absolute, but
smaller relative, aggregate rank-compute spread than legacy.

The agreed goal is therefore **not yet fully met**: the measured time is
accounted for by location, with several concrete mechanisms identified, but a
large safe recovery estimate is still missing. Do not start another small
optimization campaign on the strength of these results alone.

The next bounded diagnostic should compare the unchanged legacy and candidate
at **21, 42 and 83 tasks**, using the same J5 checkpoint, physical integration
interval, workload/numerical gates and controlled CPU/socket placement. Start
with two bare timings per side/count and minimal parent/rank work measurements;
do not repeat all 27 observer variants. Add a paired one-socket versus
two-socket 83-task placement comparison to distinguish topology sensitivity
from rank-count sensitivity. All counts fit on one 128-core socket on bb13.

Compare absolute gap, rank-mean and rank-maximum local work, calls and work per
rank, designated MPI time and each implementation's scaling. Normalize compute
by measured work, while preserving total time, because fewer tasks also change
cache working sets and ownership granularity. Scaling alone is not a unique
causal decomposition; it should discriminate the current hypotheses and select
which complete path needs a dependency or hardware-counter experiment.

The current harness enforces 83 ranks. It must be parameterized and its
placement, completeness and correctness checks adapted before this sweep;
changing only `sbatch --ntasks` is not sufficient. Keep collecting one review
tarball containing results, manifests, placement and logs.

Machine-readable evidence: `parallel-blocks-cluster-199-54687.json`. Batch times
and counters therein are inclusive; never sum parent and child entries.
