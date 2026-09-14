# Performance targets after Stage 188

Prepared 11 September 2026 against `61509429`. This is a proposed implementation
sequence, not a new timing result. The first interpolation implementation is
described in [Stage 189](parallel-blocks-restriction-columns-189.md). Stage 188 is
the comparison baseline; compile both sides with `FP_CONTRACT=off` in fresh
build directories. Preserve the common checkpoint I/O fix.

## Evidence and priority

The retained Stage 187 optimized J4/J6 four-rank profile identifies the following
regions. Values are rank-mean exclusive wall seconds summed across its two
reporting windows, not maximum-rank elapsed time or a cluster speedup forecast.

| Region | Seconds | Implication |
| --- | ---: | --- |
| Scalar replay | 25.23274 | First target: restriction execution and repeated scalar access |
| Scalar/inner inverse kernels | 9.27143 | Audit arithmetic versus gathers and temporary storage |
| Shared producer residual | 9.00194 | Later audit of necessary primitives and representation transfers |
| Block ghost pack/install | 7.37026 | Target repeated traversal and copies |
| Boundary buffer allocation/zero | 4.61096 | Audit initialized ranges and buffer lifetimes |

Source: `/private/tmp/wavetrisk-stage187-timing-01/candidate-profile/run.log`.
The profile used compiler-default contraction and experienced local memory
pressure. Refresh the matched component profile with contraction off before
accepting an optimization. These broad regions do not isolate the cost of any
individual proposed edit. In particular, the historical name "scalar replay"
includes production restriction execution; it is not evidence of an enabled
validation oracle or removable duplicate computation.

Stage 187's 83-task improvement and Stage 188's negligible observed arithmetic
policy penalty are established separately in their campaign reports. There is
no current controlled cluster ratio against legacy under the new policy.

## First candidate: direct scalar restriction columns

Stage 184 already converted divergence to column access. Extend that approach
to the remaining restriction chain in `src/parallel_block_mpi.f90`:
`restrict_block_scalar_flux_level`, `interpolate_block_scalar_flux`,
`complete_block_coarse_flux`, and their partial/coarse-flux consumers.

Today the physical-layer loop encloses horizontal-node work. Each interpolation
reconstructs restriction indices and repeatedly calls
`block_scalar_record_value`, which resolves patch/boundary/ghost storage,
computes a logical record index, and decodes compact storage. Parent edge masks
are also read for the activity maximum and again for each edge. The shared
horizontal geometry makes repeated layer-by-layer resolution a concrete source
of removable work.

Proposed changes:

1. Resolve a bounded horizontal stencil and its storage classes once per
   parent/child node invocation; gather physical-column data using the existing
   checked column storage API. Start with invocation-local addresses, avoiding
   a persistent cache or additional MPI collective.
2. Reuse invariant restriction indices, weights and parent edge masks across
   physical layers in compact production storage. Keep live flux values tied to
   the current restriction level and publication epoch. Do not assume oracle
   full-record geometry is shared across layers.
3. Evaluate interpolation, partial flux and coarse flux with the original
   ordered additions, subtractions and products. Preserve fine-to-coarse level
   ordering and audit aliasing before moving the node/layer loops. Store each
   completed parent flux column directly for the existing divergence consumer.
4. Preserve finite/poison checks for consumed values, soil/surface exclusions,
   independent oracle records, boundary and remote-node handling. Validate
   invariant addresses once; do not remove numerical validity checks to gain
   speed. Keep the original scalar implementation as the independent reference.
5. Count resolved stencils, columns and replaced scalar lookups through the
   existing profiling reduction. Separate restriction from divergence time if
   needed to attribute the gain; measure setup and scratch capacity as well.

For a 30-layer column, horizontal setup could be amortized over 30 layers.
This is a work-removal hypothesis, not a 30-fold kernel or solver speedup claim.

Expected production scope: `parallel_block_mpi.f90`, the scalar storage module
only if the existing column API needs an extension, and profiling definitions
only for targeted attribution. Avoid changing ownership or wire formats in
this candidate.

## Following candidates

**Inverse and ghost data movement.** Audit `parallel_block_inverse.f90` and the
inverse/ghost paths in `parallel_block_mpi.f90`. Use the existing compiled
addresses to batch column gathers and installs and eliminate intermediate
copies whose final consumers can use native storage. Stage 178 already caches
stencil addresses: do not recreate that optimization. Measure local copying,
allocation and MPI wait separately. Preserve boundary alias ordering and
generation invalidation after adaptation, remap and restart. Choose changes
from actual copy counts and profiles before altering MPI placement.

**Boundary initialization.** The three instrumented receive-buffer paths in
`comm_mpi.f90` already reuse capacity when possible. Audit the full-buffer zero
fills and exact-size reallocations; initialize only ranges that require defined
zeros, retaining all receive and accumulation semantics. A poison-filled-buffer
test must establish that every other consumed element is overwritten. Capacity
reuse must not increase retained peak memory enough to erase the gain. This
may also benefit legacy, so report common improvements separately from changes
that close the block/legacy gap.

**Shared producer interfaces and memory.** Inventory the remaining PV, pressure,
KE, Bernoulli, Exner and physics consumers before removing any Domain bridge.
Replace one complete expensive producer-to-consumer chain, then release storage
after its final consumer has moved. Shared primitive computation performed once
is necessary work, not automatically duplication. Use current live-allocation
and process-memory evidence; the large geometry compression in Stage 183 is
already implemented. Do not pursue another full ownership rewrite or checkpoint
codec change without evidence that it dominates the intended workload.

## Acceptance for each candidate

- Freeze Stage 188 source/binary identities and identical arithmetic flags.
  Run meaningful component comparisons against the preserved implementation,
  including boundary/ghost storage and invalid/poisoned inputs.
- Local development: evolved J4, max level 6, four ranks, the established
  checkpoint/reload plus subsequent remap window. Compare every saved checkpoint
  exactly with Stage 188; run checked production and independent dynamics,
  adaptation and remap oracles. Cover RK3 and RK4 where the changed path is used.
- Retain J5 restart coverage and the unchanged legacy comparison screen, with
  both builds using contraction off. Local legacy screen acceptance does not
  establish bitwise equality or new 83-task legacy equivalence.
- Measure optimized baseline/candidate sequentially in alternating order,
  excluding warmups, correctness runs and detailed profiling. Retain VM
  compression, swap and memory evidence. Reject a speed claim if variability or
  memory pressure prevents attribution; retain an experiment only when its
  component evidence and numerical gates justify it.
- Cluster acceptance: J5 build, max level 7, 83 tasks, immutable checkpoint 3,
  `time_end=0.315`, `dt_write=0.01`, expecting new checkpoint 4 and a subsequent
  remap. Use the ordinary baseline/candidate exact-gate campaign with three
  timing pairs in one allocation. The Stage 188 arithmetic diagnostic mode has
  a different purpose and must not substitute for this candidate gate.
- Keep checkpoint costs and ordinary-step timing separate. Retain the short
  case for iteration; use a longer interval on the same case if a small effect
  cannot be resolved. Do not reduce J5 to J4 on 83 tasks for the primary timing
  comparison.

Implement and accept the restriction candidate before combining these targets.
No new performance percentage is promised by this plan.
