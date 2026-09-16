# Stage 199: resolve producer and representation costs before optimization

Stage 199 preserves Stage 197 numerical sources, the Stage 196 historical control,
and the pinned legacy Domain executable. Stage 198 job 54667 is the preceding
measurement reference. There is no new production optimization or arithmetic
change. The question is which avoidable operations determine completion time in
the 0.689 s/step bare gap, rather than which routine has the largest positive timer.

## Four measurements

1. Split `basic_operators` into the matched pressure/step1 patch producers, scalar
   physics capture, level divergence capture, Bernoulli and Exner restriction,
   native mass restriction, and mass exchange. Patch producers retain their
   existing order. Timers surround complete patch or level operations, not nodes
   or scalar accesses. Legacy and block producer labels do not assert identical
   mathematical work: block step1 includes native mass production.
2. Record ordered batches with step, context, parent sequence, layer, grid level,
   position/component, logical work and rank-local start/end. Mass post/pack and
   wait/install have exact `(k,l,pos)` keys; scalar operations spanning all layers
   use `k=0`. Scalar restriction/divergence and ghost/edge exchanges carry their
   grid level. The capture-finalization exchange consumes all layer/level producers
   and uses `(0,0)`; do not attribute this fan-in to one arbitrarily chosen layer.
   Existing MPI/post/wait self-times are recorded as a subset of each batch.
3. Split adaptation into catalog traversal/allocation/ownership selection, source
   allocation, interior copying and validation, boundary/ghost construction,
   stencil construction, migration sizing/routing/serialization, receive allocation,
   local installation, retained-state preparation and staging release. Source-block
   extraction records the leaf-node ownership weight as its work count. Required
   validation remains in its original production locations and is timed explicitly
   with its containing construction/copy operation.
4. Remove sampled clocks from scalar read/write/address functions. A separate
   census binary performs guarded integer counting only. Enclosing batches record
   the seven access-count deltas; the exclusive operation ledger remains separate.
   Per-access nanosecond estimates and isolated kernel replay are not used.

Fine patch-producer and allocation timers contribute to the exclusive ledger but
are omitted from the event tape to limit memory and event overhead. Larger batches
retain aggregate storage counts. Thirty-five new operation IDs extend the ledger
from 101 to 136 entries per context. Hooks restore the exact preceding observer
source bytes when removed; numerical patches retain the Stage 197 hash.

## Imbalance hypothesis and its limits

Block ownership uses `subtree_weight_Domain`: a leaf contributes `PATCH_SIZE**2`
and an internal patch contributes the sum of its children. It does not explicitly
weight internal-level traversal, boundary size, route count, compatibility copies,
or scalar storage lookups. Equal leaf weights can therefore hide unequal overhead.
This is a source-supported hypothesis, not a demonstrated explanation of the full
gap. Both the shared and pinned legacy Domain load proxy count nodes and edges whose
absolute masks exceed ADJZONE (`comm_mod::domain_load`). The block leaf-node
proxy therefore differs from the Domain active-mask proxy. Neither proxy alone
establishes which rank is late in the measured workload.

Stage 198 showed a wavelet maximum-minus-mean non-MPI time of approximately 93 ms
for blocks and 39 ms for legacy. Tendency values were much closer: 148 and 131 ms.
Relative to mean CPU work, those spreads are about 49% versus 76% for wavelets
and 40% versus 66% for tendency (blocks versus legacy). Blocks are more balanced
on this aggregate relative measure; their much larger CPU workload produces larger
absolute delays. Do not describe this as demonstrated deterioration of load balance. Negative compute/wait
correlation is expected when ranks synchronize; it cannot by itself separate
arrival imbalance, MPI progress, memory effects and transport cost.

The source also contains a specific repeated-operation candidate: with
`mass_transaction` active and multiple grid levels, `basic_operators` exchanges
native mass at edges, and its caller exchanges mass at the same `(k,l)` again
before native divergence. The ordered 108/109 records will show the number and
cost of these occurrences and any intervening batches. This observation is not
permission to remove an exchange; data consumers, intervening work and the
critical-path contribution must be established first.

The batch analyzer retains each rank's measurements, computes per-key work/time
correlations, and identifies the nearest compatible preceding producer batch.
These links are chronological candidates, not proof of a dependency. The raw parent
and sequence records remain available for source-guided analysis, including
multi-producer fan-in. Clocks have independent origins; never subtract a timestamp
on one rank from a timestamp on another. Inclusive batch rows overlap and must not
be summed or added to the disjoint full-step ledger. Independent per-operation
maxima are not a completion path or a recoverable-time estimate.

## Calibration, acceptance and collection

Use seven fresh builds and 27 short runs on the same exclusive nk allocation,
with 83 ranks, J5/J7, RK4, contraction off and identical physics objects. Two
rounds reverse order. Preserve the previous checkpoint/workload gates, including
checked production/oracle runs, exact same-implementation comparison and the fixed
cross-implementation legacy screen. The latter is not bitwise legacy equivalence.

The passes are bare, minimal common clock, operation timers, timers plus counters
and ordered batches, and the separate storage-census binary. Its control runs in
mode 2 (identical batch/timer collection with access census disabled); its active
runs use mode 3. Comparing those two isolates added census overhead. Comparing the
census control with the ordinary mode-2 observer measures binary-layout effects.
Per-region calibration flags remain mandatory. The campaign reports locations and
leaves `optimization_ready=false`; source/consumer audit must establish avoidable
work and a conservative completion-time benefit without double-counting MPI waits.

`ranking.csv`/`ranking.md` retain the signed full-step ledger. `batches-summary.json`
contains rank distributions and chronological links. `batches199-rank-*.txt` are
raw ordered records; normal operation, dependency and common-clock records remain
available. The collector includes these files and build provenance in one tarball.
Records are buffered during timesteps and written after the run; no barriers or
per-event I/O are introduced. Missing ranks, inconsistent scopes, corrupt records,
missing production coverage and event-buffer overflow fail closed.

Follow [the cluster commands](parallel-blocks-cluster-199.md). Local validation
covers J4/J6 restart correctness and instrumentation; cluster performance claims
require the bb13 campaign.

## Recommended scaling follow-up

A paired strong-scaling study at 21, 42 and 83 tasks on bb13 would test how the
legacy/block ratio changes as work per rank decreases. Hold J5/J7 inputs and
endpoints fixed, use deliberate comparable placement across both sockets, and
record actual affinities and memory evidence. The present standard driver still
requires 83 tasks; changing only its command-line rank argument is not supported.

Compare per-rank mean and maximum non-MPI work, operation counts, MPI/post/wait
intervals and total step time for both implementations. A roughly persistent
CPU excess per matched unit of work supports extra computation/representation
cost; increasing communication fractions at higher counts support a scaling
limit, but repartitioning, boundary growth and cache/NUMA changes can confound it.
Normalize imbalance by mean CPU work as well as reporting milliseconds. Use
within-count exact restart/reference checks and the fixed legacy screen; do not
assume bitwise equivalence across different MPI task counts. This is a proposed
follow-up experiment, not a result or an implemented rank-sweep runner.
