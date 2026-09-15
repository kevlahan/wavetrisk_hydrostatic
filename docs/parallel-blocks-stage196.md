# Stage 196: attribute dynamics excess to boundary mechanisms and computation

Stage 195 is committed as `5dfe0db1`. It places 41.87% of the mean excess in
complete dynamics WT_after_step transactions and 41.22% elsewhere in dynamics.
Its envelopes cannot establish whether boundary handling or numerical work
is the primary cause inside either part. The new stage adds observations only;
production algorithms, arithmetic order, selector defaults and checkpoint
comparison tolerances are unchanged.

## Boundary architecture

Legacy routes field boundaries through the general `update_bdry` interfaces
and their start/finish implementations in `comm_mpi.f90`. These already have
separate scan/pack, receive routing, allocation/zeroing, MPI posting, local copy,
wait and installation phases. The abstraction itself does not make these free:
its repeated field/route scans and allocation are real work.

Blocks retain these Domain routines for shared inputs and consumers. They also
use reusable specialized paths, rather than one independent implementation in
every numerical routine: `exchange_block_scalar_restriction_ghosts`, SOL ghost
refreshes, inverse owner transport, scalar input transport, native mass boundary
plans, and Domain writeback. These paths implement different storage layouts,
ownership and dependency order. Several have their own packing and request
loops. Sharing more machinery may improve maintainability, but cannot be assumed
to improve performance without measuring the traversals and synchronization.

The source audit also found that the older fine-region timers did not isolate
all reverse-import work. Some full import scans/packing and boundary imports
with MPI_Alltoallv could remain in wavelet self time. Stage 196 explicitly
brackets those bodies and their collectives rather than treating that residual
as computation.

A concrete round trip occurs in native wavelet compression: Domain-owned
coefficients are imported to compact blocks, compressed, and published back to
Domains before inverse reconstruction. Source inspection establishes the extra
representation work; it does not quantify how much is safely removable.

## New measurements

Every sample belongs to the same rank, timestep and common Stage 195 dynamics
envelope. There are three caller contexts: inside WT_after_step, inside the
implementation-specific tendency-preparation tag, and the remaining dynamics.
For direct comparison, the last two are combined as **outside WT**. Thus the
two reported excess ledgers cover the full 83% target without adding a nested
tendency timer a second time.

Within each context, a bounded stack records **exclusive wall time**, call counts
and logical work counters for up to 77 named operation scopes. The same clock
read marks each transition. Root/unclassified self time is retained. A context
change splits any open parent scope correctly; selected child time is subtracted
from the parent by construction. All per-rank context totals must reconcile to
the independent common WT/dynamics envelopes before results are accepted.

Boundary-only observation separates:

* General Domain scan/pack, route scan, allocation/zero, local copies, install,
  MPI posting and waiting.
* Inverse transfer planning, packing, aliases, posting/self-copy, wait and install.
* Block ghost packing/install and MPI; scalar transport/install and MPI.
* Native mass begin/finish exchange, including new preparation/post/install
  scopes around the complete routines, with nested waiting separated.
* Domain writeback preparation/install and MPI, plus its publication/commit remainder.
* Full Domain-to-block import preflight/install, packing/zeroing/local traversal
  and MPI_Alltoallv; boundary-only import preflight, zero/pack, collectives and
  installation/local copies.

Full operation observation additionally resolves scalar input/geometry/capture,
inverse scaffold staging, inverse numerical kernels, scalar restriction and
traversal, wavelet/compression self work, shared producers and basic operators,
native mass/velocity, RK assembly/retention and representation refresh. Timers
around numerical loops include traversal/setup within those scopes. Mixed and
unclassified work is labelled explicitly; it must not be renamed pure arithmetic.
The exact disjoint group map is in `analyze_operations_196.py` and the kit's
`operation-instrumentation.json`.

Existing boundary work counters are captured even with the old profiler off:
field count, route candidates scanned, doubles packed/received/zeroed, requests,
geometry samples and inverse plan build/reuse counts. New native mass counters
record send/receive buffer doubles, local-copy entries and requests. Full and
boundary imports additionally count buffer doubles zeroed, patches installed
and local/remote boundary entries. Existing
block transfer message/byte records are captured by caller context and profile
ID. These are logical source counters, not measured DRAM/network traffic or a
complete cross-implementation byte census. General boundary calls count start
and finish scopes; they are not unique complete halo exchanges.

## Avoiding a repeat of the profiling ambiguity

The cluster campaign retains untouched block and legacy references. The same
instrumented binary on each side runs in three modes:

1. Minimal: common workflow observation only; new operation clocks/counters off.
2. Coarse: selected boundary/transfer operations and logical work counters.
3. Detailed: all selected operation scopes and counters.

The old detailed and block profilers remain off in all three modes. The new
observer rejects attempts to combine them. This removes the old process-CPU
clock calls, profile reductions and detailed-report collection from the new
operation measurements. The additional hooks can still perturb execution;
that is measured rather than assumed negligible.

Each of six rounds rotates measurement modes and alternates block/legacy order.
The summary reports overhead both for the ordinary whole-step wall statistic
and separately for WT and outside-WT envelopes. A large or asymmetric observer
effect invalidates a literal decomposition of the bare gap; simplify the
instrumentation before drawing causal conclusions. No uniform correction is
applied to operation times.

All records are buffered until after the simulation loop. No measurement MPI
calls, barriers or file writes are inserted into timed dynamics. The source
transformer proves that removing its additions restores every input numerical
source byte. Existing profile wrappers have different disabled-path guards, so
Stage 196 adds independent paired hooks before those guards; it does not enable
the old profiling path just to obtain nested regions.

## Acceptance and interpretation

Five builds preserve untouched references, optimized observers and a checked
block observer, with common external physics objects and contraction off. The
58-run campaign retains all 146 Stage 195 exact/workload/fixed-legacy gates,
83-rank affinity records and swap evidence. Each diagnostic run must match its
own untouched implementation exactly at the checkpoint. Every rank must have
complete operation records, valid nesting, counters and conserved context time. Required boundary, import
and kernel scopes must have coverage; full-mode WT scope calls must equal the
independent common WT call count. This detects silently disabled hooks.
The subtraction of paired operation ledgers must also equal the matched common
dynamics excess.

Interpret signed differences, including negative legacy offsets. Determine
whether boundary/transfer preparation, waits or named numerical work dominates
**the excess**, not merely the block runtime. Inspect counts per field/route
alongside times to distinguish extra work from slower execution of the work.
A named MPI wait can include late peer arrival caused by computation; it is
not isolated network service time. There are no hardware performance counters
or globally synchronized arrival timestamps in this portable measurement kit.

Residual wavelet work still needs finer observation if it remains large. This
stage supplies measured locations and boundary mechanisms, not a promise that
all remaining work will be explained by one category or that any particular
ownership redesign will recover the whole excess.

Use [the cluster procedure](parallel-blocks-cluster-196.md) for the build,
58-run campaign and the agreed single-tarball log collection process. Local
validation evidence and final kit hashes are recorded separately.

## Completed local validation

The final kit sources byte-match the locally built block and legacy source
inventories. Five four-rank J4/J6 cases passed: optimized block and legacy in
both boundary-only and full-operation modes, plus a checked block build in
full-operation mode. Every case completed 12 steps, checkpoint 2 -> 3/reload
and two remaps, conserved every rank/context, covered the required scopes and
matched its Stage 195 implementation reference exactly. The checked case also
matched the optimized block reference; the fixed legacy numerical screen passed.
Local checked validation used bounds checks and floating-point traps with the
numerical oracles off. The cluster campaign retains checked-off and checked-on.

The final unit/integration suite ran 162 tests: 157 passed and five optional
cases were skipped. Prototype records predate the final import scopes and are
excluded from acceptance. Local runs are correctness/coverage tests, some of
which overlapped; their timings are not evidence about the cluster gap.

See `parallel-blocks-stage196-local-results.json` for comparisons, source/build
identities, scopes and work counts; `parallel-blocks-stage196-kit.json` records
the final package hashes. Cluster review is recorded below.

## Cluster review — job 54659

The 83-rank cluster campaign and archive audit passed. See
[the review](parallel-blocks-cluster-196-54659.md) for conserved signed
operation accounting and calibration. Boundary/ownership transfer routines
account for 73.2% of dynamics WT excess and 44.1% outside WT. Stage 196
remains an observer-only change; its implementation and cluster evidence are
ready to preserve together as the Stage 196 commit.

Correction to the work-counter coverage described above: generic boundary
scan/zero counters (12–15) and local geometry/sample counters (4–5) stayed
behind old-profiler call-site guards. They are unavailable, not zero work.
Timing coverage remains valid. The frozen kit is unchanged; a future
measurement revision must repair these call sites and enforce coverage.
