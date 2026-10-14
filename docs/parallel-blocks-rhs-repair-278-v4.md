Current reviewed result: [kit-v4 job55501](parallel-blocks-stage278-job55501-results.md) passes science/ownership and improves traffic/time, but still fails absolute performance, memory and scaling acceptance. Stage266 remains accepted.

Historical preparation record follows; its pending-campaign statements describe the pre-submission state.

# Stage278 kit-v4: dependency, round and context repair

Prepared 2026-10-13. Stage266 remains accepted; this is an unpromoted repair of
kit-v3's measured performance failure. No broader ownership migration or partition
change is stacked onto it. [Authoritative scorecard](/Users/kevlahan/.codex/worktrees/26d8/wavetrisk_hydrostatic/docs/parallel-blocks-scorecard.md).

## Evidence motivating the changes

Job55500 passed science and connected RHS/RK ownership, but candidate ordinary
steps were 46–56% slower than266 and 109–143% slower than main at83/160. At256,
derived exchanges contributed105,795,480 sends (99.70% of recorded connected
stage messages) and52.91GiB (85.82% of its payload). Mean derived payload was537
bytes. Summed rank peak RSS was118.76GiB versus43.07GiB for266; these are sums of
individual peaks, not a simultaneous allocation peak. These counters locate a
structural target; they do not establish critical-path attribution or identify
the cause of MPI/internal memory growth. [Verified review](parallel-blocks-stage278-results.md).

## Implemented changes and lifetime contract

1. `parallel_block_rhs_routes.f90` replaces full-program ghost refreshes with a
   field-specific union of numerical readers and writers across levels. It keeps
   canonical owned slots coherent for later fallback providers, active gradient
   and restriction inputs, direct/SW/pentagon writes, signed aliases and current
   final writers. This is conservative cross-level demand, not minimal per-phase
   liveness. Compilation occurs at topology/mask plan creation; execution retains
   direct fused arithmetic and compact precompiled exchange indices.
2. `parallel_block_rhs_dynamics.f90` sends Bernoulli and Exner together before
   injection and together after injection: four completed halo transactions become
   two per level/layer. Both fields are ready at each boundary and use the same
   canonical producer identities. Original numerical/injection ordering remains.
   A reusable two-column buffer costs16 bytes per retained node; no30-layer
   replication or runtime allocation is introduced.
3. `parallel_block_rhs_distributed.f90` owns one isolated communicator per live
   executor plan; all derived routes borrow it using distinct level/phase/channel
   tags. Previously23 communicator contexts were held per active grid level.
   The new schedule holds21 routes per level but one shared context in total.
   Requests complete before returning; all routes are released before their
   context is freed on rebuild/teardown. Separate executors cannot share tags
   accidentally because their contexts remain isolated.
4. `parallel_block_rhs_exchange.f90` distinguishes borrowed and owned contexts,
   validates the portable tag range and rank-consistent signature, and releases
   only owned contexts. Legacy callers retain their existing independent contexts.
   Input packing precedes signed local alias writes, preserving snapshot semantics.
5. `parallel_block_rhs_backend.f90` includes the added thermodynamic buffer in
   the component memory ledger. Compiler temporaries are scoped to preparation;
   plans and buffers are reused within an epoch and released before replacement.

The five-file [patch](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/repair-v4/solver-changes.patch) is relative to sealed kit-v3.
The complete candidate still contains the explicit23-file overlay over266.
No claim is made that all retained Domain/inverse/physics workspaces are removed.

## Qualification and quantitative results

Fresh optimized, checked and independent RHS/RK oracle replays pass at4 ranks:
12 steps, 2 reloads, 2 remaps, exact accepted-family checkpoint/step states and
62 exact VTK members each. These are source-bound to all251 build-source files
and23 overlay files. Current P>D ownership must be rechecked in the campaign;
the four-rank J4 fixture has40 Domains and cannot establish that gate.

The8 new communicator tests cover1/4 ranks, checked/optimized modes,8 rebuild
epochs, repeated payload exchanges, distinct tags, repeated release and a surviving
borrower, plus invalid/rank-inconsistent tag rejection. Existing34 exchange
invocations also pass, including large-message rendezvous and protocol failures.
A fresh serial J5 checked build and link passed on bbserv:72 compilation commands
and1 link with GCC13.2/MVAPICH4, unchanged23 external physics objects. This was
compiler qualification only, not a Slurm campaign or J5 runtime test.

For the same4-rank J4 validation fixture:

| Metric | Prior278 v3 | Repair v4 | Change |
|---|---:|---:|---:|
| Derived payload bytes | 4,068,589,440 | 3,301,517,760 | -18.85% |
| Derived send messages | 781,440 | 697,200 | -10.78% |
| All recorded stage send messages | 783,264 | 699,024 | -10.75% |
| Maximum route payload allocation, MiB | 109.55 | 109.22 | -0.31% |
| Maximum summed component memory, MiB | 437.286 | 437.367 | +0.018% |

The new buffer offsets route storage savings on the largest local rank. There is
**no demonstrated overall memory reduction** from this local ledger; MPI internal
storage and compilation peaks are excluded. Cluster RSS and preparation must be
measured. Input/RK adapter traffic and primitive event counts are unchanged.
[Raw costs](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/repair-v4/cost-comparison.json).

Bare complete12-step runs use identical fixtures, matched external physics objects,
4 MPI ranks, oracles/audits off and a forward/reverse order266→v3→v4→v4→v3→266.
Timesteps are not counted as independent repetitions.

| Source | Two run medians, seconds | Mean of run medians, seconds |
|---|---:|---:|
| stage266 | 5.504607, 5.108444 | 5.306526 |
| prior278 | 4.948468, 5.068136 | 5.008302 |
| candidate | 4.937739, 5.029329 | 4.983534 |

Candidate change is -0.49% versus prior278
and -6.09% versus266 on this local screen.
The 0.49% difference from prior278 is smaller than the observed run spread;
no repeatable speedup is established. Mean complete-trajectory time is70.332s
versus68.214s for prior278 (+3.11%) and69.503s for266 (+1.19%); the ordinary-step
comparison therefore does not demonstrate a whole-trajectory gain.
Full trajectory and launch times for every run are preserved in the
[bare comparison](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/repair-v4/bare-comparison/results.json). This screen cannot forecast
J5 performance or establish MPI scaling. Prior local near-parity failed to predict
the cluster regression; the same limitation applies here. In particular, the
observed message/payload reductions do not justify promising the32–36% J5 runtime
cut required to recover266, or the52–59% needed to reach main at83/160.

## Shared-memory assessment

`MPI_Comm_split_type(..., MPI_COMM_TYPE_SHARED, ...)` can create a group whose
processes can share memory; it does not itself create shared storage or eliminate
a numerical dependency. A direct-load path needs `MPI_Win_allocate_shared` and
correct producer/consumer visibility and completion synchronization. `MPI_Win_sync`
is a memory ordering operation, not by itself a peer notification protocol.
See the [MPI communicator specification](https://www.mpi-forum.org/docs/mpi-5.0/mpi50-report/node190.htm),
[shared window allocation](https://www.mpi-forum.org/docs/mpi-5.0/mpi50-report/node313.htm)
and [window synchronization](https://www.mpi-forum.org/docs/mpi-4.1/mpi41-report/node331.htm).

Such a path may reduce within-node matching/copy costs, but it cannot fix excessive
rounds, uneven block work or inter-node scaling. No measured comparison currently
supports its addition here. Keep it as a separate transport A/B after measuring
this reduced schedule: shared-communicator/window lifetime once per topology plan,
local first-touch, immutable packed publication buffers, explicit epoch completion,
portable inter-node MPI, and within-NUMA/cross-NUMA/inter-node comparisons. Avoid
node-wide barriers per phase or shared placement that concentrates traffic on one
NUMA node. The current bb13 test is one node, so any shared-memory gain there would
need independent multi-node validation before a scaling claim.

## Next campaign and acceptance

The sealed kit-v4 retains the proven exclusive nk/bb13,256-task,mem=0,one-hour
allocation:4 serial builds,8 validation runs and8 bare timing runs. It directly
compares repair v4 to freshly built266 and main at83/160, and266 at256. Historical
v3 timings are contextual, not a same-allocation A/B control. Keep exact science,
checked/oracle, all96 zero-Domain ranks doing authoritative work at256, retired
Domain entry guards, preparation, component memory and rank RSS checks.

The decision is whether this specific reduction in dependencies/rounds and context
multiplicity recovers absolute performance and resource cost. If still slower,
retain266, quantify the remaining deficit and repair the remaining dominant
schedule before another migration. Fewer messages alone do not pass performance.
Main parity and useful >=640 execution remain separate unmet dual-goal gates;
256 is intermediate, and the640-block ceiling/imbalance remain unchanged.
No new Slurm campaign has been submitted. User launch/status/collection:

```bash
bash ~/rk278-connected-rhs-kit-v4/launch.sh
```

```bash
bash ~/rk278-connected-rhs-kit-v4/status.sh
```

```bash
bash ~/rk278-connected-rhs-kit-v4/collect.sh
```

Reproducibility evidence is under `/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/repair-v4`. The original incremental checked-build
failure (stale transitive module interfaces after exported derived types changed)
is preserved; fresh builds passed. A test-only exact-real warning was corrected
using exact bit comparisons. Neither failure changed the numerical gates.

## Delivery and optional local scaling decision

Sealed archive: `/Users/kevlahan/Desktop/rk278-connected-rhs-kit-v4.tar.gz`. SHA256:
`dc95a4d3723376efb63b9fac94be24701eb0f4f8d8cae066acd6c685d9e6763f`. Fresh extraction, immutable sources, packaged grid,
12 campaign helper tests and bbserv fixture preflight pass. Staged remotely at
`/home/kevlahan/rk278-connected-rhs-kit-v4`; no Slurm submission.

A fixed-J4 1/2/3/4-task comparison of266, prior278 and repair278 is worthwhile
to distinguish fixed overhead, memory-bandwidth effects and within-machine
scaling, given the inconclusive four-task repair gain. Record ordinary-step and
complete-trajectory times, per-rank memory, speedup T(1)/T(P) and efficiency
T(1)/(P*T(P)), consistent placement and exact science. Repeat only when observed
variation changes the decision; keep run order balanced. This matrix has not
been run. All counts are below J4 D=40, so it cannot test zero-Domain execution
or replace J5 and inter-node qualification.
