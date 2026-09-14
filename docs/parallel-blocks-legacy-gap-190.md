# Why Stage 190 is still slower than legacy

Stage 190 is committed as `cdd2eb64`. Campaign 54553 measures 1.43 s versus
0.558 s ordinary-step medians: 2.56272 times legacy, or 0.872 s more per step.
The evidence points to costs of a mixed Domain/block implementation, including
reconstruction, data transfers, retained validation work and synchronization.
It does not support an inherent 2.56x penalty for block decomposition.

The entire excess cannot yet be assigned causally: the campaign contains
block profiles but no legacy phase profile. Below, measured block phase costs
are distinguished from source-confirmed extra work and unmeasured hypotheses.

## Ordinary-step attribution

The new `analyze_stage190_gap.py` uses the existing tested detail parser. All
12 steps have complete 83-rank records; the checkpoint/restart step (8) is
excluded. Remap steps (5 and 10) remain ordinary, matching the timing harness.
The maximum discrepancy between the sum of exclusive regions and a rank's
whole-step time is 1.564 microseconds. Region groups below are disjoint.

These are means over the 11 profiled non-checkpoint steps and ranks, not the
unprofiled timing medians used for the 2.56x ratio. They sum to 1.40714 s.

| Block execution region group | Seconds/step | Share of profiled step |
| --- | ---: | ---: |
| Retained Domain boundary handling | 0.28835 | 20.5% |
| Native inverse, including its transfers | 0.25287 | 18.0% |
| Scalar input, restriction and associated MPI | 0.20194 | 14.4% |
| Block reconstruction, installation and migration | 0.17812 | 12.7% |
| Native mass/velocity, including mass MPI | 0.09983 | 7.1% |
| Wavelet/compression work | 0.09297 | 6.6% |
| Block ghost pack/install/MPI | 0.08285 | 5.9% |
| Shared Domain producers | 0.08215 | 5.8% |
| Domain/block publication | 0.04946 | 3.5% |
| RK and tendency assembly | 0.03472 | 2.5% |
| Other timestep work | 0.04389 | 3.1% |

These rows are NOT all overhead relative to legacy. Some perform necessary
mathematics or boundary work that legacy also performs. In particular, do not
subtract a whole inverse or boundary row from the gap and claim that saving.

## 1. The optimized run does not execute two complete tendencies

`prepare_block_multistage_tendency` in `src/time_integr.f90:139` runs
`block_tendency_compatibility_ml` to produce the shared inputs. The complete
Domain tendency callback is guarded by `validate_oracle`. The assertion in
`assert_multistage_legacy_tendency_calls` rejects any complete legacy tendency
call in optimized RK3/RK4.

The production compatibility pass in `src/multi_level.f90:230` still computes
shared primitives, pressure, diffusion inputs, mass and velocity work. But the
old Domain scalar tendency, velocity source and velocity gradient comparisons
are guarded by validation. The cluster producer oracle-record counts and the
Domain-source/Domain-gradient counters are zero. The named oracle profile
region has zero time. Its absence alone would not prove the absence of all
checks: the serialization checks below are separate and still active.

There is limited arithmetic duplication. `RK_sub_step_compatibility` forms
mass/velocity RK results in Domain-indexed native workspaces and copies them
into compatibility fields; `begin_local_block_accumulated_tendency_trial` also
forms the authoritative block RK state. The ordinary-step RK compatibility
region costs only 0.00836 s, and native RK assembly 0.00136 s. Even these whole
regions are much smaller than the 0.872 s observed median difference. Removing
another complete legacy RHS is therefore not the opportunity here.

## 2. Block reconstruction and validation still run every timestep

`src/main.f90:780` opens the grid-change transaction after RK, synchronizing
block state into Domain storage for physics/adaptation. After Domain adaptation
and inverse reconstruction, `src/main.f90:832` calls catalog construction,
source-block extraction, block migration/installation and retained-state setup.
This occurs each block timestep, not only on checkpoint steps.

`build_source_blocks` deallocates the previous source store and extracts every
current source-local block. `install_local_blocks` clears the installed store,
allocates it, and performs Fortran deep copies of source/received blocks.
Plans and workspaces are then rebuilt or revalidated against the new generation.
This construction phase is additional to the retained Domain adaptation path.
Its ordinary-step group costs 0.17812 s, including 0.03652 s of named migration
MPI. It excludes other block ghost and publication work charged elsewhere.

There are three definite production validation costs:

- `src/parallel_block_build.f90:260` unconditionally invokes migrating-source
  pack/unpack/pack comparison.
- `src/parallel_block_mpi.f90:5613` invokes received-block storage verification
  with serialization enabled, before the `run_full_validation` branch.
- `src/parallel_block.f90:692` and `:729` pack both sides of local deep copies
  and compare the byte streams.

`full_validation=.false.` therefore does not turn these proofs off. They cost
0.04419 s per ordinary step in exclusive rank-mean time, within the 0.17812 s
construction row. Across the complete reporting windows, the counters cover
11.33 GB source proof payloads, 11.33 GB received proof payloads, and 13.86 GB
copy proof payloads. Those are counted payload sizes, not total memory traffic;
round-trip packing and comparison touch the data several times.

This is a concrete first removal target, but the direct proof cost alone is
only about 3.1% of the profiled step. It cannot explain or eliminate the whole
gap. Keep cheap shape/identity/ownership checks in production and retain full
round-trip proofs in explicit validation runs.

For a larger improvement, retain unchanged block subtrees and plans through
adaptation. The test changes its grid almost every step, so merely skipping
reconstruction when the entire global grid is unchanged is unlikely to help
much. Incremental reuse must operate below the whole-grid level and refresh
changed physics fields without retaining stale numerical inputs.

## 3. Two ownership models create extra communication and copying

The current code retains Domain owners for the shared physics/primitive pass
and Domain-indexed native mass/velocity workspaces, while authoritative compact
blocks have separately balanced final owners. This creates a real bridge:

- Scalar capture transports Domain-produced inputs to block owners through
  `MPI_Alltoallv` in `src/parallel_block_mpi.f90:21739` and related branches.
- Scalar restriction has its own block ghost/edge exchanges.
- Native inverse has gather, alias/boundary, and scatter phases between these
  representations; see `src/parallel_block_inverse.f90:361`.
- Compatibility publication supplies Domain consumers again, and adaptation
  triggers another Domain-to-block reconstruction.

The retained shared producers also continue to use Domain boundary exchanges.
Across both whole profile windows, 2.22065 s of rank-mean boundary time is
attributed to the shared-producer caller, plus 0.46928 s to its basic operators.
These cross-cut caller totals are already included in the phase timings and
must not be added again.

Native mass communication loses one piece of legacy overlap. In the legacy
`trend_ml`, a tendency exchange is started before velocity-source work and
finished at the next coarse level. Production `exchange_native_mass` in
`src/multi_level.f90:692` posts, performs local copies, and waits before returning;
its AT_NODE exchange occurs at the next-level consumer. There is no equivalent
velocity-source overlap inside that schedule. Named native mass wait costs
0.07825 s per ordinary step. This is a justified scheduling target, not a claim
that all that wait would disappear with an asynchronous API.

The broader change is to reduce Domain/block round trips and execute more
shared producers on the owners of their consumers. Optimize the bridge before
assuming that faster scalar arithmetic alone can reach legacy performance.

## 4. Load balancing does not reflect the mixed workload

The catalog weights come from `subtree_weight_Domain`, which counts leaf-patch
nodes (`src/domain.f90:424`). The block-owner assignment does not model measured
kernel time, retained Domain work, halo cost, packing, or cross-owner transfers.
Rebalancing blocks does not move the existing Domain-side producer workload.

Named MPI regions account for 0.49604 s per ordinary step on average (35.3%).
They include MPI spin and arrival imbalance, not just data transmission. The
same 83 ranks remain bound to distinct cores throughout the campaign, with no
swap activity.

| Rank | Mean whole step | In named MPI regions | Outside named MPI regions |
| --- | ---: | ---: | ---: |
| 52 | 1.40701 s | 0.24381 s | 1.16320 s |
| 82 | 1.40694 s | 0.86501 s | 0.54193 s |

The time outside named MPI regions differs by over 2x while whole-step times
are almost equal. This establishes a substantial difference in work/wait
exposure across ranks. It does not identify every delayed dependency or prove
that all differences come from the weight model. Outside-named-MPI time also
contains checks, copying and potentially unclassified MPI; it is not a pure
useful-computation measurement.

A latency diagnosis must therefore include phase-by-rank costs and the separate
Domain and block ownership maps. Pinning everything to one socket may alter
NUMA/communication behavior, but it does not remove the duplicated storage or
ownership bridge and is not the primary remedy indicated by these sources.

## 5. The storage API still imposes avoidable work

The compact scalar execution workspace retains a record-oriented layout:
`field(17,samples)` in `src/parallel_block_scalar_storage.f90:87`, behind an API
that decodes the former 50-slot record addresses. A fixed field across adjacent
samples has stride 17 doubles; across a physical column the increment is
17 times the sample stride. `scalar_read_column` gathers those values into its
result. Legacy scalar fields provide direct arrays for each field/layer.

Stage 190 removes repeated geometry/address calculations from the restriction
closure, but other kernels still use checked address resolution and temporary
record/column values. Column batching does not by itself make the underlying
column data contiguous. This is a specific implementation tradeoff worth
measuring, not evidence that block storage is inherently worse. The scalar
replay arithmetic row itself is 0.09657 s per ordinary step; the larger
0.20194 s scalar group includes setup, capture/installation and MPI.

Domain boundary buffers also zero their full allocated capacity rather than
only the active receive prefix (`src/comm_mpi.f90:1087`). Whole-window counters
record 208.8 billion doubles zeroed versus 2.014 billion doubles received, a
ratio of about 104. Boundary storage/zeroing costs 0.02387 s per ordinary step.
This is concrete wasted work, but the routine is shared with legacy, so it is
not established as block-specific excess. Both comparators must receive any
common-boundary optimization in subsequent comparisons.

## Recommended next experiments

1. Isolate serialization/deep-copy proofs behind explicit validation and run
   the unchanged exact checkpoint gates. Measure this small, well-defined
   ablation separately. Do not remove numerical safety checks wholesale.
2. Profile legacy with the same region accounting and common I/O fix, retaining
   the authoritative unprofiled legacy binary for timing and checking that the
   instrumented output matches it. An older instrumentation helper exists, but
   it assumes an unpatched legacy archive and needs adaptation for this kit.
3. Test production lifetime/ownership changes one at a time: reuse unchanged
   subtrees/plans, restore mass-exchange overlap, and reduce the Domain/block
   bridge. Keep source identities, work trajectories, exact block gates and
   fixed legacy tolerances. Do not change tolerances to accommodate a speedup.
4. Use the resulting paired profiles to choose between moving shared producers
   to block owners and a field-oriented scalar workspace. A further small
   restriction-only optimization is unlikely to close the dominant gap.

No new numerical implementation was changed during this investigation. The
Stage 190 commit remains the validated baseline for these experiments.

## Reproduction

Run the helper with the preserved candidate profile and campaign results:

```sh
python3 test/parallel_block_profile/analyze_stage190_gap.py \
  /private/tmp/wavetrisk-stage190-cluster-54553/run_candidate.log \
  /private/tmp/wavetrisk-stage190-cluster-54553/results.json
```

The compact output is recorded as `parallel-blocks-legacy-gap-190-analysis.json`,
including file hashes, all rank averages, group membership and limitations.
