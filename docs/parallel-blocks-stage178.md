# Stage 178: inverse transport and stencil execution overhead

Base: accepted Stage 177, commit `70176e4a`. Stage 178 accepted after
cluster RK3/RK4 debug oracle and optimized RK4 restart tests passed.

Cluster production total: 19.858 s; profiling total: 20.404 s. First-eight-step
inverse rank-average time: 3.4363 s (kernels 0.69136 s). The historical Stage
177 comparison is not a controlled repeated timing experiment. All RK4 printed
state/mass records agree across debug, production and profiling runs; cache
and final-owner oracle mismatch counters are zero for both integrators.

## Changes

1. Compile block-side inverse transport addresses once per topology generation.
   Pack/unpack directly between current block payloads and the persistent wire
   buffer. This removes per-node `sample`/`reshape` temporaries and moves
   family/storage dispatch outside the physical-layer/scalar loops. Addresses
   are integers, not pointers retained across RK payload swaps.
2. Combine the outer-vector patch-interior and boundary scatters into one
   peer exchange and one completion. They write disjoint storage families;
   order within each family/peer is unchanged. Both ends derive the combined
   subplans from existing metadata without an additional setup collective.
3. Post alias receives/sends after packing the remote-source snapshot and
   before ordered local copies. Complete communication before installing
   remote aliases. Local copies touch only the native workspace, not the
   in-flight buffers. MPI progress/overlap benefits depend on the MPI runtime.
4. Cache off-patch stencil addresses used by scalar lifting/reconstruction and
   inner-velocity interpolation. Previously those kernels repeatedly searched
   patch/boundary/ghost tables for each field, layer and geometry lookup.
   The cache resolves each encountered stencil coordinate once per generation.
   It does not cache numerical geometry, masks or field values. Changing
   masks/RK payloads therefore remain visible. Remap addressing is excluded.

The fourth change removes integer computation from the numerical kernels;
it does not reduce the number of physical interpolation operations or change
their arithmetic. The first change removes repeated packing/indexing work.
This is not a claim that all remaining inverse communication or Domain
compatibility has been removed, or a prediction of a particular speedup.

## Preserved contracts

- No level/phase dependency is skipped or reordered. In particular, all-level
  local alias copies, fixed-coarse resets, plus-side/pentagon writes and ghost
  refreshes remain in place.
- Remote alias inputs are snapshotted before local copies, and remote install
  permutations remain unchanged. No values are consumed before MPI completion.
- The same native path runs with and without the oracle. Independent Domain
  phase checks and final tendency checks remain; no tolerance is relaxed.
- Generation changes after adaptation, ownership changes and restart rebuild
  transfer routes and stencil caches. Neither cache contains numeric payloads.
- Mass/velocity/temperature production chains from Stages 175–177 are unchanged.
  Their removed Domain producer counters must remain zero in production.

## What the profile should show

The first-eight-step 83-task Stage 177 reference has 496 gather calls,
568 alias calls and 568 scatter calls. This stage combines one pair of scatter
calls for each outer-vector level transition. Gather/alias call counts and
all three byte totals should remain unchanged on identical inputs/partition:

| Native transfer | Stage 177 bytes |
|---|---:|
| Gather | 4,835,316,480 |
| Aliases | 7,694,283,360 |
| Scatter | 19,475,907,360 |

Scatter messages and completion points should decrease. The new line
`inverse stencil addresses: resolved reused = ...` counts actual off-patch
searches and searches avoided while native inverse execution is active.
Counters are reduced only at the existing optional profile report, not at
each inverse call. Watch kernel, gather/scatter, alias and setup times as well
as total dynamics time. These timers are inclusive; do not sum overlapping
rows. Additional integer route/cache metadata and the combined exchange
buffer have a memory cost; there are no per-phase buffer allocations.

Stage 177 cluster production time was 21.284 s; first-eight-step inverse
rank-average time was 3.8343 s, including 0.84914 s local kernels and
2.0146 s boundary synchronization. Use repeated paired runs on the same
allocation for speed claims. Four-rank local gains do not establish gains
on the 83-task cluster configuration.

## Transfer and build

Only these three source files change from Stage 177:

```sh
scp src/parallel_block.f90 src/parallel_block_inverse.f90 \
  src/parallel_block_mpi.f90 bbserv:~/wav/src/
```

Compile on a CPU compatible with the execution partition (`DEBUG=false`
uses `-march=native`). Use fresh build directories and serial builds.
Select RK3 in `src/shared.f90` for the first build, then restore RK4.

```sh
# timeint_type = "RK3"
make -j1 DEBUG=check PARAM=param_J5 TEST_CASE=climate \
  BUILD_DIR=build/stage178-rk3-check BIN_DIR=bin/stage178-rk3-check
# timeint_type = "RK4"
make -j1 DEBUG=check PARAM=param_J5 TEST_CASE=climate \
  BUILD_DIR=build/stage178-rk4-check BIN_DIR=bin/stage178-rk4-check
make -j1 DEBUG=false PARAM=param_J5 TEST_CASE=climate \
  BUILD_DIR=build/stage178-rk4-opt BIN_DIR=bin/stage178-rk4-opt
```

## Complete test settings

From the established cluster test directory, use the identical checkpoint-3
restart input for each run. Copy the selected binary only after a successful
build; the commands above do not update `bin/climateJ5`.

### RK3 checked oracle

```sh
cp ~/wav/bin/stage178-rk3-check/climate ./climateJ5
export WAVETRISK_VALIDATE_BLOCK_DYNAMICS=1
export WAVETRISK_VALIDATE_BLOCK_ADAPTATION=1
export WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS=0
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=0
set -o pipefail
srun -n 83 ./climateJ5 simple.in 2>&1 | tee RK3_debug_oracle.log
```

### RK4 checked oracle

```sh
cp ~/wav/bin/stage178-rk4-check/climate ./climateJ5
export WAVETRISK_VALIDATE_BLOCK_DYNAMICS=1
export WAVETRISK_VALIDATE_BLOCK_ADAPTATION=1
export WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS=0
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=0
set -o pipefail
srun -n 83 ./climateJ5 simple.in 2>&1 | tee RK4_debug_oracle.log
```

### RK4 optimized production timing

```sh
cp ~/wav/bin/stage178-rk4-opt/climate ./climateJ5
export WAVETRISK_VALIDATE_BLOCK_DYNAMICS=0
export WAVETRISK_VALIDATE_BLOCK_ADAPTATION=0
export WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS=0
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=0
set -o pipefail
srun -n 83 ./climateJ5 simple.in 2>&1 | tee RK4_production.log
```

### RK4 optimized profiling

```sh
cp ~/wav/bin/stage178-rk4-opt/climate ./climateJ5
export WAVETRISK_VALIDATE_BLOCK_DYNAMICS=0
export WAVETRISK_VALIDATE_BLOCK_ADAPTATION=0
export WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS=0
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=1
set -o pipefail
srun -n 83 ./climateJ5 simple.in 2>&1 | tee RK4_profile.log
```

## Local validation

Test root: `/private/tmp/wavetrisk-stage178.eoSNHa`. The standard user test
directory and `bin/climateJ5` are untouched. Local checks use four MPI ranks.

- Fresh serial RK3/RK4 checked builds and optimized RK4 build: PASS, normal
  strict 132-column limit retained. Checked builds have bounds checks, FPE
  traps and signaling-NaN initialization. The pre-existing optimized
  `air_temperature` descriptor warning remains; no warning suppression added.
- `sh test/parallel_block_inverse/run.sh`: PASS at checked `-O0` and `-O2`.
  Bulk pack/install is compared bit-for-bit with the original node routine
  for both scalar/vector and solution/wavelet families, mixed patch/boundary
  routes, nonzero field-level offsets, multiple blocks and empty plans.
- RK4 debug oracle: PASS, one restart timestep/all four substages plus final
  inverse transform. Printed state/mass records agree with Stage 177.
  Stencil searches: 57,678 resolved, 34,605,522 reused. Debug tests overlap
  and their timings are not a performance comparison.
- RK3 debug oracle: PASS, one restart timestep/all three substages plus final
  inverse transform. Printed state/mass records agree with Stage 177.
  Stencil searches: 57,678 resolved, 25,939,722 reused. All cache/final
  shadow reports in both oracle tests have zero mismatches.
- Optimized RK4 non-oracle: PASS, ten timesteps including vertical remapping,
  checkpoint-4 write/reload, adaptation and two post-restart timesteps.
  All printed numerical state/mass records match Stage 177. Checkpoint 4 is
  byte-identical, SHA-256:
  `6086e1fab875369095b1dea9800af7a39261ab751d6f10916718ed3c081179f3`.
  Domain mass-flux/velocity compatibility rows remain absent; native velocity
  Domain-source/gradient counts remain zero in both profile windows.

First-eight-step local work comparison (four ranks):

| Metric | Stage 177 | Stage 178 |
|---|---:|---:|
| Gather calls/messages/bytes | 496 / 806 / 427119360 | 496 / 806 / 427119360 |
| Alias calls/messages/bytes | 568 / 5568 / 1692104640 | 568 / 5568 / 1692104640 |
| Scatter calls/messages/bytes | 568 / 1605 / 1972149120 | 496 / 1426 / 1972149120 |

There are 475,386 resolved and 285,527,334 reused inverse stencil addresses
in this window (99.83% reuse). The numerical interpolation operation count
does not change. Local work counts are not the 83-rank cluster counts.

Sequential same-fixture optimized/profile comparison, Stage 178 followed by
the committed Stage 177 binary, with no builds or other test runs overlapping:

| Reported time (s) | Stage 177 | Stage 178 |
|---|---:|---:|
| Full ten-step total | 306.15 | 270.64 |
| First-eight-step inverse, rank-average | 62.851 | 51.214 |
| Inverse local kernels | 31.722 | 25.577 |
| Inverse setup/import | 9.7250 | 8.5442 |
| Native gather | 4.5125 | 3.6271 |
| Native aliases | 4.7269 | 4.5403 |
| Native scatter | 8.8352 | 6.8609 |

Total reported time is 11.6% lower; inverse time is 18.5% lower in this
single local pair. It is not a repeated, order-balanced benchmark and does
not establish the 83-task benefit or attribute every timing difference to
these changes. Historical Stage 177 local timing was different again;
do not mix that older run with this pair to claim a controlled result.
The fresh baseline and candidate have identical printed numerical records
and byte-identical checkpoint 4. Inclusive rows above must not be summed.

Local binary hashes:

```text
c7f4d566dd4aea9a5470b76440fd5622c408e959a33ca74a9b09f6eb1d363ffb  bin/stage178-rk3-check/climate
c6ea4ae1090277ac3323452aa568c0989c7dbb74eef4baa5611da0c67c5a24f8  bin/stage178-rk4-check/climate
c6091757cff79826e1f94451fde06a6dfcbc4c33534b896fc8a8864c802492ad  bin/stage178-rk4-opt/climate
```

Logs under the test root:

```text
rk3-oracle/RK3_debug_oracle.log
rk4-oracle/RK4_debug_oracle.log
production/RK4_optimized_restart_profile.log
baseline177/RK4_baseline177_profile.log
```
