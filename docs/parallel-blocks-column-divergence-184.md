# Stage 184 candidate: direct physical-column divergence

Prepared from Stage 183 checkpoint `99571e6129261c36804f41d98ebf1dc41bd4f626`.
This is a bounded first consumer of direct column access, not completion of the
column-stencil roadmap. Stage 183 full cluster acceptance remains pending.

## Change and hypothesis

`recompute_block_scalar_divergence_level` previously resolved eight logical
record values and wrote one logical value for each node at each physical layer.
It now resolves the four horizontal stencil locations once per node/scalar,
reads live strided columns, evaluates the same ordered six-edge expression,
and writes the resulting native divergence column directly. For 30 layers,
location calls fall from 240 to 4 per column. This is a work-count reduction
in this consumer, not a forecast of whole-step speedup.

The storage API checks column extents, variable boundaries and geometry classes
once per column operation. Compact field reads directly address the 17-slot
array; compact geometry is fetched once and broadcast across the physical
column. Oracle mode reads its independent full records at every layer, including
geometry. The old scalar divergence/accessor path remains available to oracle
comparisons and other consumers.

There is no persistent plan, pointer or numeric cache: addresses are local to
one consumer invocation after flux publication. The loop changes from
layer/node to node/layer; it writes only native divergence and reads finalized
flux/geometry, so no newly written divergence feeds this loop. Numerical
parenthesization is unchanged. Every consumed flux, including inactive-node
flux, retains finite/poison checks. Activity-based read skipping is deferred.
Soil/surface columns are excluded exactly as before. Native mass remains the
production mass path; the mass scalar path here runs only in oracle mode.

This changes only two production files, both required on a current Stage 183
remote tree:

- `src/parallel_block_scalar_storage.f90`
- `src/parallel_block_mpi.f90`

The existing profile-report reduction carries three added counters without a
new collective. `divergence columns: columns locations replaced scalar reads`
reports executed columns, column location calls, and the equivalent old scalar
read count. The last is computed from the executed layer count, not a measured
baseline counter. Existing restriction-kernel exclusive timers provide the
broader timing context. No production Domain bridge or oracle is added.

## Cluster timing to collect now

Yes: preserve an optimized **Stage 183** binary before replacing either source
file. Use the standard J7 restart fixture and **83 tasks**, retaining the usual
RK scheme, CPU placement and timestep interval. Build on the execution CPU or
with a compatible target (`-march=native` has previously caused cross-CPU SIGILL).

Run three repeats in fresh sibling output directories from the same immutable
seed, preferably within the same allocation. Do not advance the seed between
repeats. Keep full RK3/RK4 checked production/oracle restart acceptance separate.
For basic production timing, set:

```sh
export WAVETRISK_VALIDATE_BLOCK_DYNAMICS=0
export WAVETRISK_VALIDATE_BLOCK_ADAPTATION=0
export WAVETRISK_VALIDATE_BLOCK_REMAP=0
export WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS=0
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=0
export WAVETRISK_PROFILE_BLOCK_DETAIL=0
unset WAVETRISK_RESTART_PROBE WAVETRISK_PROBE_DOMAIN WAVETRISK_PHYSICS_PROBE_ALL WAVETRISK_ALLOCATION_CENSUS
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
```

In each already prepared fresh run directory, copy the preserved Stage 183
optimized executable as `climateJ5`, then use the established launch:

```sh
set -o pipefail
sha256sum climateJ5 simple.in
/usr/bin/time -p srun -n 83 --partition=bb ./climateJ5 simple.in 2>&1 | tee timing.log
```

Retain complete logs, input/checkpoint hashes, executable/source identities,
compiler flags, RK scheme, job/allocation identifiers, node/CPU placement and
exit status. Report per-step solver timing and whole-job elapsed separately;
startup, I/O and restart costs affect the latter. If using independent `srun`
allocations, report placement/load confounding. One separate optimized detailed
profile (both profile switches set to 1, other switches unchanged) is useful;
its timing is not interchangeable with the unprofiled repeats.

After candidate acceptance, alternate Stage 183/candidate runs in the same
allocation where possible. First require exact candidate-versus-block-baseline
checkpoint fields. The original unmodified legacy fixed-screen failure remains
visible: Stage 182 explains above-screen errors on the local fixture but does
not make that gate pass. Do not relax tolerances, bypass the strict legacy
timing campaign gate, or deploy diagnostic replay binaries. No clean speedup is
yet established. Authoritative main is unchanged; final-owner placement remains
deferred.

## Local validation

Results and identities are recorded in `parallel-blocks-column-divergence-184-results.json`. This
candidate's bounded local scope is the helper suite, an optimized J6 interval
through checkpoints 5/6 with exact Stage 183 comparison, checked RK4 first-step
production/oracle runs, and a separate first-step work-count profile. Full
checked RK3/RK4 cluster restart/remap coverage remains required before acceptance.

The initial source-only build lacked the untracked physics package's generated
`Makefile.inc` and stopped before solver compilation. The dependency file was
copied from the verified Stage 183 archive after confirming all physics `.F90`
sources matched, and the build resumed successfully. Solver optimized and
checked objects use separate fresh directories. The physics package retains
its existing `-std=f2003` build flags; solver DEBUG checks do not instrument that
external package. The pre-existing optimized air-temperature descriptor warning
remains; checked solver compilation passes `-Werror` and the 132-column limit.

### Results

- Optimized and checked RK4 solver builds pass. The 49-test helper suite passes;
  the final two focused storage tests additionally verify independent oracle
  geometry and rejection of invalid column extents/classes/writes.
- Optimized four-rank J6 execution through time 0.3350 completes both
  checkpoint/reload cycles and remaps. Compared fields, headers, thresholds and
  topology at checkpoints 5 and 6 match Stage 183 exactly.
- Checked RK4 first-step production and all-three-oracle cases pass. These do
  not establish full checked restart coverage or RK3 coverage for this candidate.
- The separate optimized first-step profile reports 122,240 columns, 488,960
  location calls, and 29,337,600 equivalent former scalar reads. The location
  reduction is 60-fold in this consumer. Removed velocity Domain-source and
  Domain-gradient calls, and producer full-oracle record counts, remain zero.
- Restriction kernels/replay (a broader region than this consumer) reports
  per-rank accumulated wall min/mean/max 1.1675/1.3848/1.5396 seconds. This single
  profile has no matched Stage 183 sample and establishes no runtime speedup.

Raw runs: `/private/tmp/wavetrisk-stage184-column-runs-02/`.
Builds: `/private/tmp/wavetrisk-stage184-column-opt-01/`.
Reproduction driver: `/private/tmp/wavetrisk-stage184-validate.py` (fresh output
directory required on reuse). The original checkout and authoritative main
were not modified. These results were recorded before the retained changes were committed.
