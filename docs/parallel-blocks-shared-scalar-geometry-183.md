# Stage 183: shared horizontal geometry in scalar execution storage

## Prerequisite and scope

Stage 181 is committed as `78ba3574`. Stage 182's controlled experiments
identified two float32 rounding-midpoint crossings in the physics interface
that explain the checkpoint errors above the original screening limits.
See `parallel-blocks-mixed-precision-diagnosis-182.md` and its JSON evidence.
The user authorized proceeding on that diagnosis while keeping the original
legacy-screen failure explicit. No numerical tolerance or production physics
input is changed here. **No replay code is linked into production.**

This stage replaces the storage behind the scalar kernels' existing logical
50-slot addressing interface, across interiors, boundaries, and ghosts:

- 17 values remain per scalar/vertical sample.
- 33 geometry/metadata values are stored per horizontal node, not per sample.
- A separate inactive-level geometry class preserves existing local
  soil/surface zeros versus transported geometry. This costs one additional
  33-value record per horizontal node, not another vertical/scalar expansion.
- Oracle runs keep independent full 50-slot records. Production neither
  allocates those records nor reconstructs them for normal physics capture.

It also removes repeated work: local producer geometry is installed once per
patch; remote interior geometry is no longer expanded over every field;
boundary/ghost geometry is installed once per node; and their field buffers
are initialized separately. Three-/four-value live transfers use contiguous
or indexed field accesses rather than rebuilding a 50-value record.

Numerical formulas, operation ordering within the kernels, ownership, MPI wire
contracts, adaptation/restart invalidation, and mass/velocity compatibility
are retained. The 50-slot extent now describes a **logical address space**,
not allocated bytes. Detailed memory output and the allocation census count
actual field/geometry/full-oracle allocations.

## Production files to transfer

All four files are required, including the new module and Makefile dependency:

```sh
scp /Users/kevlahan/wavetrisk_hydrostatic/Makefile bbserv:~/wav/
scp \
  /Users/kevlahan/wavetrisk_hydrostatic/src/parallel_block_scalar_storage.f90 \
  /Users/kevlahan/wavetrisk_hydrostatic/src/parallel_block_mpi.f90 \
  /Users/kevlahan/wavetrisk_hydrostatic/src/parallel_block_profile.f90 \
  bbserv:~/wav/src/
```

Neither `multi_level.f90` nor `time_integr.f90` changes in this stage. The
authoritative `main` branch/archive remains unchanged. The new stage is
committed as a locally validated checkpoint at the user's request; full
cluster acceptance remains pending. The earlier Stage 181 commit is separate.
See `parallel-blocks-handoff-183.md` for continuation and remaining work.

## Cluster test settings

Keep your existing RK3/RK4 selection and standard restart fixture. Use separate
fresh build directories when switching checked/optimized modes; copy the newly
built executable into the run directory before each run.

For **both RK3 and RK4 checked oracle runs** (`DEBUG=check`):

```sh
export WAVETRISK_VALIDATE_BLOCK_DYNAMICS=1
export WAVETRISK_VALIDATE_BLOCK_ADAPTATION=1
export WAVETRISK_VALIDATE_BLOCK_REMAP=1
export WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS=0
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=0
export WAVETRISK_PROFILE_BLOCK_DETAIL=0
unset WAVETRISK_RESTART_PROBE WAVETRISK_PROBE_DOMAIN WAVETRISK_PHYSICS_PROBE_ALL WAVETRISK_ALLOCATION_CENSUS
```

Also run **checked production** (`DEBUG=check`, particularly important here
because production uses compact storage and oracle mode uses full records):

```sh
export WAVETRISK_VALIDATE_BLOCK_DYNAMICS=0
export WAVETRISK_VALIDATE_BLOCK_ADAPTATION=0
export WAVETRISK_VALIDATE_BLOCK_REMAP=0
export WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS=0
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=0
export WAVETRISK_PROFILE_BLOCK_DETAIL=0
```

For **optimized production timing** (`DEBUG=false`), use those same six zero
settings. For a **separate optimized detailed profile**, keep all oracles and
detailed diagnostics zero but set:

```sh
export WAVETRISK_PROFILE_PARALLEL_BLOCKS=1
export WAVETRISK_PROFILE_BLOCK_DETAIL=1
```

Use the established launch pattern and distinct output names, e.g.:

```sh
cp ~/wav/bin/climateJ5 .
srun -n 83 --partition=bb ./climateJ5 simple.in 2>&1 | tee RK4_production.log
```

Do not interpret oracle-mode timing as compact-storage production timing.
Use the matched legacy/block experiment protocol for speedup claims; do not
compare independent cluster allocations without reporting placement/load
confounders. The original legacy checkpoint screen remains separately visible.

## Local validation

The final optimized build is
`/private/tmp/wavetrisk-stage183-storage-opt-02/bin/climate`.
Its full J6 run is `/private/tmp/wavetrisk-stage183-storage-run-02`:
both checkpoint/reload cycles and remaps completed, and **all compared
checkpoint fields, thresholds, topology, and headers match the block baseline
exactly** at checkpoints 5 and 6. An earlier compact implementation also passed
the same exact comparison before the remaining receive/ghost copy loops were
removed.

The storage module is compiled in a focused Fortran regression with bounds
checks, signalling-NaN initialization, and FPE traps. Tests cover patch versus
boundary/ghost indexing, two geometry classes, canonical installation,
independent scalar/vertical field updates, allocation reuse, empty storage,
release, remote geometry propagation, and capacity reduction. All 48 helper
tests pass. A clean `DEBUG=check` RK4 build and an isolated RK3 checked build
also compile successfully.

All four checked first-step cases completed successfully: RK3 production,
RK3 all-oracle, RK4 production, and RK4 all-oracle. These exercised bounds/FPE
checking in both compact and full-record storage modes. They are **first-step
checks**, not full checked restart/remap tests. The longer optimized production
and census runs cover both checkpoint/reload cycles and remaps. The checked
evidence is `/private/tmp/wavetrisk-stage183-checked-runs-01/checked-results.json`.

The allocation-census run also matches the final production checkpoint fields
exactly. At the same logical sample (42, `scalar-storage-ready`), summing the
four ranks gives:

| Counted allocated capacity | Stage 181 | Stage 183 |
| --- | ---: | ---: |
| Scalar workspace, including covered arrays/inline storage | 6.517 GiB | 2.329 GiB |
| Total across the four instrumented modules | 10.487 GiB | 6.300 GiB |

The reduction is **4.188 GiB**. The expanded `%full` scalar arrays have zero
allocated capacity throughout all sampled production phases. All four ranks'
55-sample sequences align. This is a comparison at a common logical phase,
not a synchronized wall-clock measurement. The census excludes other modules,
automatic temporaries, allocator/MPI/OS overhead, and unsampled transients;
these numbers are not total process memory or RSS.

Census evidence is
`/private/tmp/wavetrisk-stage183-storage-census-run-01/census-results.json`.
These are validation/diagnostic runs, not clean timing experiments. The local
production observations still show substantial system memory compression;
no local or cluster speedup is claimed from their wall times. Full cluster
checked production/oracle restart tests remain the next acceptance gate.

A persistent summary of binary/input identities, checkpoint comparisons,
allocation evidence, and all four checked cases is saved in
`parallel-blocks-shared-scalar-geometry-183-results.json` alongside this note.

## Tested source identity

The four files in the workspace match the final tested optimized archive.
Verify these SHA256 values after transfer:

```text
c52f7d17797fd6f34e279fc874f457f47aa4d1eca133dc6962850007d218bf2a  Makefile
79d0b06bbbd903b4170282ac82ef631d3a978aed4ab74b8c2996e5f393166bb9  src/parallel_block_scalar_storage.f90
d4cdab696d6d26341916d482bcec8b7ec5f6d17dc58c8807cd5dcec95fa2a852  src/parallel_block_mpi.f90
85317ae10bd3f71b2d5700ca9eddca62a12b3fca67aea0f3299a4a245d4cde93  src/parallel_block_profile.f90
```

The locally tested optimized binary SHA256 is
`912eea24082b2612d3a60456478140bf1df79b4f5934cc6c8bd6c9231d19fb2a`.
Cluster-built executable hashes will differ. Production contains no causal
physics replay; those isolated Stage 182 executables must not be deployed.
