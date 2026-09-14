# Stage 189 cluster campaign 54549

The supplied campaign passes its numerical gates and shows a modest, consistently
positive timing difference on `nk`, node `bb13`, using 83 ranks in one allocation.
Compilation excerpts supplied on 14 September resolve the CPU flag question:
both optimized builds used `-march=native` and `-ffp-contract=off`; the checked
build used `-O0 -ffp-contract=off` with runtime checks.

## Correctness and identity

All 12 runs complete 12 timesteps, reload checkpoints 3 and 4, and report a
vertical remap after checkpoint 4. The input is the J5/J7 RK4 fixture with
`time_end=0.315`, `dt_write=0.01`, 30 atmospheric and 10 soil layers.

All 14 recorded gates pass: 12 independent run comparisons plus two reference
self-checks. The reports contain zero differences in headers, topology,
thresholds and all 16 field groups, with no nonfinite values. They include
checked oracle-on/off, optimized baseline/candidate, repeated/profile output,
and checked/optimized comparisons. These conclusions verify the supplied
comparison reports; raw checkpoints were not supplied for a new comparison.

Both source manifests match the local reference manifests. Every run binary
hash matches its recorded build role. The build record reports 23 unchanged
external physics objects, GCC 13.2.0, the MVAPICH 4.0 module and explicit
`FP_CONTRACT=off` for all three builds. No new legacy comparison was performed.

## Performance

Warmups, checked runs and profiles are excluded from the timing medians.

| Metric | Stage 188 baseline | Stage 189 candidate | Reduction |
| --- | ---: | ---: | ---: |
| Median campaign-run elapsed | 33.95563 s | 33.53563 s | 1.24% |
| Median of ordinary-step medians | 1.53 s | 1.50 s | 1.96% |
| Scalar restriction region, exclusive rank-mean wall sum | 1.688068 s | 1.517936 s | 10.08% |

All three candidate elapsed times (33.4920--33.7065 s) are below all three
baseline elapsed times (33.9151--33.9628 s). The median reduction is 0.4200 s.
This is consistent evidence within this short campaign, not a statistical bound
on future runs or a multi-node performance result.

The scalar restriction profile improves in both reporting windows. Its mean
self-CPU sum also decreases, supporting removal of local computational work.
These profile times are sums of rank-mean exclusive times, not a decomposition
of the solver critical path. Broad inverse kernels and ghost copying are
approximately unchanged in the supplied profiles.

The candidate reports 107,960 interpolation columns and 6,477,600 location
calls, replacing 194,328,000 equivalent scalar reads. Both windows satisfy
60 locations per column and 30 equivalent layer reads per location. Existing
native mass/velocity work counts match; removed velocity Domain source/gradient
counts and producer oracle record counts remain zero.

The isolated local interpolation gain was not a predicted whole-solver gain.
This campaign measures the combined effect in its actual workload. No cluster
memory/swap measurements or independent CPU-binding dump were supplied; the
recorded launch uses `--cpu-bind=cores` and a single exclusive node allocation.
Do not compare these elapsed times directly with earlier two-node `bb` results.

## CPU flags confirmed, 14 September 2026

The user supplied the expanded compilation commands for
`src/parallel_block_mpi.f90` from all three build logs:

| Build | Optimization | Native CPU targeting | FP contraction |
| --- | --- | --- | --- |
| Stage 188 baseline | `-O2` | `-march=native` | off |
| Stage 189 candidate | `-O2` | `-march=native` | off |
| Stage 189 checked | `-O0` with bounds/FPE/initialization checks | absent | off |

These excerpts agree with the recorded Makefile SHA256
`8c22b84bcbdc74089abfba8c8b78e46c24a0c331b35000ddec71a947d85b1f1a`.
The user subsequently confirmed that the edit removing native targeting was in
the main checkout's `wavetrisk_hydrostatic/Makefile`. The builder runs `make` in
the separate `$HOME/wt-stage189` worktree and verifies that worktree's Makefile
hash. The main-checkout edit therefore did not affect these builds. The checked
build's different optimization flags are intentional; its exact checkpoint
comparisons pass.

The comparison therefore uses matching explicit CPU-targeting and arithmetic
options for baseline and candidate. The build-host CPU model is not identified
by these excerpts; do not describe this as proof that compilation occurred on
`nk`. Full build logs were not supplied, so this confirmation is based on the
changed translation unit's expanded commands together with the existing build
manifest and executable hashes.

Recommendation: retain Stage 189 and commit the implementation and validation
records. No repeat campaign is needed solely to resolve this flag question.
Acceptance covers this single-node `nk` configuration and tested window; the
observed 1.24% elapsed reduction is not a universal or two-node speedup claim.

## Evidence

- `parallel-blocks-cluster-189-54549-campaign.json`: byte-preserved supplied results.
- `parallel-blocks-cluster-189-54549-analysis.json`: verified gates, file hashes,
  extracted profile regions and explicit provenance limitations.
- Raw copies of all four supplied files: `/private/tmp/wavetrisk-stage189-cluster-54549/`.

No production source or compiler settings were changed during this analysis.
