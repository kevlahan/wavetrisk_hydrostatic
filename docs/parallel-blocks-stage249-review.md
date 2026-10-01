# Stage249 direct J5 comparison with main — job 55239

The direct bb13 comparison completed and confirms a large regression in the
current migrated executable relative to the read-only `main` export. At 83
ranks, the mean of two independent ordinary-step run medians was **0.540682 s**
for main and **1.677614 s** for the migration: **+1.136932 s per step**, or
**3.103 times** main's elapsed time. The A-B-B-A run medians were 0.542203,
1.686057, 1.669172 and 0.539162 s. Adjacent-pair regressions were +210.96%
and +209.59%. This is well beyond the 3% review guard; no second A-B-B-A
batch was needed and there is no promotion.

The 16-rank validation pair was 1.478569 s for main and 5.740728 s for the
candidate. It is a single pair used for correctness coverage, not a separate
confirmed performance result. Main's current 83-rank result is close to the
0.542883 s legacy result in Stage229; the current candidate is close to its
1.674537 s result in Stage248 v7. These cross-campaign observations support
stability of the broad comparison, but the Stage249 A-B-B-A is the accepted
current gap. The last pre-migration block candidate measured 0.758268 s in
Stage229, a 215.386 ms gap to legacy; Stage230 measured 0.754893 s for its
post-cleanup candidate against the prior block implementation, without a new
legacy run. Thus the current candidate's roughly 0.923 s increase relative to
Stage230 is a lead for source investigation, not a paired causal estimate.

## Correctness and provenance

All eight runs completed the 12-step J5 restart with two remaps and 62 VTK
members per run. The 16- and 83-rank placements were verified on bb13, with
83 ranks split 42/41 across its two sockets and one physical core per rank.
Every run has the same adaptive state sequence. Main's checkpoint hash is
stable across all its runs; the candidate's checkpoint hash is stable across
all its runs and matches the previously qualified Stage248 v7 J5 hash. All
four main-to-candidate checkpoint comparisons pass the established numerical
policy, including identical topology and thresholds; the atmospheric values
are numerically equivalent within tolerance, not bitwise identical. The
candidate's storage absence checks pass on all ranks.

Both binaries were built on bb13 using serial compilation, the same 23
external Simple-physics objects and FP contraction disabled. The main source
was a read-only export of `6a245db7f1632b0205c2bebe95e5f0cd535c2582`
with only the established Stage190 checkpoint I/O hardening. The pinned J5
input, seed and grid hashes matched prior campaigns. The review archive
contains 1,091 files and verifies at SHA-256
`2faf3b9e2664c43d8f230754d42909ebddd379e268a2b7ca5348e21cfd0c9295`.
The independent [reviewer](../test/parallel_block_profile/review_cluster_249.py)
recomputed the run medians and paired result and checked archive hashes,
placements, rank exits, numerical reports, storage logs and run order. Its
machine-readable output is [here](parallel-blocks-stage249-review.json).
The return archive omits checkpoint, VTK and binary payloads, so those payload
hashes and numerical reports were checked as recorded rather than regenerated.

## Interpretation and next step

The campaign has no phase profiler, but it is not free of candidate-side
validation work. `time_step` calls `assert_domain_prognostic_absent` in the
mirror-free build; that assertion scans Domain field descriptors and writes
one line per rank per step. The archived logs confirm 83 such lines for every
timed step in each candidate run, and none in main. This is part of the
current executable's measured cost. The confirmation proves its current
end-to-end regression, but cannot allocate the gap between the RK migration,
other block operations, storage assertions, MPI waiting and logging. Whole-run
rank user CPU was also higher in the candidate; it is not phase attribution.

Before another broad removal, calibrate the timed assertion/logging cost on
bb13 with the same J5 83-rank placement, source, inputs and A-B-B-A run
protocol, suppressing only diagnostic output while retaining the storage
assertion. Then localize the residual post-Stage230 cost with bounded phase
and work-count evidence or narrow feature ablations. Keep main unchanged and
retain the full bbserv results directory until this review is accepted.
