# Stage250: calibrate candidate storage logging before further Domain removal

Job 55240 has completed. The verified [cluster review](parallel-blocks-stage250-review.md)
finds no measured savings from suppressing the diagnostic write and a
1.123003 s/step quiet-candidate gap to `main` at 83 ranks. The following is
the preparation protocol used for that run.

Stage249 job 55239 measured 1.677614 s/ordinary step for the current
mirror-free candidate versus 0.540682 s for read-only `main` at 83 ranks on
bb13, a 1.136932 s end-to-end gap. The candidate's timed `time_step` calls
`assert_domain_prognostic_absent`, which scans Domain field descriptors and
writes one diagnostic line per rank per step. Main does not. This job measures
the write's cost with the assertion retained, and then remeasures the quiet
candidate's gap to main. The investigation takes precedence over any further
Domain removal or promotion.

The sealed source patch changes only `src/domain.f90`: it caches the value of
`WAVETRISK_QUIET_STORAGE_ABSENCE` once per process and conditionally omits the
diagnostic write. The storage assertion still runs. Loud (`0`) and quiet (`1`)
runs use the **same patched binary**, same full 12-step J5 restart, same bb13
83-rank physical-core map and same external physics objects. Main is the
unchanged Stage249 read-only export with Stage190 checkpoint I/O hardening.

The campaign starts with one full validation run each for main, loud candidate
and quiet candidate. It requires the established adaptive trajectory,
checkpoint/VTK reproducibility and numerical gates. It then runs two
loud-quiet-quiet-loud batches (four independent full runs per candidate mode)
to estimate logging savings, followed by one main-quiet-quiet-main batch (two
independent full runs per version) to measure the residual regression.
Ordinary steps are 2–7 and 10–12; each full run contributes one median.
There is no automatic promotion.

The prior Stage248 feature-off baseline was 1.629455 s/step at 83 ranks,
whereas the Stage230 block candidate was 0.754893 s/step. That cross-campaign
contrast suggests the large cost may be shared by early owner-column,
primitive, diffusion or owner-tendency features used in both Stage248
variants; it does **not** identify a culprit. After logging calibration,
narrowly ablate those shared feature groups on the same bb13 fixture and
placement, then profile the residual critical-rank work and MPI wait. Do not
assume the Stage229 provisional phase breakdown still applies to this source.

Packaging and both J5 source builds passed locally with `make -j1`, respecting
the local four-core limit. The pinned patch applied to the frozen candidate,
both packaged builds compiled with FP contraction disabled, and they shared
23 identical external Simple-physics objects. Source/package hashes, Python
and shell syntax, fixture-rejection preflight and fake-Slurm fresh extraction
passed. These local checks are build and workflow evidence, not speed evidence.
The [kit README](../test/parallel_block_profile/logging_calibration_250/README.md)
has the no-edit transfer, launch and collection commands. Cluster timing is
pending; `main` is unchanged.
