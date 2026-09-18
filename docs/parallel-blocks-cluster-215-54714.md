# Stage215 job54714 — accepted performance result, pending promotion

**The candidate meets the 100 ms target.** Corrected ordinary-step medians show
117 ms saved per step at 83 tasks (12.57% less time; 1.144× speedup), confirmed by
both reversed pairs. Correctness gates pass. No solver branch or main was changed
by this review; recommend carrying this candidate forward as the next baseline.

|Run set|Baseline|Candidate|Time saved|
|---|---:|---:|---:|
|16-task screen, one pair|2.860s|2.570s|290ms (10.14%)|
|83-task early transfer, one pair|0.930s|0.814s|116ms (12.47%)|
|83-task confirmation, mean of two run medians|0.931s|0.814s|117ms (12.57%)|

The fresh confirmation order was A1=0.932s, B1=0.816s, B2=0.812s, A2=0.930s.
Adjacent-pair savings are **116ms and118ms**. Both exceed 100 ms; their mean is
outside the predeclared ±10 ms near-threshold band. One ABBA batch suffices under
the agreed rule. No extra measurements are needed to establish this campaign's
acceptance. Two runs per version are independent run-level observations; the
nine retained timesteps per run are not nine independent repetitions.

## Timing correction discovered during review

The frozen kit reported115.75ms saved (0.92725→0.81150s). Its parser incorrectly
consumed the restart-exclusion flag on the checkpoint transaction's timing row.
`main_mod::write_checkpoint` calls `restart` before returning to the outer step's
log print. Thus the first full integration step after reload is the following
row, not the checkpoint row itself. This was already identified in the earlier
evidence reconciliation, but the campaign inherited the old Stage210 parser.
That was a preparation oversight.

Raw logs confirm this ordering in all 12 runs. Correct selection excludes
**steps1,8,9**, using one-based indices: initial startup, checkpoint/save/reload
transaction, then the first full post-reload step. Nine ordinary steps remain,
including ordinary remaps. The corrected parser also handles field/data saving
independently and does not discard the step after a checkpoint with no reload.

`protocol_215.py` now preserves the restart flag across the checkpoint row.
Two regression tests cover the real nested ordering and checkpoint-without-reload;
all 13 protocol tests pass. The original kit, supplied review archive and recorded
results were preserved. The audit reproduces both their original decisions and
the corrected decisions. Every advance/stop decision and the executed12-run
schedule remains valid; a cluster rerun is unnecessary.

The printed timing precision contributes at most approximately1ms uncertainty
to each paired difference here: the rounding intervals are115–117ms and117–119ms.
This does not threaten the 100 ms result. These are rounding bounds, not statistical
confidence intervals or bounds on machine-to-machine performance.

## Verification

- All 1,236 archive inventory entries pass their byte counts and SHA-256 hashes.
- All three build inventories and compiler-log hashes match the frozen Stage215
  kit. Stage202 batching is on, Stage203 is off and FP contraction is off. Only
  the checked build has the route-audit macro. The build record reports the same
  23 external physics objects.
- All 12 run binary identities match their recorded builds. Raw logs reconstruct
  every timestep and the original stored medians. Timing runs have profiling and
  numerical oracles disabled; separate checked and oracle validations pass.
- All 12 stored same-count checkpoint comparison gates pass reassessment, with
  zero differing fields, topology/header differences or nonfinite values.
  This includes baseline self-comparisons as well as actual candidate comparisons.
- All 12 adaptive-workload comparisons reproduce from the logs. Both reloads
  (3 and4) and a remap after the second reload are present. The cross-count
  comparison passes its fixed tolerance; cross-count bitwise equality is not
  claimed.
- All 594 rank/run placement, resource and exit records agree with the prescribed
  CPU maps on bb13:16 tasks on one socket,83 split42/41 across two sockets.
  No swap-in/out is recorded during the timing runs. Confirmation peak rank RSS
  is approximately380–381MiB for both variants; no material peak-memory increase
  is apparent in these records. Sums of per-rank peaks are not simultaneous RSS.

The review archive deliberately omits checkpoint and executable payloads.
Consequently numerical acceptance was reassessed from the stored field reports
and binary identities, not by rehashing absent executables or rerunning checkpoint
comparisons. No new solver run was performed for this review.

## What this establishes

This is a useful example of the earlier successful pattern: amortize repeated
horizontal/address preparation across fields and layers, and pack directly into
retained storage. The native block plans were already cached; generic Domain
halo scans were a separate remaining path. Treating every route-cache proposal
as already exhausted would have missed it.

The117ms is the **net complete-step** benefit, including route construction,
retained storage and remaining communication. The candidate leaves MPI message
counts, payloads and numerical operations intact. There was no component ablation,
so do not assign117ms solely to packing, metadata scans, zeroing or reduced MPI
waiting. The result supports work removal without proving a detailed causal split.

Relative to the historically quoted385ms gap,117ms is about 30%; a subtraction
would leave roughly 268 ms. That is only contextual arithmetic: no contemporaneous
legacy timing was included, and the shared generic Domain path could also benefit
legacy if applied there. This campaign establishes faster blocks relative to the
accepted block baseline, not a newly measured block-versus-legacy gap.

Recommend retaining Stage215 as the next investigation baseline without another
broad profiling campaign. Any next candidate should be evaluated incrementally
against it. Do not re-count this saving or automatically infer that the remaining
Domain preparation can deliver another 100 ms. Main remains
`b82647b36eb825aed0c0aae4e92adcf4ab46d608`.

Archive SHA-256:
`4ce8982cac81835a980bc766c0cdb5adcf68e04575f32a20607ae07f6574367e`.
Machine-readable audit: `parallel-blocks-cluster-215-54714.json`.
Reproduce with `review_cluster_215.py --archive ARCHIVE --kit FROZEN_STAGE215_KIT
--out REVIEW_JSON`. Use the original frozen kit to reconstruct archived medians;
the repository parser contains the correction.
