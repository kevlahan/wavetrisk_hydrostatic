# Stage230 cluster review — job 55089

Stage230 passes qualification as an architectural cleanup, with no demonstrated
material speedup. All 15 runs completed on bb13 (nk). The accepted comparison is
the committed Stage229 block implementation versus the seven Stage230 source
changes, with STAGE202/203 enabled in both and FP contraction disabled.
This campaign did not run legacy, so it does not remeasure the legacy gap.

## Correctness and storage

- Checked production and full dynamics/adaptation/remap oracles pass at 16 and
  167 tasks, with 160 geometric Domains. Optimized production also passes at
  167 tasks, covering ranks without owned Domains.
- Every run completes 12 matching step records, reloads 3 and 4, and two remaps.
- Same-rank checkpoint reports are exact; all 62 VTK member hashes match each
  same-rank baseline. Cross-rank numerical gates pass as well.
- Production and checked candidate logs report zero Domain tendency-storage
  allocations on all ranks, on both initializations. Oracle runs retain that
  storage on all ranks, as required by the independent reference path.
- Source/build manifests match the sealed kit; the three builds record identical
  23 external physics objects. The archived build logs and run binary identities
  match the recorded hashes. Affinity evidence verifies all actual rank placements.

This qualifies removal of `trend` for the restarted RK4/Simple production path.
Other configurations and independent oracle paths deliberately retain it.
It does not establish that Domain `sol`/`q1`, primitive/diffusion scratch,
geometry-indexed native workspaces or wavelet storage have been removed.

## Ordinary-step performance

| Comparison | Baseline (s/dt) | Candidate (s/dt) | Observed reduction |
| --- | ---: | ---: | ---: |
| 16-task screen, one pair | 2.334464 | 2.317022 | 17.44 ms / 0.75% |
| 83-task transfer, one pair | 0.760555 | 0.753845 | 6.71 ms / 0.88% |
| 83-task A–B–B–A, mean of two run medians/version | 0.757403 | 0.754893 | 2.51 ms / 0.33% |

The confirmation order was A1=0.761783088, B1=0.755185961,
B2=0.754599438, A2=0.753023279 seconds. Adjacent-pair savings were +6.60 ms
and -1.58 ms. Their opposite signs preclude treating the small average reduction
as a demonstrated speedup. Both comparisons are comfortably within the 3%
regression guard. The automated `disagreement=false` refers to agreement about
that regression threshold, not agreement about the sign of a speedup.
There is no minimum gain requirement for this architectural deletion, and no
reason to repeat this campaign solely to resolve a few milliseconds.

The reviewer recomputed high-precision ordinary medians from each raw `test_log`
and stdout event order. Steps 1, 8 and 9 are excluded: startup, saving/checkpoint,
and first post-reload step. Nine ordinary steps remain per run. Timesteps are not
independent repetitions. Checked/oracle timings and single-pair screens are not
included in the A–B–B–A mean. The 167-task production validation was 0.582699 s
versus 0.577331 s for baseline (+0.93%); that is exploratory correctness-run
context, not a replicated high-rank performance conclusion.

## Evidence and next step

All 2,590 inventory files pass hash verification. Archive SHA256:
`100257973da5f70b29a3fccb3255c20d3d5cf0017f596d272fdffb6bd884c9d2`.
The independent reviewer is `test/parallel_block_profile/review_cluster_230.py`;
its results are in `parallel-blocks-stage230-review.json`. Raw extracted evidence
is retained under `.stage230-validation/review-55089/`.
The small review archive omits executables and raw checkpoint/VTK payloads;
their collected comparison reports and hashes were checked, not recomputed from
omitted payloads. Retain full results on bbserv until the work is committed.

Recommend committing this qualified cleanup on `codex/performance-investigation`
with its evidence. Main remains unchanged. The next implementation unit is the
complete inverse → shared primitives → tendencies → RK → forward/inverse cycle
on block-owned state, using existing intermediate oracles and zero-Domain-rank
tests. Establish which producers/consumers still require Domain `sol`/`q1`, move
those consumers together, reject fallback access, and delete the corresponding
storage only once the complete cycle passes. Do not return to local-halo cache
variants or another broad measurement campaign. This review itself makes no
source changes, commits or cluster submissions.
