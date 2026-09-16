# Handoff after Stage205 — 2026-09-16

## Objective and next action

Identify the concrete causes of the remaining legacy/block performance gap before
implementing another optimization. The user wants evidence for a coherent change
capable of saving at least **0.10 s per ordinary time step**, not a succession of
roughly 1% gains. Elapsed locations, inclusive timers, MPI waits and logical value
counts are not automatically removable work or predicted savings.

Start with a source-led comparison of the complete forward/inverse dynamics
wavelet paths against legacy, using the existing Stage204/205 measurements.
Match mathematical outputs, traversals, temporary preparation, indirect accesses,
local copies and necessary boundary dependencies. Investigate the numerical work
inside the RK tendency-capture envelope next. Do not schedule another broad
lifecycle campaign or implement an optimization without a supported mechanism.
The user has authorized the present commits; a new implementation campaign has
not yet been selected.

## Repository and preserved state

- Working repository: `/Users/kevlahan/.codex/worktrees/8778/wavetrisk_hydrostatic`.
- Branch: `codex/stage187-validated`.
- Accepted Stage202 source: `7c406f4d`; accepted comparison configuration uses
  `STAGE202=on FP_CONTRACT=off`. The Makefile still defaults STAGE202 to off:
  explicitly select the accepted configuration for experiments.
- Stage203 preservation commit: `b9847944`. Its numerical/build files match the
  frozen cluster-tested candidate. It remains experimental, `STAGE203=off` by
  default; candidate comparisons explicitly use `STAGE202=on STAGE203=on`.
- Stage204 preservation commit: `7f9effa5`.
- This handoff accompanies the Stage205 measurement/evidence commit. Resolve its
  exact hash with `git log`; start a new checkout from this branch's completed
  Stage205 commit, not from an older `parallel-blocks` branch tip.
- No push or merge was requested with these commits. Generated physics objects,
  caches and build-created `src/test*.f90` symlinks remain outside the commits.

## Read these first

1. `docs/parallel-blocks-cluster-205-54699.md` and its JSON companion: latest
   audited evidence, qualifications and decision.
2. `docs/parallel-blocks-cluster-204-54697.md` and JSON: broader operation ledger.
3. `docs/parallel-blocks-cluster-203-54695.md`: prototype results and limits.
4. `docs/parallel-blocks-stage205.md` and
   `test/parallel_block_profile/consumer-map-205.json`: observer scope and
   producer/consumer contracts. Contracts are hypotheses, not removal proofs.

## Latest results and conclusions

Stage205 job 54699, exclusive bb13, 83 ranks distributed 42/41 across two sockets,
J5/J7, contraction off, identical physics objects. All 19 runs and 41 gates pass.
The workload resumes checkpoint 3, completes 12 steps (11 ordinary), writes and
reloads checkpoint 4, and remaps after restart. Legacy uses a fixed tolerance
screen; block/oracle/observer comparisons are exact where required.

| Version | Bare ordinary-step median | Excess over legacy |
|---|---:|---:|
| Stage202 | 0.9275 s | 0.3845 s |
| Stage203 | 0.8795 s | 0.3365 s |
| Legacy | 0.5430 s | — |

Stage203's 0.048 s gain is repeatable but below the predeclared target. The current
ratios are about 1.71x/1.62x, not the historical 2.25x.

Both block calibration pairs pass. One legacy wavelet phase changes 2.0047%
against a 2% threshold; keep `observer_effects_require_review`. Do not waive the
threshold, discard stable structural evidence, or repeat the whole campaign just
to make the status green. Calibrate any new targeted observer against its own
adjacent controls.

Useful exclusions and limits:

- Scalar capture/restriction/temperature plans each rebuild once and reuse three
  times per ordinary rank-step. All selected plan work costs about 33 ms/step.
- One full scalar capture is followed by three compact captures. Geometry is not
  reseeded within a generation. The narrow scalar geometry preparation costs
  0.3 ms; this is not a measurement of every geometry calculation in the solver.
- Direct Domain/block imports plus publications total about 40 ms/step. Repeated
  stage/generation is not proof of unchanged field values; the RK marker also
  persists into later contexts.
- Checksum scans plus 13 validation Allreduces total about 19 ms/step. First
  checksum-reduction duration anticorrelates with scan duration, consistent with
  producer waiting. Do not add that waiting as independent recoverable savings.
- The 77 ms exclusive / 122 ms inclusive RK tendency-capture envelope contains
  actual scalar restriction and physical tendency computation. It is not simply
  copying. Full legacy tendency/oracle duplication is disabled in production.
- Signed instrumented rank-mean excess: dynamics wavelets 173 ms, other dynamics
  112 ms, adaptation/reconstruction 89 ms, other phases -3 ms; total 371 ms.
  This is a different statistic from the bare median gap of 384.5 ms.
- Capacity/swap pressure is absent. Cache/bandwidth inefficiency is not ruled out.
- No coherent, demonstrated avoidable 100 ms maintenance budget was established.

## Source entry points

- `src/time_integr.f90`: `prepare_block_multistage_tendency`.
- `src/multi_level.f90`: `block_tendency_compatibility_ml`, including oracle guards.
- `src/parallel_block_mpi.f90`: wavelet activation/inverse routines,
  `capture_block_domain_multistage_candidate_tendency`,
  `evaluate_candidate_block_velocity_recomposition`,
  `evaluate_candidate_block_exner_shadow`, field writeback/import routines,
  migration and plan-generation predicates.
- `src/parallel_block_inverse.f90`: native inverse workspace and optional
  Stage203 execution. Keep Stage202 and Stage203 paths distinct in comparisons.

## Artifacts and reproduction

Raw review archives are on the Desktop:
`/Users/kevlahan/Desktop/stage203-review-54695.tar.gz`,
`/Users/kevlahan/Desktop/stage204-review-54697.tar.gz`, and
`/Users/kevlahan/Desktop/stage205-review-54699.tar.gz`.
These contain logs/comparator reports, not executable or checkpoint payloads.
The Stage205 report records the archive hash and all 6,528 inventoried files
were checked. Raw analyses, lifecycle counts, clock conservation and summary
calculations were independently reconstructed.

The corresponding frozen kits and archives are currently under `/private/tmp/`
as `wavetrisk-stage20N-cluster-kit[.tar.gz]`; their hashes are in the committed
`docs/parallel-blocks-stage20N-kit.json` files. Do not assume temporary artifacts
will persist indefinitely or delete them. Preparation scripts 204/205 depend on
the exact Stage203 kit, which itself uses the Stage202 kit. All are based on the
pinned Stage196 bundle `9bc3a19556dab79b3cba3944cf5df53517389da3`.
Do not silently regenerate a different kit and call it the tested artifact.
One trailing-space-only cleanup in `analyze_batching_203.py` was made when
committing; the historical frozen kit is unchanged.

Before these commits, 19 focused tests passed, including actual Fortran field,
inverse-helper and observer tests plus real patch/hash/submit protocols using
a compiler stand-in. Earlier local restart/oracle runs and cluster evidence are
preserved in the stage reports; full solver builds were not rerun for committing.

Keep the user's single review-tarball collection workflow for future campaigns:
use the stage-specific `collect_cluster_logs_N.py`, then transfer that one archive
from bbserv to `kevlahan@murier.math.mcmaster.ca:~/Desktop/`.
