# Stage229 — recoverable work and the current legacy gap

Requalify the existing Stage203 combined implementation on accepted Stage226,
then locate its net effect and the remaining excess with a small common observer.
Production remains commit `1d7bc050c0f787ad3748c7a66fed6767fd6f640f` with Stage203
off. Main remains `b82647b36eb825aed0c0aae4e92adcf4ab46d608`. Nothing is promoted
automatically. Stage227's deletion and Stage228's publication-copy target are
excluded.

## Why this candidate

Stage203/204/205 bare comparisons reported 54.5/40.0/48.0 ms per-step gains.
Those older runs predate Stage215, use the historical timing selection, and have
not been reprocessed here; these figures justify requalification, not a current
speedup prediction. The original prototype missed the former 100 ms threshold.
It did not have zero benefit. The historical instruction evidence also showed
about 13% fewer wavelet instructions.

The candidate already exists in source behind `STAGE203=on`. It runs scalar and
inner-vector inverse kernels alongside the outer-vector kernel in one workspace,
omitting intermediate gather/scatter and compact ghost refreshes. Required entry,
exit, signed-alias and consumer-version work remains. Batched scalar field access
is enabled with the inverse change, as in the successful historical combined
trial. This is not a new owner migration or another local-halo cache variant.
No feature ablation is included; the campaign cannot assign the net benefit to
one of these two features independently.

The historical Stage207 signed location ledger was approximately 179 ms dynamics
wavelets, 81 ms tendency context, 24 ms other dynamics, 92 ms adaptation and
−3 ms other phases. Those are old locations, not today's removable budgets.
Stage215 already recovered substantial work. Stage228 put the proposed vector
publication copy below 0.5 ms; that target stays closed.

## Sources and builds

Ten fresh binaries on one nk node:

- J4: Stage226 baseline, Stage203 optimized, Stage203 checked.
- J5: baseline, Stage203 optimized, Stage203 checked, pinned legacy, and a
  separate coarse-observer binary for baseline, Stage203 and legacy.

Baseline and candidate numerical source are identical; compile-time selectors
are the intervention. All block builds retain Stage202 on, Stage197 resident
wavelets on and normal accepted runtime selectors. FP contraction is disabled
in every solver build, including legacy. All binaries use the same 23 external
physics objects. Compiler commands, hashes, source manifests and logs are sealed.

Legacy comes from a read-only `git archive` of the pinned main commit, with only
the established checkpoint I/O hardening and dependency include as necessary.
The actual main branch/worktree is never checked out or modified. Observers are
applied only to disposable source trees; stripping their added statements restores
all original source bytes.

## Correctness and performance funnel

1. J4 at **41 tasks / 40 Domains**: baseline, checked candidate, optimized
   candidate with all numerical oracles. This gates zero-Domain participation,
   one-Domain senders, restart and remap before performance work.
2. J5 at **16 tasks**: baseline, checked candidate, optimized candidate with
   oracles, optimized candidate without oracles. Then one fresh bare A/B screen.
   Advance for at least 1% gain; a 3% slowdown is a clear regression. Otherwise
   stop this candidate's performance path.
3. For an advancing candidate, early **83-task** baseline/checked validation and
   one fresh bare A/B pair. Require at least 1% gain and at least 20 ms to justify
   confirmation of the 25 ms target.
4. Fresh **A–B–B–A at 83 tasks**, independent of all screening and diagnostic runs.
   The target is **25 ms per ordinary step**, chosen for this already-implemented
   but architecturally nontrivial candidate. Both adjacent gains must support the
   target. Add only one further ABBA batch if paired gains disagree across zero
   or the threshold, or the aggregate gain is within 5 ms of the target. Continued
   ambiguity is inconclusive; no indefinite repetitions.
5. Run the bounded location comparison below. No timing from checked or oracle
   runs enters acceptance. No automatic promotion, additional variants or deeper
   instrumentation follows the result.

The local and cluster 41/40 gates cover above-Domain execution for this workload;
they do not newly certify every task count, physics setting or maximum-rank
configuration. Any eventual promotion must review the full retained capability
contract. This campaign prioritizes normal 16/83 performance over scaling.

## Location comparison, including negative results

Use fresh bare controls and identical coarse observers on **83 tasks**. A is
current blocks, B is Stage203, L is legacy. After a candidate reaches the 83-task
pair, run this fixed 12-run sequence:

`A B L AO BO LO LO BO AO L B A`

Each implementation has its own bare–observer–observer–bare calibration. These
runs do not enter candidate confirmation. If the 16-task screen stopped the
candidate, omit B/BO and run the eight-run current-gap sequence:

`A L AO LO LO AO L A`

Thus a negative candidate screen still produces a current A/L location ledger,
without sending the rejected candidate on to larger-task optimization tests.
A failure of solver correctness, source identity or placement stops the campaign
immediately and preserves partial logs; it does not trigger this fallback.

Calibration for each implementation is fixed before submission:

- Mean of observer run medians minus mean of bare run medians: absolute shift
  no greater than max(3 ms, 1% of that implementation's control).
- Each adjacent pair: absolute shift no greater than max(5 ms, 2%).
- Within-version run-median spread: no greater than 2% of control.

If any implementation fails, the location comparison is explicitly unqualified;
raw evidence remains available and no automatic observer reruns occur. A calibration
pass concerns whole-step resolution and is not proof that every individual region
is unbiased. Tiny phase differences are not optimization budgets.

The observer adds six outer phase boundaries and two batch transaction envelopes,
with fixed buffers and records written after the complete run. It calls MPI_Wtime
but introduces no MPI communication or synchronization, point/layer clocks,
PMU sampling or timed-step file I/O. Local clocks are never subtracted across ranks.

The disjoint report comprises setup, dynamics wavelets, tendency context,
remaining dynamics, physics, adaptation/reconstruction, remap and ordinary output.
It verifies four WT and four tendency transactions per ordinary step, no overlap
between these transactions, rank/step completeness, phase conservation and event
flags. Tendency contexts are implementation-specific complete preparation calls,
not claimed identical arithmetic kernels. Wavelet/other-dynamics locations retain
all their required publication, transfer and preparation costs.

For each pair report **signed** A−L, B−L and A−B location differences, keeping
negative offsets. Compute means across the same ranks and ordinary steps so they
sum to the observed whole-step mean difference. Do not sum region medians or
independently selected rank maxima. Retain explicit descriptive bridges between
bare run-median differences, bare ordinary-step means and observer means.

Only an uninstrumented intervention establishes a measured recoverable portion.
The A−B location changes show where that net effect appears; B−L indicates where
excess remains. A large positive remaining location does not prove redundancy.
Any later candidate must identify legacy's avoided work and the full replacement,
including dependencies and costs outside its immediate region. MPI wait is not
independent savings on top of the work causing that wait.

## Workload, timing and placement

Use the frozen adaptive 12-step restart fixtures. Exclude one-based steps **1,8,9**
(startup/data saving, checkpoint/reload, first full post-reload step). Retain steps
2–7 and 10–12, including ordinary remapping. Validate actual logged events; collect
high-precision `<run_id>_log`. Do not change `dt_write` or relabel this as a continuous
no-output trajectory. Timesteps are correlated observations, not independent runs.

J4/J5 runs retain identical printed adaptive trajectories. Same-count block A/B,
checked and observer results must match exactly. Legacy observers must match their
own bare reference exactly; legacy versus blocks uses the unchanged fixed numerical
screen, including exact topology/threshold/soil requirements. All runs complete
both checkpoint loads and both remaps.

nk has two 128-core sockets. Slurm allocates the complete node, one physical core
per rank, with the established 16-rank one-socket map and 41/83-rank maps alternating
sockets. Build and execution use the same native-compilation host. Exclude bb04;
no `--exact` or step-level `--cpus-per-task` restriction. The run records/proves every
rank placement. bb scaling remains deferred.

Normal successful path: **29 runs**, or **33** only if one extra confirmation batch
is required. A failed 16-task performance screen uses 17 runs including the current
A/L diagnostic; a failed 83-task screen uses 25 including A/B/L diagnosis. These
counts include correctness runs. Performance failure never causes another candidate
variant or an unbounded measurement campaign.
