# Stage 200: discriminate causes of the legacy/block performance gap

This is a diagnostic stage. The Stage 197 block candidate and pinned legacy
numerical sources, common checkpoint I/O fix, RK4, J5/J7, contraction-off flags
and physics are frozen. There is no new numerical optimization. Checked solver flags do not extend to the common external
physics objects, which are held identical across builds. Stage 196 is
only the reconstruction base for the supplied candidate patch, not a third
performance contender.

The objective is evidence supporting removal of a substantial part of the
approximately 0.69 s/step deficit at 83 tasks, not another one-percent change.
Use **0.10 s/step at 83 tasks** as a proposal-screening threshold (about 15% of
the deficit). This is neither a promised saving nor a statistical threshold.
A proposal must connect a measured excess to a source-level mechanism, show
repeatable supporting evidence, and state what is necessary work versus
avoidable work. MPI waits and producer work must not be double-counted.

## Experiment

| Configuration | Tasks | Placement |
|---|---:|---|
| n21-s1 | 21 | One socket |
| n42-s1 | 42 | Same socket, nested CPU subset |
| n83-s1 | 83 | Same socket, nested CPU subset |
| n83-s2 | 83 | Alternating ranks across the two sockets (42/41) |

Discover CPU/core/socket/NUMA/last-level-cache topology on the compute node;
never assume Linux CPU numbering equals Slurm core numbering. Within each
socket, spread cores round-robin across discovered last-level caches. Use
nested prefixes for the one-socket sweep. Pin each rank to one physical CPU
before MPI initialization. Both implementations use exactly the same rank map.
Use `numactl --localalloc` for consistent first-touch policy. Disable MVAPICH's
independent affinity setting. Observer binaries verify affinity after MPI
initialization as well as recording launcher affinity for every run.

Allocate all of bb13 exclusively. Run steps sequentially. The step launcher
omits `-c`, `--exact` and step-level `--exclusive`, so a smaller task count does
not restrict the step to an unintended CPU subset. The wrapper verifies that
every prescribed CPU is allowed before binding. This follows the Slurm srun
resource semantics: https://slurm.schedmd.com/srun.html#OPT_exact . Node/job
exclusivity remains required. NUMA/cache behavior still changes with placement;
a two-socket change alone does not prove a network or memory-bandwidth cause.

Keep the same checkpoint 3, `time_end=0.315`, `dt_write=0.01`, adaptive workload,
checkpoint 4 reload and post-restart remap in every run. Do not shorten the
physical integration interval at larger task counts. This is strong scaling,
not equal-work-per-rank scaling.

## Minimum campaign and calibration

Five fresh builds: candidate, checked candidate, legacy, and a minimal observer
for each implementation. No old detailed-timer or storage-census build.

42 short solver runs:
- Two checked runs, oracle off/on, at 83 tasks across two sockets.
- For each of four configurations: one warmup per implementation; two bare
  timing runs per implementation; one coarse-only observer run per implementation;
  one MPI/work observer run per implementation.

Bare timing pairs reverse implementation order on the second round; configuration
order reverses too. Coarse/probe order reverses implementation order. Repeats are
for detecting gross inconsistency, not claiming precise confidence intervals.
Additional repeats should be confined to an ambiguous configuration.

The same observer binary has two modes. Mode 0 records common whole-step,
phase and transaction clocks only. Mode 1 additionally times the existing
selected MPI scopes and accumulates integer work/volume counters. It records
no per-element clocks, per-patch/kernel timers, dependency tape or ordered
batch events, and inserts no barriers or I/O during timed steps. Source hook
calls are retained for a reversible observer-only patch, but unselected hooks
return without clocks. The disabled and active observer are calibrated against
the bare binary and against each other's parent envelopes. A 2%/1ms parent
change is a review flag, not an error bar or an overhead correction.

Record per-rank, per-step phase durations; MPI subset and remaining elapsed
time; MPI calls/volumes; available work counters; CPU affinity; solver RSS and
faults; and before/after node VM, CPU/load and available frequency-policy snapshots.
These snapshots surround the subprocess, not individual timed steps. Retain the actual slowest full-step
rank rather than summing independently chosen phase maxima. Counter 35
(temperature direct-flux nodes) is one useful work proxy, not all computational
work; all 42 counters and raw per-rank records remain available. No hardware
cache-miss or bandwidth measurement is claimed.

Same-configuration repeats and observers must match checkpoints exactly.
Checked/optimized comparison is exact at 83 two-socket tasks. Same-rank
one-/two-socket candidate warmups must match exactly. Across different task
counts, keep the established fixed numerical screen and identical timestep,
dt, Jmax and degree-of-freedom sequence; rank-dependent reduction rounding must
not be silently treated as exact equivalence. Cross-legacy comparisons use the
same existing screen, never a newly relaxed tolerance. Any failed gate rejects
the timing result.

## Decision table

| Observation | Supported next investigation | What it does not prove |
|---|---|---|
| Large excess remains at 21 tasks and outside-MPI cost per measured work stays high | Whole representation/capture/packing paths or repeated work | Which individual getter is expensive |
| Gap grows with task count alongside frequent MPI scopes and less work per rank | Ownership granularity, synchronization frequency and fixed per-rank cost | That message volume is the main cause |
| Large absolute and relative rank spread tracks measured work | A concrete ownership/load metric mismatch | That all MPI waiting can be removed |
| 83-task result changes markedly with socket placement | NUMA/cache/fabric sensitivity; follow with targeted hardware evidence | Pure network latency or bandwidth attribution |
| Detailed phase attribution changes with observer mode | Reduce measurement effects or repeat that case | A universal correction factor |

These are discriminators, not automatic optimization prescriptions. Assess both
implementations' scaling, absolute gaps and work-normalized cost. Fewer ranks
change cache footprint, active cores/clock behavior and decomposition; report
those confounders. `optimization_ready` remains false until human/source review
supports a substantial recoverable path. No production optimization is bundled.
