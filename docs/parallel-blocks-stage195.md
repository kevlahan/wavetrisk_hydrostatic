# Stage 195: explain the complete legacy/block performance difference

Status: cluster diagnostic campaign 54657 passed and has been independently
audited from the supplied logs. See [the cluster review](parallel-blocks-cluster-195-54657.md)
and its analysis/campaign JSON records. The matched phase ledger closes the
measured excess; complete causal attribution remains unfinished. Production
numerical algorithms are unchanged. This stage
supersedes the earlier proposal to select scalar replay or a particular ownership
redesign immediately. Its deliverable is a causal explanation of excess time,
with explicit uncertainty, before choosing another performance implementation.

## What the existing campaign establishes

Use Stage 194 job 54653, whose provenance, 28 completed cases and 64 numerical
and workload gates were audited separately. The repeated uninstrumented ordinary
step medians are 1.260 s for blocks and 0.558 s for legacy: a 0.702 s excess,
2.258x ratio. The following accounting is a different statistic: paired rank
wall-clock profiles, not printed CPU-time medians.

The new `compare_complete_gap_195.py` reparses all 83 ranks, requires matching
step/restart/remap coverage and printed numerical workload, and checks that all
exclusive regions sum to each rank's step. It forms signed differences rather
than adding block-only costs. It separately uses rank means and each
implementation's actual slowest rank on each step, never independent region
maxima. Remap and non-remap ordinary steps are reported separately.

Ordinary profile rank means are 1.258938 s blocks and 0.606897 s legacy, an
excess of 0.652041 s. Exclusive group differences sum to this excess with a
mean rounding residual of approximately 1e-8 s. This closes *time accounting*,
not *causal attribution*. Shared timer names do not guarantee identical work.

| Existing profile group | Blocks minus legacy, s/ordinary step |
|---|---:|
| Domain boundary handling | -0.159229 |
| Block reconstruction and migration | +0.102189 |
| Native inverse | +0.165299 |
| Scalar input and restriction | +0.184251 |
| Block ghost handling | +0.078704 |
| Native mass and velocity | +0.071268 |
| Shared Domain producers | +0.036847 |
| Domain/block publication | +0.046988 |
| Wavelet and compression | +0.086125 |
| RK and tendency assembly | +0.026456 |
| Other timestep work | +0.013144 |
| Total, before displayed rounding | +0.652041 |

These groups are a descriptive partition, not a list of independently removable
costs. For example, legacy boundary work can appear under a separately named
transport or ghost phase in blocks. Its negative boundary difference cannot be
ignored when considering the positive block groups.

The existing matched complete WT_after_step envelope, over the full profile
window including checkpoint, takes 4.653590 s in blocks and 1.283609 s in legacy.
Both execute 4,980 rank-summed calls. The excess there is 3.369981 s. All recorded
step wall times sum to rank means of 20.590256 and 12.610064 s, respectively:
7.980192 s excess. The remaining 4.610211 s lies outside the matched wavelet
transaction. These whole-window values include restart and must not be used as
ordinary-step percentages or reconciled directly to the 0.702 s timing median.

## A measurement issue to fix first

The profiled run's printed ordinary-step median is 1.310 s for blocks versus
1.260 s uninstrumented, a 3.97% increase. Legacy is 0.627 versus 0.558 s, a
12.37% increase. There is only one detailed profile per implementation, so these
increases include run variation and cannot serve as overhead corrections.
Nevertheless they make it unsafe to treat the detailed profile ratio as the
uninstrumented performance ratio. Profiling distorts the two paths differently.

Do not interpret the difference between 0.652041 s and 0.702 s as uninstrumented
missing work: they differ in instrumentation, aggregation, clock and run.

## Source-boundary audit

The production source confirms additional lifecycle work, but does not determine
its recoverable cost. In `main.f90`, the RK path materializes block state before
Domain physics/adaptation, then rebuilds the catalogue, source blocks and migrated
local state after adaptation. In `multi_level.f90`, the compatibility tendency
uses shared Domain producers and native mass/velocity, while scalar restriction
runs on final block owners. Thus representation and ownership transitions are
real. A full duplicate legacy timestep is not established: oracle branches and
several legacy numerical producers are explicitly disabled in production.

Existing coarse IDs also need boundary alignment. Legacy's dynamics timer wraps
`dt_step`, while the block dynamics envelope includes preparation/publication
around the multistage transaction. Legacy's adaptation timer is inside `adapt`;
the block adaptation envelope includes inverse reconstruction and block lifecycle
work. The same ID is therefore insufficient for a direct causal subtraction.

## Required next cluster measurements

Instrument both versions at the same physical workflow boundaries, using one
small common measurement implementation and the same monotonic wall clock.
The mutually exclusive outer phases should be:

1. Advance from the entering solution to the state ready for split physics,
   including RK work and required state publication.
2. Split physics, including representation refresh required by that consumer.
3. Adaptation and reconstruction until the new grid and solution are usable,
   including block catalogue/build/migration where required.
4. Remap decision and remapping, with event flags.
5. Output/checkpoint, including restart, separated from ordinary steps.
6. Explicit outside-phase residual, with total step conservation.

Nested matched transactions such as WT_after_step and a complete tendency
producer/consumer call belong in a second, non-additive report. Preserve phase
call counts, step identity, grid workload and state checks. Individual legacy
and block kernels may legitimately differ in invocation structure; compare the
complete physical transaction before treating component differences as excess.

Use three modes in each implementation: measurement disabled, coarse envelopes
only, and detailed profiling. Interleave repeated block/legacy pairs within each
mode on the same 83-rank nk allocation, retain fixed bindings, and buffer records
until after timed work. Do not add synchronization barriers to the normal timed
path. Calibrate the observer effect using the same wall-clock output in both
versions before using profiles for attribution.

For the dominant matched differences, collect call/work counts and classify
actual production operations as numerical arithmetic, allocation/initialization,
representation copy/packing, topology/address preparation, or communication
progress/wait. Include Domain and block work in both builds where applicable.
Track producer/consumer identity to establish whether a value is computed twice;
a timer called replay is not that evidence. Distinguish logical byte counters
from measured DRAM or network traffic.

Rank event chronology should compare phase durations and matched dependency
arrivals without assuming MPI clocks are synchronized across ranks. Whole-step
rank means close to maxima do not rule out phase imbalance masked by waiting.
Hardware-counter collection, if available on bb13, must be calibrated separately
and distinguish memory bandwidth, cache/TLB effects and instruction overhead.
Do not request a different MPI stack or infer a bandwidth limit from RSS alone.

Only after those comparisons identify a leading causal cost should a controlled
ablation change ownership, persistence, layout or redundant preparation. Preserve
full numerical work and the existing exact/oracle/legacy gates. Experiments that
skip work cannot enter the performance comparison.

## Acceptance and remaining uncertainty

A successful diagnostic stage reconciles the paired timestep excess within
measurement uncertainty, explains dominant differences using matched boundaries
and operation evidence, and explicitly reports unresolved causes. An additive
time partition alone is insufficient. There is currently no validated estimate
for how much of the 0.702 s can be removed by an ownership redesign, persistent
blocks, or a particular scalar optimization. No speedup target is promised.

Initial accounting is in `parallel-blocks-stage195-initial-accounting.json`.
Reproduce it with:

```sh
python3 test/parallel_block_profile/compare_complete_gap_195.py \
  /private/tmp/stage194-review-54653/results
python3 -m unittest discover -s test/parallel_block_profile \
  -p test_complete_gap_195.py
```

Six focused tests pass, covering signed conservation, actual slowest-rank
selection, missing ranks, invalid times, unmatched events and restart exclusion.
The tool has also completed against the supplied full Stage 194 archive.

## Implemented diagnostic campaign

`prepare_cluster_195.py` freezes observer-only patches on committed Stage 194
and the pinned legacy reference. The instrumentation transformer proves that
removing only its imports, calls and build-list addition restores every original
source byte. Neither the working production source nor its default selectors is
modified. Source manifests cover the climate driver and Makefile as well as
numerical source files. The common observer module is byte-identical in both
instrumented trees.

The five builds retain untouched block and legacy timing references, plus
instrumented block, checked instrumented block and instrumented legacy. All
observer modes on each side use the same diagnostic binary. Stage 193 scalar
columns and Stage 194 inner geometry are enabled for every block case, including
the untouched reference; Stage 192 sparse scatter remains off.

The common observer wraps the same climate-driver `time_step` call with a
monotonic wall clock. Minimal mode records step wall time and restart/remap flags
without reading clocks in phase or nested hooks. Coarse mode adds six mutually
exclusive physical workflow phases and nested transaction totals. Detailed mode
adds local event intervals and the established 67-region block / 65-region
legacy profiles. All common records are bounded, with explicit overflow errors,
and written after the simulation loop. No observer MPI calls or barriers are
added. Restart initialization does not reset common per-run records.

The campaign executes two checked modes, eight warmups and six interleaved
rounds of four block/legacy pairs: 58 runs. Mode order rotates; pair order
alternates. The full campaign has 146 numerical/workload gates. Every diagnostic
run must match its untouched implementation's checkpoint exactly. Checked runs
also match the optimized workload; cross-implementation comparisons retain the
existing fixed legacy limits. Rank affinity and memory/swap records are retained.

`analyze_workflow_195.py` rejects missing/duplicate rank records, inconsistent
schema/mode, incomplete phases, invalid times, mismatched event flags, phase
conservation errors, events outside their phase, and transaction count/time
disagreement. Paired phase differences must sum to the paired rank-mean excess.
The report includes per-round differences and observer overhead relative to
minimal mode, using the same wall-clock statistic. Bare printed/elapsed timing
remains a separate reference.

Transaction 1 is the matched complete WT_after_step envelope. Transaction 2 is
explicitly implementation-specific: block preparation includes capture and stage
preparation, whereas legacy measures trend_ml. Matching call counts do not make
those scopes numerically equivalent. Both are nested diagnostics and must not be
added to the mutually exclusive outer phases.

Hardware counters and globally synchronized MPI arrival tracing are not part of
this portable kit. Existing detailed self times/counters and the new local event
records support follow-up attribution; they do not establish DRAM bandwidth or
network service time by themselves. The campaign result is labelled
`measured_phase_accounting_requires_causal_review`, never a claim that every
performance cause has already been identified.

Local build/run evidence is recorded separately in
`parallel-blocks-stage195-local-results.json`; kit identities are in
`parallel-blocks-stage195-kit.json`. Follow `parallel-blocks-cluster-195.md` for
the build, run, and single-archive log collection protocol.

Local validation completed with fresh optimized block and legacy observers and
a fresh checked block observer. Five four-rank J4/J6 runs cover block detailed,
block minimal, block checked/coarse, legacy detailed and legacy coarse modes.
Every run completed 12 steps, checkpoint 2 -> 3/reload and two remaps, matched
its untouched reference checkpoint and printed workload exactly, and passed the
fixed cross-implementation numerical screen. Block checkpoints also matched the
saved independent oracle. Common phase/event coverage and conservation passed
for every rank. Local checked execution uses phase oracles off; full fresh
phase-oracle coverage remains in the cluster checked-on case.

The full helper suite passed 149 tests with five optional MPI tests skipped.
After final coverage checks, all 16 focused Stage 195 tests passed, including
actual observer execution at O0/O2, source/build guards, campaign failures and
paired accounting. Both frozen kit patches reproduce their source manifests,
the observer is identical on both sides, and standalone kit imports and shell
syntax pass. Local source identity also matches, allowing already-materialized
physics symlinks with identical target bytes in the block dependency tree.
Local timing values are not accepted performance evidence.
