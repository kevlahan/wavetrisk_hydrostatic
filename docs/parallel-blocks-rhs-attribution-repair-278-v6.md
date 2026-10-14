# Stage 278 v6: attribution-led repair preparation

Prepared 2026-10-13. **Stage 266 remains accepted.** Two separately applicable
patches are prepared and locally tested: an attribution overlay for unchanged v5,
and a small route-construction repair. The repair does not establish a substantial
performance gain or a solution to 256-task scaling. A standalone performance
campaign for this repair is **not justified**.

## Why the repair direction changed

Communication is a candidate bottleneck, not an established dominant cause.
V3→v4 removed 9–18% of messages and 29–40% of payload for only 2.5–8.4% observed
ordinary-step improvement; its combined changes prevent assigning that gain to
communication alone. V4→v5 changed no numerical traffic and improved ordinary
steps 22–28% through preparation/storage changes. Stage 261's asynchronous path
regressed 4.97%; Stage 262 attributed most of its regression to the local-first
split/prefill design, with further cost from the asynchronous path. Those older
sources differ from v5, so this is design history, not current cost attribution.

V5's existing RHS timer includes input delivery, arithmetic, validation checks,
workspace copies, packing and MPI. Its 256-rank validation also enables an oracle.
Neither that timer nor the message count identifies the cause of the bare
160→256 slowdown. The revised work therefore measures these distinctions before
selecting communication batching, placement, layout or broader ownership changes.

## Prepared code

**Attribution-only patch:** extends the existing bounded detail profiler; it
does not change v5 numerical plan construction, ownership, arithmetic or exchange
schedule. Seventeen additional regions cover geometry compilation/distribution,
route and input-plan construction, retention, resident inputs, RHS work, RK
adapters, complete exchanges and optional pack/post/wait/install/accounting splits.
The RHS work residual includes arithmetic, workspace initialization/copies and
checks; it is not a pure floating-point kernel timer. Route construction includes
its metadata MPI as well as local work.

The default is coarse attribution when the existing detail profiler is enabled.
`WAVETRISK_PROFILE_RHS_FINE=1` enables the expensive exchange splits. With detail
profiling disabled, the new hooks do no clock reads, allocation, MPI or output.
Fine mode is not a bare-performance measurement. There are no added hot-path
collectives or barriers; output uses the existing end-of-window report. The
existing 256-step buffer limit and overflow reporting remain; the expanded static
step buffer adds about 68 KiB per rank. No communicator or solver-buffer lifetime
is extended.

The reader verifies complete twelve-step/rank coverage, rejects overflow,
non-finite/inconsistent records and oracle timing runs, and checks that self
regions partition each step. It reports regions from the **same slowest whole-step
rank**, plus per-rank ordinary means and remap/non-remap groups. It does not sum
independent regional maxima. MPI completion includes possible late producers;
process CPU includes MPI spinning. These observations do not alone distinguish
transport delay from producer lateness.

**Separate numerical preparation repair:** seed canonical fallback entries
directly in the existing owner directory with the exact old `[1,1]` ordering.
Previously every phase created explicit identity definitions, walked them in
count/fill passes, copied their keys/orders/slots/home addresses and then used
them to initialize the same directory. The generic compiler retains explicit
definitions as an independent reference. Remote overrides, final-writer selection,
signed aliases, cycles and epoch checks retain their contracts.

The demand mask already includes owned canonical components, so omitting the
explicit fallback records does not remove readers. This removes at least seven
integer entries per canonical fallback definition from temporary definition,
order and slot arrays (28 bytes with the current four-byte integers), plus their
processing and home-address materialization. It adds no persistent representation.
Actual numerical sends, payloads, exchange rounds and phase order are unchanged.
The preparation signature includes one additional integer to reject disagreement
about fallback mode; no extra collective round is introduced.

Both patches apply independently, in either order, to exact sealed-v5 files.
The working candidate contains both; the separate diagnostic overlay contains
only attribution. Neither is promoted to accepted production.

## Local attribution and measured limits

Four-rank J4 optimized fine-profile means on each ordinary step's slowest rank,
before the numerical repair:

| Region | Seconds/step |
|---|---:|
| RHS work residual: arithmetic, workspace/copy/checks | 0.822 |
| Route construction, including metadata MPI | 0.421 |
| Exchange packing and installation | 0.294 |
| MPI completion wait inside RHS exchanges | 0.203 |
| Input-plan construction residual | 0.187 |
| Per-exchange traffic accounting | 0.0028 |

Exchange regions also include input/RK/metadata exchanges; they are not solely
the derived numerical stream. The profile rules out treating the tiny accounting
scan as a worthwhile campaign driver locally. It identifies route construction
as a material region and the repeated fallback records as removable work, but
does not establish a cluster bottleneck. Other substantial inherited costs remain
visible, including wavelet/inverse work, source extraction, boundary copies and
split physics. No additional migration has been stacked.

Observed route-construction time falls from 0.421
to 0.377 s, about 10.5%, in the single before/after
fine-profile pair. That difference is only about 1% of complete profiled step time.
Other regions also vary, so it is not a precise isolated causal saving.

Uninstrumented ABBA ordinary-step run medians:

| Execution order | Variant | Seconds |
|---|---|---:|
| 1 | sealed v5 | 4.193454 |
| 2 | repair plus attribution, profiling off | 4.156123 |
| 3 | repair plus attribution, profiling off | 4.209652 |
| 4 | sealed v5 | 4.366745 |
| Separate control | attribution-only, profiling off | 4.152945 |

The repair mean is 4.182888 s versus 4.280100 s for v5 (−2.27%),
but v5 controls span 4.13% and the single attribution-only off control is faster
than the repair mean by about 0.7%. This is **inconclusive**, not a defensible
substantial speedup. Timesteps are not independent repetitions. No confidence
interval or J5/256-rank extrapolation is claimed.

The initial fine-profile run was 8.99% slower than its same-binary off control. This
motivated the separate coarse/fine selector. Final repair off/coarse medians are
4.345157/4.240529 s; their reversed ordering is additional run variation, not
negative profiling overhead. Use explicit off/coarse/fine comparisons on bb13
before interpreting a small region difference. Checked and oracle timings are
science evidence only.

## Qualification completed

- Serial optimized and checked builds; at most four local MPI ranks and no
  overlapping builds/replays.
- Final optimized ordinary, coarse-profile and independent-oracle replays,
  plus a final checked fine-profile replay: exact accepted-family science,
  all twelve states, restart/remap and VTK gates pass.
- All integer ownership, primitive/RK work, traffic and persistent-buffer ledger
  fields match v5 row by row in the oracle replay. The memory ledger excludes
  compiler temporaries and total RSS, so no RSS improvement is established.
- The separately built attribution-only overlay on unchanged v5 passes its own
  complete optimized coarse-profile science replay.
- Checked/optimized generic, canonical, shared-context and variable-width route
  tests pass, including signed aliases/cycles, remote writer overrides, ranks
  without canonical owners, bad epochs and communicator lifecycle checks.
- Fine/coarse/disabled profiling, invalid options and bounded-buffer overflow
  tests pass; six reader regression tests cover same-rank attribution and rejection
  of incomplete, invalid, overflowed or oracle records.

Raw failed development compile evidence is preserved alongside corrected tests.
The sealed v5 kit and all previous campaign evidence remain unchanged. No source
was transferred to bbserv, no new Slurm job was submitted, and no fresh GCC 13/
MVAPICH qualification or high-P performance result is claimed for this work.

## Decision and next evidence needed

Retain the canonical fallback change as a tested proposal, but do not run a
standalone performance campaign for its small mixed result. The next substantial
repair requires **cluster attribution**, not additional guidance about the goals.
The proposed 10–15% improvement remains an engineering target without a supported
forecast. No final scaling-efficiency threshold has been invented.

Use the attribution-only overlay on unchanged v5 for a bounded matched bb13
diagnosis at 83/160/256: ordinary off/coarse pairs first, with a separate fine
exchange run only where its result changes the repair decision. Preserve exact
science, fixture and placement identities and separate remap/non-remap results.
This is a measurement specification, not a sealed or submitted campaign kit.

The decision is whether the 160→256 loss is in preparation, RHS work, copy/layout
cost, inherited inverse/wavelet/compatibility work, rank imbalance or exchange
completion. A large wait region alone does not authorize a claim of network
dominance or guarantee that batching will help. If producer lateness and transport
remain indistinguishable, that uncertainty must remain explicit.

Choose a repair with a measured budget capable of recovering the remaining
28–48 ms versus 266 and improving the 160→256 slope. Legacy still requires 241–260 ms
at 83/160; most of that gap predates 278. Preserve connected RHS/RK ownership and
paired retirement while removing demonstrably redundant production work. Keep
ordinary-step absolute time and scaling as separate gates, retain memory limits,
and qualify subdivision/placement separately before the >=640-rank goal. Whole-
timestep independence and shared-memory transport are not established shortcuts.

## Evidence and reviewable patches

- [Attribution-only patch](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/repair-v6/attribution-only.patch),
  [separate canonical fallback repair](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/repair-v6/canonical-fallback-repair.patch),
  [overlay identities](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/repair-v6/patches.json).
- [Final qualification](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/repair-v6/final-qualification.json),
  [attribution-only qualification](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/repair-v6/probe-qualification.json),
  [traffic equivalence](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/repair-v6/traffic-equivalence.json), [local metrics](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/repair-v6/metrics.json).
- [Decision](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/repair-v6/decision.json), [bounded diagnostic specification](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/.stage278-validation/repair-v6/diagnostic-plan.json),
  [attribution reader](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/test/parallel_block_profile/connected_rhs_278/attribution.py).
- [Job 55502 review](/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic/docs/parallel-blocks-stage278-job55502-results.md),
  [Stage 261](/Users/kevlahan/.codex/worktrees/26d8/wavetrisk_hydrostatic/.stage261-validation/review-55272.md),
  [Stage 262](/Users/kevlahan/.codex/worktrees/26d8/wavetrisk_hydrostatic/.stage262-validation/review-55281.md),
  [authoritative scorecard](/Users/kevlahan/.codex/worktrees/26d8/wavetrisk_hydrostatic/docs/parallel-blocks-scorecard.md).
