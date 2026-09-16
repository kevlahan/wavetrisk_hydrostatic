# Performance-gap diagnosis: current objective and evidence requirements

This document records the user's instruction after review of Stage 197 job 54663.
It supersedes optimization recommendations in the Stage 197 proposal and review.

## Objective and scope

Find and fully document the main causes of the legacy/block performance gap
before starting another optimization. An optimization must have evidence that
it can recover a substantial part of the roughly 0.69 s ordinary-step excess.
A slow routine, a large byte count, a redundant-looking expression, or a list
of named timers is insufficient evidence by itself.

Work is diagnosis only: source/caller audits, calibrated observation, analysis
and documentation. Keep the existing numerical executables and arithmetic
unchanged. Preserve Stage 197 as a completed experiment; do not enable another
production selector, rewrite a kernel, alter ownership scheduling, change
arithmetic flags, or commit an optimization as part of this diagnosis.
Observational builds may add measurements in disposable copies, with exact
reference/workload checks and explicit observer-effect checks.

The comparison pair must be identified in every report. The validated Stage 197
candidate and pinned legacy in job 54663 give 1.250 versus 0.5565 s, or 0.6935 s
bare median excess. Retain the untouched Stage 196 executable as the historical
control; it gives 1.265 versus 0.5565 s in the same campaign. Do not silently
switch between those pairs when attributing savings or quoting the remaining gap.
The standard workload remains J5/J7, RK4, 83 ranks on exclusive nk/bb13, contraction
off, identical physics objects, fixed fixture and core placement.

## Existing ledger: measured locations, not established causes

The most thoroughly calibrated detailed location ledger is Stage 196 job 54659:
means over the same 11 ordinary steps, 83 ranks and six repeats. Its total is
0.659377 s, not the bare median excess. Preserve its statistic and scope rather
than rescaling it to force agreement with 0.69 s.

| Disjoint location | Measured mean excess, s | What is still unknown |
|---|---:|---|
| Dynamics WT: boundary/ownership routines | 0.205488 | Required versus redundant routing; local preparation versus dependency delay; recoverable critical-path cost |
| Dynamics WT: other kernels/staging/residual | 0.075154 | Matched arithmetic, reconstruction, staging and representation costs |
| Dynamics outside WT: boundary/ownership routines | 0.117894 | Net specialized-path cost, including savings in the shared Domain path; reasons for mass/scalar waits |
| Dynamics outside WT: other kernels/storage/residual | 0.149248 | Matched fused producer work versus block storage/addressing/checking and numerical work |
| Adaptation/reconstruction | 0.114188 | Mask/topology work versus compact-state construction, migration, installation and repeated setup |
| Other outer phases, signed net | -0.002595 | Small signed contribution; retain in conservation |
| Total | 0.659377 | Location accounting is complete; causal accounting is incomplete |

The first four rows total about 83% of this measured mean excess. This percentage
never established that a particular optimization would remove 83%, nor that the
costs were predominantly redundant arithmetic or network traffic.

Stage 197 locates its own small effect: about 16.8 ms less rank-mean WT time,
consistent with the 15 ms observed bare-step reduction. It leaves inverse
dependency counts, requests and values unchanged. It does not explain most of
the WT rows above or the outside-WT/adaptation rows.

Stage 197's new legacy profile adds 8.7% to its bare printed-step median, compared
with 1.6% for the candidate. Its raw profile excess of about 0.622 s cannot be
presented as an exact decomposition of the 0.6935 s bare median gap. Stage 195/196
and197 statistics and observer effects must remain visible in the report.

## Evidence already established and its limits

- All Stage 197 correctness/workload gates pass; numerical failures do not explain
  this observed slowdown. This is the existing fixed legacy numerical screen,
  not bitwise block/legacy equivalence.
- Rank placement is consistent and node swap deltas are zero. The larger block
  memory footprint is measured; a bandwidth/cache/NUMA bottleneck is not yet
  demonstrated by that fact.
- Both implementations have four dynamics WT calls per rank per step, plus the
  adaptation-phase WT. The excess is not explained by extra calls to the whole
  dynamics WT driver. Its internal work and storage lifetimes differ.
- The complete legacy tendency is invoked through the oracle branch only in
  the inspected production block RK path. Blocks retain shared Domain primitive/
  physics producers and use native mass/velocity/temperature paths. This rules
  out claiming a wholesale duplicate legacy timestep from routine names; it
  does not establish that every retained producer or intermediate is necessary.
- The shared generic Domain boundary path has lower measured traffic and traversal
  counts in blocks. Additional specialized paths must be accounted alongside
  those savings, rather than counting only positive block costs.
- Selected temperature-flux and divergence evaluation counts are lower in blocks.
  The scopes are limited and legacy producers are fused, so these are not a full
  equivalence census. They do not establish the cost per matched computation.
- A production scalar call requests a 50-slot compatibility record before using
  one value in its mass shortcut or six in its temperature flux. Its measured
  source-level sample count implies about 8.29 million requested slots per rank
  per ordinary step. Generated code and actual exclusive time have not established
  the instruction, cache or critical-path cost. The containing direct-flux scope
  is about 55 ms rank-mean time; this call alone cannot stand in for the whole gap.

## Diagnostic work packages

These are investigations, not an optimization implementation sequence. Keep all
three major regions in scope instead of repeatedly selecting one local routine.

### 1. Repair the comparison and measurement foundation

Use the same completed step set, masks, levels, field families and physical
outputs on both implementations. Define matched work at the fused-operation
level where legacy cannot be separated without changing its execution. Keep
shared work in its own row; do not arbitrarily allocate it to temperature.

Batch counter reporting outside leaf loops. Compare observer-disabled/minimal
and focused observations against untouched executables on the same allocation.
Check parent-region distortion as well as whole-step distortion. If distortion
is comparable to the cost being attributed, that measurement is insufficient;
reduce hooks, use sampling, or isolate unchanged kernels for diagnostic replay.
Do not subtract one global 8.7% correction from every legacy region.

Retain numerical expressions and control flow when inserting observation; prove
source changes are observational where possible and require exact reference
checkpoints and matching workloads. Numerical tolerances are unchanged.

### 2. Trace complete dynamics WT transactions

Map transform, compression, scalar/vector reconstruction and final publication
against the corresponding legacy transaction. For every intermediate, record:
producer; authoritative owner; required levels/fields; consumers; invalidation;
local copies; remote transfers; repacking; and repeated initialization.

Source entry points: `adapt_mod::WT_after_step`, the native provisional/final WT
calls in `time_integr_mod`, the native wavelet/compression/inverse routines in
`parallel_block_mpi_mod`, and gather/scatter/alias operations in
`parallel_block_inverse_mod`. Compare with the pinned legacy adapt/inverse path.

Separate exclusive numerical/traversal work from representation preparation,
installation and MPI completion intervals. Attribute the actual 86–90 ms inverse
completion/posting region and 55 ms preparation region; do not assume the waiting
portion is independently removable or caused by message bytes. Gather/scatter/
alias record counts currently identify dependencies, not the cause of delay.

### 3. Trace the complete tendency/RK pipeline

Compare legacy `trend_ml` and its shared/fused producers with block
`prepare_block_multistage_tendency`, `block_tendency_compatibility_ml`, native
mass/velocity, scalar capture/restriction/divergence, and RK assembly/publication.
Inventory which equations are computed, which intermediates are copied or
reconstructed, and which retained Domain outputs have actual consumers.

Count work by physical field, level, active/masked/boundary case and RK stage.
Measure record reconstruction, address resolution, repeated geometry access,
checks, traversal, numerical operations and output installation at batch scope.
Audit all production users of the scalar storage interface; the 50-slot call is
one hypothesis in this inventory, not a selected optimization. Inspect generated
code before treating requested logical reads as executed loads or function calls.

Explain native mass/scalar waits through producer work and dependency sequences.
Use per-rank distributions, matching dependency keys and validated timing scope.
Unsynchronized rank-local timestamps cannot establish global arrival times;
leave that causal split unresolved when the evidence cannot distinguish it.

### 4. Account for adaptation/reconstruction separately

The roughly 0.114 s excess must not disappear because the dynamics rows were 83%.
Separate threshold/mask updates, Domain adaptation/inverse reconstruction,
block lifecycle completion, catalog construction, source-block construction,
migration/installation, plan retirement/recreation and retained reconstruction.
Compare equivalent legacy work and count actual topology changes and rebuilds.

The block `main_mod::time_step` adaptation section explicitly includes
`build_parallel_block_catalog`, `build_source_blocks`, `migrate_blocks` and
`retain_post_grid_change_block_reconstruction` when a block rebuild is pending.
This identifies work to measure; it does not prove a rebuild can be skipped
or that rebuilding accounts for the entire adaptation excess.

## Required deliverable before an optimization proposal

Produce one reconciled ledger with these columns:

| Cause or matched task | Legacy work/time | Block work/time | Excess | Evidence | Required versus avoidable work | Recoverable step time and uncertainty |
|---|---|---|---|---|---|---|

Every major location row above must map to measured causal rows or an explicitly
sized unresolved remainder. Report workload counts and scope definitions beside
timings. Include negative differences, observer effects, and the difference
between additive rank means and elapsed critical-path/slow-rank measurements.
Never sum independent per-region rank maxima or add a parent and its children.
A zero unclassified timer does not imply a zero causal unknown.

For a proposed change, establish all of the following before implementation:

1. A demonstrated mechanism, not just a hot routine or suspected inefficiency.
2. A measured amount of unnecessary or avoidably expensive work in matched tasks.
3. Evidence that this work delays step completion, with dependency constraints.
4. A conservative estimate of a substantial recoverable portion of the remaining
   excess, including new costs and uncertainty. Avoid double-counting producer
   savings and the resulting reduction in other ranks' waiting.
5. A correctness-preserving producer/consumer design and fixed validation plan.

If the recoverable cost is small, cannot be separated from observer effects,
or remains an upper bound that assumes required work disappears, it does not
qualify. Continue diagnosis. No new optimization is authorized by a hypothesis
appearing in this document.

## References

- [Stage 195 location ledger](parallel-blocks-cluster-195-54657.md)
- [Stage 196 boundary/operation ledger](parallel-blocks-cluster-196-54659.md)
- [Stage 197 measured effect and limitations](parallel-blocks-cluster-197-54663.md)
