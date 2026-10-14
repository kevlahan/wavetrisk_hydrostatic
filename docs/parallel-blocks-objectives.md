# Block refactor objectives — reset 2026-10-09

The user authorized this reset after the project-wide review and Stage277 return.
This document supersedes earlier next-step recommendations where they conflict.
Historical plans/results remain evidence, not the current work queue. Maintain
the [scorecard](parallel-blocks-scorecard.md) after every reviewed campaign and
every accepted production change.

## Required outcome

Complete independent block ownership of the supported production timestep while
removing non-essential Domain numerical pathways and repeated representation
transfers. Achieve performance comparable to or better than legacy main at
P<=D and useful efficient execution at P>=4D, where D is the geometric Domain
count. J5 has D=160: 256 ranks is an intermediate check, not the final target;
at least 640 ranks is required for fourfold qualification.

Geometric Domain identifiers, required external I/O adapters and independent
validation oracles may remain. Substantial interior numerical work restricted to
Domain owners, duplicate production computation, or recurring whole-state
publication merely to feed legacy numerical readers is unfinished migration.
Renaming storage or extracting Domain-free arithmetic does not complete ownership.

The final deliverable is a qualified implementation, not a diagnostic framework,
a count of removed references, or a runnable-but-slower P>D configuration. No
numerical efficiency threshold has yet been agreed: report speedup and relative
efficiency honestly and agree the threshold before final scaling acceptance.

### Workload continuity — user clarification, 2026-10-14

Keep the current J5/D=160 workload as the deciding absolute-performance and
above-Domain scaling track. The user may later consider J6 coarsest/J8 finest
with four times the Domains and more numerical work. This is complementary:
with D=640,160/256 tasks are below the Domain count,640 only reaches it, and
P>=4D requires2,560 tasks. A larger-workload result cannot replace breaking the
Domain barrier on the original workload. No J6/J8 campaign is requested now.

## Immediate substantial deliverable: connected production RHS ownership

Start from accepted Stage266 `ff740744608731a05d53771c32a3055309f3a586`.
Reuse the verified direct numerical kernels, Stage274/275 dependency evidence,
Stage276 resident-input evidence and Stage277 assignment checks. Do not adopt
the Stage277 geometry-owner cache or Stage273 partition candidate implicitly.

The next implementation must connect resident block state to authoritative
mass/velocity RHS results and their real RK consumers on block owners. Its
boundary includes the shared flux/PV/KE/vorticity producers, necessary pressure/
Bernoulli/Exner work, restriction and halo phases, source/gradient, diffusion,
required physics dependencies and the affected RK state. The implementation may
use multiple direct kernels, but must preserve numerical fusion, column order,
final shared-vorticity write order and consumer epochs.

This is an architectural milestone, not a promise that all of it fits a tiny
two-file stage. Develop numerical components in reviewable local steps; qualify
each independently activatable connected production slice before building the
next on it. Do not ship an isolated block producer that must return all its
primitives to unchanged Domain consumers, or stack several unmeasured cutovers.

Before changing placement, maintain a concrete contract for each field:

| Required entry | Purpose |
|---|---|
| Producer, current owner and replacement owner | Identifies actual work being migrated |
| All production readers, including physics and boundaries | Prevents discovering a mandatory bulk return after implementation |
| Topology, mask, ownership and state/halo epochs | Defines valid cache and operand lifetimes |
| New communication, preparation, storage and publication | Makes replacement cost explicit |
| Removed computation and transfers | Proves the new path replaces work rather than adding another representation |
| Correctness and full-step acceptance evidence | Defines when this slice is finished |

Use the existing source and evidence to complete this implementation contract;
another broad timing audit is not a prerequisite. A new diagnostic is justified
only by a named unresolved dependency or decision. If a retained physics adapter
is necessary during a slice, account for its inputs/outputs and Domain-owner work,
and name its eventual removal. It cannot silently become the permanent design.

Acceptance requires actual production execution on ranks without Domains,
authoritative outputs consumed by the real solver, absence of the displaced
Domain computation for that slice, and retained scientific/oracle gates. Include
all preparation, communication, boundary, physics, restriction and publication
costs in bare complete-step timing. Reference packet streams and duplicate
arithmetic belong only in validation. Do not recreate the Stage251–255 per-value,
dependency-depth owner interpreter.

## Remaining architectural milestones

1. **Resident timestep and inverse ownership.** As readers move, retire tendency
   exports, RK compatibility materialization and provisional Domain interior
   publication for the affected fields. Migrate the remaining geometry-owner
   inverse dependency/alias chain while preserving the locality that made
   Stage229 faster. Required boundary/scaffold refreshes stay until replaced by
   an equivalent block contract. Adaptation, remap and lifecycle readers must
   be included before calling the whole timestep independently owned.
2. **Sufficient granularity and balanced work.** Extend subdivision beyond the
   current maximum of four blocks per Domain and validate memory/halo costs.
   At 4D ranks the present catalogue permits at most one block per rank, which
   cannot balance unequal indivisible work effectively. Evaluate partition
   changes separately: minimizing leaf weights is not a runtime guarantee.
3. **Real multi-node performance.** Qualify matched J5 160/320/640-rank execution
   on bb hardware, and higher counts when block capacity and allocation permit.
   Keep the same hardware/placement policy for comparisons; never splice nk and
   bb measurements into one scaling curve. Plan a bounded multi-node readiness
   check when it can expose implementation constraints, rather than postponing
   all inter-node testing until the end.
4. **Supported product completion.** Validate longer representative adaptive
   workloads and enumerate remaining model-mode restrictions. Current P>D support
   is restarted compressible, non-mode-split RK4 climate with Simple physics,
   without NCAR topography/SSO. Cold starts and other required modes need explicit
   qualification; they cannot be silently counted as complete from J5 evidence.

The primary numerical ownership work and granularity are both necessary. They
can be designed alongside one another, but production changes must be measured
separately. Explicit MPI shared memory, isolated micro-optimizations and cosmetic
Domain-reference cleanup are deferred unless a concrete critical-path dependency
requires them.

## Performance and correctness discipline

### Production timing priority — user clarification, 2026-10-13

Prioritize absolute performance and scaling of ordinary complete timesteps,
including their preparation, computation and communication. Production checkpoint
or data saves occur only every O(10^3) steps; remapping occurs approximately every
five steps. The short fixture deliberately exercises saving/reload often, so its
whole-trajectory elapsed time is a secondary lifecycle/overhead diagnostic and
must not outweigh ordinary-step performance as a production throughput proxy.

Keep the established ordinary-step exclusions for comparisons: omit checkpoint/
data-save steps and the first step after startup/reload; retain ordinary remaps.
Report remap and non-remap steps separately as supplementary evidence so a median
does not hide a recurring remap penalty. Assess recurring remap overhead at the
production frequency; only construct a frequency-weighted estimate when the
incremental event costs can be identified from measurements. Preserve all science,
reload/remap, ownership and memory gates and the existing dual goals.

- Compare each connected change to its exact accepted predecessor and read-only
  main on the same allocation at supported counts. Include above-Domain execution.
  Initially use the short 83/160/256 campaign; extend to multi-node milestones.
- Do not promote a repeatably slower change and defer recovery. An inconclusive
  small result is not an established gain. Repeat only when resolving uncertainty
  changes the acceptance decision. Non-regression is an intermediate gate; it
  does not satisfy the final requirement to close the existing legacy gap.
- Keep optimized/checked/oracle validation separate from uninstrumented timings.
  Preserve exact accepted-family comparisons, established main-relative limits,
  adaptive topology, masks, reload/remap, boundary epochs and zero-Domain ranks.
- Use ordinary-step exclusions consistently; do not treat timesteps as independent
  repetitions. Longer workloads supplement the short correctness fixture before
  final acceptance. Report memory and preparation costs as well as runtime.
- Retain short, proven delivery scripts with immutable sources and fixture hashes.
  Capture failed-job evidence and real Slurm accounting automatically. Avoid
  changing campaign plumbing merely because a stage number changed.

### Attribution-led repair discipline — user clarification, 2026-10-13

Choose the next substantial repair from bottleneck attribution. Message counts,
payload reductions and time inside MPI completion do not by themselves establish
transport as the dominant bottleneck. Use existing evidence first; add bounded
diagnostics only for a named unresolved repair decision. Distinguish preparation,
useful work and memory/copy costs, rank imbalance, retained adapters and exchange
costs, and quantify profiler perturbation with instrumentation-off controls.
Keep small mixed local optimizations from driving standalone cluster campaigns.
An absolute-time improvement that leaves 256 slower than 160 does not resolve the
scaling failure. Preserve both original goals and all existing science/memory gates.

## Measures of progress

Update the scorecard with accepted source, pending candidates, retired production
work/transfers, actual owners of remaining substantial work, matched main gaps,
P>D speedup/efficiency, next deliverable and its concrete blocker. Every campaign
must either qualify a production change or resolve a decision needed for one.

Stage277 is closed as numerically validated but unpromoted: small mixed timing
changes, unchanged ownership, no further standalone cache campaign planned.
Its reusable numerical evidence is a means to the ownership deliverable, not a
replacement for that deliverable.
