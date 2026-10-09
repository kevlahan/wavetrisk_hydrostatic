# Stage277: production integration of the fused primitive kernel

Prepared from accepted Stage266 `ff740744608731a05d53771c32a3055309f3a586`.
Main control is `6a245db7f1632b0205c2bebe95e5f0cd535c2582` with the established
checkpoint I/O patch. The separate Stage273 partition candidate is excluded.

## Scope correction following the implementation trace

The initial recommendation after Stage276 was a complete block-owned producer,
restriction and consumer replacement. Stage276's nine initial interior outputs
do not yet supply that replacement: retained physics and Bernoulli/Exner paths
still need primitive fields on geometric Domain owners, and shared vorticity
has four ordered writes rather than just the two initial centre values tested
by the prototype. Returning bulk block-produced primitives would add an
intermediate route of the kind the new architecture program is meant to avoid.

Stage277 therefore makes a smaller production integration first. It **does not
migrate execution ownership**, retire the native velocity export, or qualify
efficient execution above the Domain count. That larger requested migration is
still outstanding. This limitation must accompany the kit and its results;
do not present this candidate as a completed producer/consumer ownership slice.

The concrete deliverable is a fused array kernel used by the real production
consumers, a reusable geometry/address cache, complete assignment-point checks,
and a short full-step timing campaign. No new MPI route or collective is added.
Explicit MPI shared memory remains separate.

## Source change

The isolated two-file overlay modifies `ops.f90` and `multi_level.f90`:

- Cache 37 horizontal geometry values and 16 operand addresses per positive
  interior node, once per existing writeback topology/ownership generation.
  Cache content is independent of masks and vertical layer. Current mask
  predicates continue to control the retained physics and restriction paths.
- Compute one patch batch using direct array arithmetic: mass-advection flux,
  edge PV, kinetic energy, four vorticity values, two PV values for the retained
  southwest boundary calculation, and six dual velocities reused by scalar work.
- Install results in the original traversal order, including later overwrites
  of shared vorticity. Add mass physics residuals with the original operation
  order. Bernoulli, boundary/pentagon processing, restriction and consumers read
  these actual production results.
- Skip the corresponding old interior arithmetic in production. The old formulas
  remain active for boundary work, the independent full oracle, and explicitly
  requested assignment checks. There is no steady-state reference packet stream.

The cache uses approximately 360 bytes per positive node (37 doubles and 16
default integers with the tested compiler), excluding allocation metadata. It is
not multiplied by 30 layers. New preparation, memory traffic and patch-buffer
costs are included in full-step timing; a speedup is not assumed.

## Verification and promotion gates

Local work uses serial builds and at most four MPI ranks. The required replays
are optimized assignment checks, checked assignment checks, retained full oracles
plus assignment checks, and optimized production with checks disabled. Each must
match the source-bound Stage266 checkpoint, step states and 62 VTK members.

Assignment checks compare 19 values bitwise per interior call, including mass
flux after physics. Per-rank evidence must cover all 48 native RHS transactions
and all 30 layers. Zero-Domain ranks must report zero primitive nodes, accurately
reflecting this candidate's unchanged execution ownership. Timing runs require
the production selector marker and forbid assignment logs/check mode.

The bbserv schedule uses four serial builds and sixteen solver runs on exclusive
nk/bb13 with established rank/core placement and local allocation:

| Runs | Purpose |
|---|---|
| 16-rank main and Stage266 | Science anchors |
| 16-rank optimized, checked and oracle candidate | Full assignment and scientific gates |
| 83/160-rank candidate; 256-rank oracle candidate | Ownership counts, all-call assignment checks and science |
| Main, Stage266 and candidate at 83/160; Stage266 and candidate at 256 | Eight single uninstrumented full-step timing controls |

All controls use the established fixture and ordinary-step exclusions. There
are no repeated identical timing configurations and no automatic promotion.
Compare both main and Stage266 below/equal to 160 Domains and Stage266 at 256.
Repeat only a marginal or unexpected result as needed. Resolve a repeatable
regression within this stage before stacking another change.

## What remains

Stage277 does not remove pressure/physics ownership, Domain-indexed restriction
and source workspaces, or native velocity tendency export. Stage276's resident
block input evidence remains useful for moving this integrated kernel later,
but does not authorize assuming that final field publication or every new halo
epoch is already correct. The next ownership migration still needs a closed
producer/restriction/consumer data path without a bulk intermediate return,
retirement of displaced production work, and immediate performance gates.
Useful >=640-rank execution remains a later multi-node qualification.

## Runtime and delivery

Allow approximately **20–30 minutes after allocation**, with uncertainty from
the full-call assignment checks and three extra timing runs relative to Stage276.
Stage276 measured 16:11 for four builds and thirteen runs. Queue and collection
are extra. Collection preserves exact Slurm accounting and failed-job evidence.

Use `test/parallel_block_profile/primitive_277/README.md` for copy/paste commands.
The preparation record records the actual local validation and archive digest;
no kit is ready until the production replay and fresh delivery checks pass.
