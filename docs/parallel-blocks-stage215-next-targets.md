# Remaining targets after Stage215

Source audited at `98d940c2` on `codex/performance-investigation`.
Stage215 is promoted: the applied `src/comm_mpi.f90` is byte-for-byte the
validated candidate. Its corrected 83-task confirmation is 0.931 to 0.814 s
(117 ms, 12.57%). All 22 campaign helper tests pass. Main remains
`b82647b36eb825aed0c0aae4e92adcf4ab46d608`.

The useful lesson is to compile repeatedly interpreted routes and remove
unnecessary staging on actual production paths. Stage215 does not establish
how much each individual change contributed, or the post-change component
budget. No further optimization is implemented by this note.

## 1. Closest remaining analogue: ordered same-rank Domain halo plans

`src/comm.f90:1266` (`cp_bdry_inside_0/1/2`) still scans every local
source/destination Domain pair, resolves global Domain IDs, loads pack/unpack
lists and decodes signed destinations for each field/layer. The Stage215
remote paths still call these original routines (`src/comm_mpi.f90:1231`,
1334 and 1441).

A concrete candidate is to compile the nonempty local Domain-pair groups and
their source/destination indices, with sign information, when communication
routes change. Reuse these groups across fields and layers, just as Stage215
does for remote packing and installation. Keep the values live: this is an
address plan, not a cached halo result.

Correctness obligations:

- Retain original source-Domain, destination-Domain, field/layer and element
  ordering. Ordered local copies can read values written by earlier copies.
- Retain all-level local copying even for a single-level remote request.
- Preserve signed edge negation, including signed zero, and all overloads.
- Invalidate at every relevant topology/route rebuild and initialization;
  validate actual lists in the checked build.
- Keep remote snapshots before local copies and remote installation afterward.

This is the most direct extension of the successful intervention, but **not
yet a demonstrated 100 ms target**. At 16 tasks there are roughly ten Domains
per rank; at 83 tasks only one or two in this 160-Domain workload. The eliminated
pair scans therefore have a particularly strong risk of failing to transfer
to 83 tasks. Copies themselves remain necessary. The historical whole Domain
boundary timer includes remote work already improved by Stage215 and must not
be presented as this candidate's budget.

Recommendation: first use a bounded local implementation/replay comparison
with representative one/two-Domain and ten-Domain layouts. Stop if execution
cost simply moves into plan storage or index loads. Only a convincing result
would justify the existing 16-task screen and immediate 83-task transfer check.
Do not schedule a broad observer campaign merely to obtain another timer.

## 2. Real repetition, limited scope: shared primitive/stencil geometry

The production caller (`src/multi_level.f90:315`) loops over 30 vertical layers.
Stage215 run logs confirm `zlevels=30` and all three Laplacian settings equal
2; these diffusion paths are active in the tested workload.

- `basic_operators` (1027) calls the generic pressure traversal and `step1`
  for every patch/layer. `step1` (`src/ops.f90:35`) recomputes `comp_offs3`
  offsets/dimensions and boundary classification for each call.
- `cal_Laplacian_scalars` (1582), `cal_Laplacian_vector_rot` (1635), and
  `cal_Laplacian_divu` (1670) revisit the same patch geometry for diffusion.
- `apply_onescale_to_patch` (`src/domain_ops.f90:352`) rebuilds offsets,
  dimensions and boundary extents, then invokes a procedure for each point.
  Callbacks such as `cal_Laplacian_rotu` resolve multiple `idx` addresses.
- Bernoulli and Exner each call `cpt_or_restr_scalar` (1741), traversing the
  same parent/child geometry and recomputing injection indices.

Possible removal: share patch offsets, boundary extents and ordered restriction
addresses across these calls. Start with compact per-patch metadata; do not
assume that a large per-point stencil table improves cache behavior. Preserve
arithmetic order and vertical pressure dependencies. Masks can change without
ownership changing: retain live mask predicates or explicitly validate the
mask state. Topology changes require rebuilding addresses.

**Do not make this the next standalone 100 ms campaign.** The corrected
Stage207 RK shared primitive/operator group has 247.007 ms inclusive time,
but only **45.315 ms exclusive local time** (legacy 47.428 ms). Its children
include boundary exchanges, capture and native mass work, which cannot be
credited to geometry caching. This is a real potential modest improvement,
not evidence for removing a quarter second or explaining a block-only excess.
These are historical pre-215 rank-mean measurements, not current bounds on
whole-step effects through synchronization.

## 3. Conditional exchange elimination: native mass edge halo

`basic_operators` ends its native mass restriction with
`exchange_native_mass(AT_EDGE,l,k)` when multiple grid levels are active
(`src/multi_level.f90:1121`). Its caller immediately performs another exchange
at 371 before divergence; the intervening optional rotational Laplacian uses
other storage.

This is worth proving redundant, not simply deleting. A second local copy
pass can propagate values further along alias chains; remote snapshot timing
can also matter. Demonstrate that the first exchange leaves every source and
destination in the exact state required by the second, across topologies and
local/remote ownership, before removing either call. Native mass transport
totaled 60.294 ms in the corrected Stage207 accounting, including all its
node and edge exchanges. Only part belongs to the potentially redundant call.
This therefore lacks support as a standalone 100 ms intervention.

## Paths checked and ruled out as fresh large targets

- Native mass and velocity already compile their numerical and boundary
  address programs and compare actual mask predicates. Their historical
  plan costs are about 0.314 and 0.260 ms, respectively.
- Stage202 inverse aliases already resolve route addresses outside the
  field/layer loops and bypass intermediate buffers for self copies.
- Scalar producer interior recursion is oracle-only. Production capture
  uses narrowed boundary records; full geometry is not resent every RK stage.
- Further scalar/inverse field batching is the tested Stage203 proposal;
  its sub-threshold gain is not a newly discovered opportunity.
- Direct local block ghost copies could remove intermediate scratch writes,
  but their whole historical scopes are only tens of milliseconds.
- `comm_nodes3_mpi` and `comm_nodes9_mpi` retain per-element buffer extension
  and receive clearing, but their call sites are coordinate/mesh refinement
  paths, not the repeated 30-layer RK field halo. Do not infer Stage215-sized
  gains merely from similar-looking code.

## Decision and evidence

There are concrete remaining analogues, but this audit has not established
another 100 ms saving. Give the same-rank halo extension a small, bounded
prototype budget; keep geometry reuse and proven exchange elimination as
secondary opportunities. Do not revive the rejected inverse/reconstruction
changes or claim that their previous costs are still recoverable.

Any eventual candidate must compare against committed Stage215, retaining
Stage202 on, Stage203 off and FP contraction off. Preserve correctness first,
then 16-task A/B, early 83-task A/B, and fresh A–B–B–A confirmation. Apply the
corrected timing rules: exclude saving and first full integration steps after
restart; keep ordinary remapping. Timesteps are not independent run repeats.

Evidence used: `parallel-blocks-cluster-215-54714.md/json`, Stage215 raw run logs,
and `parallel-blocks-evidence-reconciliation-210.json`, specifically
`campaigns.207.operation_differences.common` (context 2, operations 2, 3, 5,
6 and 23). The historical reconciliation uses steps 2–7 and 10–12 and charges
MPI once to its immediate owner. Its exclusive times and nested inclusive
times must not be added together. Earlier Stage214's broad dismissal of a
generic cache opportunity is superseded by the demonstrated Stage215 result;
future decisions should name and inspect the actual production path.
