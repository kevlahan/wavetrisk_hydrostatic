# Stage274: execution contract before a new ownership migration

This begins the new architecture program authorized on 2026-10-08. It is a
source-derived ownership and data-movement diagnostic, not a numerical migration
or a timing candidate. It starts from accepted Stage266 `ff740744608731a05d53771c32a3055309f3a586`.
Stage273's unqualified partition change is excluded from all three source trees.
Its separate kit and working source are preserved.

## Why this is the first change

The old Stage251 ladder identified a confirmed 386.426 ms velocity-owner penalty.
Stage252 located about 89% of the measured all-step regression in source execution.
Stage254 exposed substantial repeated producer traversal and positive-depth
MPI_Alltoallv work; MPI duration included waiting. Skipping empty collectives
and indexing the traversal recovered only part of the penalty. The historical
reports remain on `codex/performance-investigation`; these measurements do not
prove that every MPI millisecond was recoverable.

Current Stage266 already has direct, Domain-independent numerical routines in
`src/parallel_block_velocity.f90`. The problem is not the absence of such kernels:
`prepare_native_velocity_programs`, `compute_native_velocity_source` and
`compute_native_velocity_gradient` in `src/multi_level.f90` allocate and execute
those kernels by `size(grid)`. Flux, PV, physics, Bernoulli, Exner and density
inputs still belong to Domain-owner workspaces. Shipping these inputs to a
velocity-only block consumer could cost more than the current tendency export.
A new generic owner program would repeat the failed mechanism.

The immediate implementation therefore exports actual prepared stencils, target
ownership and ordered restrictions. An offline compiler tests the proposed cut
in the dependency chain before numerical replacement. This supplies a concrete
before/after data-movement budget; it deliberately supplies no estimated speedup.

## Complete target chain and retained dependencies

| Transaction | Current execution/data | Replacement contract | Work retired only after replacement |
|---|---|---|---|
| Pressure and shared primitives | Domain `integrate_pressure_up`, fused `step1`; produces Bernoulli, Exner, PV, vorticity, kinetic energy, flux inputs | Direct patch/column batches on block owners, preserving useful fusion, pole/boundary rules and exact arithmetic epochs | Corresponding shared production and capture, after every live consumer has a replacement |
| Mass flux and restriction | Domain-indexed native mass arrays; level-specific edge/node exchanges | Block-owned flux plus bounded multilevel halo/restriction schedule | Native Domain mass staging/exchanges superseded by that schedule |
| Velocity physics and source | Domain callback inputs plus ordered direct/restrict program | Batch physics adapter or native physics interface, direct source kernel, explicitly scheduled restriction boundaries | Domain source storage and repeated physics bindings no longer required |
| Gradient and final tendency | Domain Bernoulli/Exner/mass/temperature/source; exported completed tendency | Direct local stencil gradient, writes block-owned tendency consumed by block RK | Native tendency export and corresponding Domain numerical storage |

These rows form the architectural target, not a promise to migrate all of them
in one patch. A safe smaller cut must prove that the replacement does not add a
bulk intermediate return route. Stage267's canonical publication and boundary
refresh contracts remain mandatory until all their live readers are replaced.
Necessary geometric identity and explicitly bounded scaffold adapters may remain;
the diagnostic exposes them rather than inventing a block owner.

## Implemented diagnostic and planning prototype

`test/parallel_block_profile/ownership_274/instrument.py` changes only a disposable
copy of Stage266's `src/multi_level.f90`. It exports the first three actual native
velocity plan builds for each geometric Domain. It records the actual topology
and mask-selected actions, generation, ordered direct/restriction writes, gradient
stencils and node-to-current-block/rank mapping. It reads no numerical fields and
adds no MPI collective. Each rank, including ranks without Domains, writes a
header. Existing numerical computation and all science/oracle gates remain.

Subtree interiors are mapped recursively from the real block catalogue. Nodes
outside those interiors keep owner `-1`; these include boundary/scaffold storage
and require a separately proved adapter or canonical ownership map. Diagnostic
files cannot be interpreted as a globally deduplicated physical-node catalogue.

`contract.py` is an offline design compiler, not a numerical runtime interpreter.
It follows the actual write sequence by source component. It records each remote
restriction operand's producing write/version, flags same-level cross-owner
read-after-write dependencies, and flags source reads not initialized by the
sampled sequence. A direct write followed by partial restriction is preserved;
merging all writes by target would hide a real dependency.

The report contains:

- Actual prospective direct/restriction/gradient work by rank, including adapters.
- Full direct-kernel input packing budgets, deduplicated per destination and level:
  flux and PV use the declared three-by-nine arguments; physics uses three values
  per active target. These are conservative protocol budgets, not executed-load counts.
- Gradient input budgets for Bernoulli, Exner, density and density-temperature.
- Remote input bytes if upstream producers remain on Domain owners, and a clearly
  hypothetical budget if mapped upstream fields also become block-owned.
- The current completed-tendency export budget: three doubles for each mapped
  remote interior node per layer. This makes a velocity-only relocation's added
  input traffic visible against the output traffic it could retire.
- Versioned restriction transfer counts, same-level scheduling blockers, unresolved
  source reads and all unmapped adapter work.

Budgets exclude headers, transport startup, plan construction and bandwidth/cache
cost. Boundary aliases use separate Domain-local identities; canonical alias
resolution could change the budget. The hypothetical upstream placement is not
an implemented producer. Sample ordinals across Domains need not be simultaneous,
and samples are not weighted by RK/level invocation frequency. Multiply neither
sample totals nor per-layer budgets into claimed whole-step savings.

The retained-Domain-input scenario is not a lower bound on migration traffic.
Some inputs may already be block-resident or cheaply reconstructed from block
state. The prototype must prove their value epochs before crediting that reuse;
this conservative first budget does not model it.

Absence of the detected source blockers is a necessary test only. It does not
prove a complete communication schedule, correct boundary epochs, or a speedup.
The report always marks production readiness and performance qualification false.

## What changes after this diagnostic

Use the observed dependency cut to choose a bounded direct block prototype.
Prefer a complete fused producer/consumer slice whose large intermediates remain
local. Compile compact integer stencils and communication lists when topology,
mask or ownership changes; validate all three before reuse. Keep numerical array
access direct and preserve within-target arithmetic/write order.

Before production activation, the prototype must compare exact intermediate
results against the existing independent reference on real adaptive inputs.
Measure construction, pack/unpack, numerical loops, communication and final
publication together. Do not infer net benefit from a fast isolated kernel.
Retire the displaced production path in the same qualified change; keep the
reference path test-only. If a cut requires expensive bulk copies or repeated
same-level remote stages, redesign the cut before implementing a general scheduler.

The scientific gate includes optimized/checked/oracle runs, restart/remap,
boundary epochs, checkpoint/state/VTK comparisons and ranks without Domains.
The later performance gate uses uninstrumented complete steps against Stage266
and supported-count main, with 83/160/256 as an intermediate check and >=640
ranks/multiple nodes as the strategic target. A repeatable regression blocks
further migration. Finer block granularity must be qualified separately.

## Campaign

Four serial builds share the sealed 23 physics objects: accepted Stage266,
legacy main with established checkpoint I/O hardening, diagnostic and checked
diagnostic. Eight full J5 runs use one exclusive nk/bb13 allocation:

1. 16-rank main and Stage266 correctness anchors.
2. 16-rank diagnostic, checked diagnostic and diagnostic oracle.
3. 83- and 160-rank diagnostic topology samples.
4. 256-rank diagnostic oracle, including ranks without Domains.

There are no timing comparisons and no repeats. Printed solver times are retained
only as raw diagnostics; file output perturbs them. Runtime/accounting collection
is automatic. Estimate 12–18 minutes after allocation, using Stage272's 18:52 for
four builds/fourteen runs as the nearest measured basis. Export/postprocessing
cost is new, so this estimate is less certain than the earlier timing campaigns.
The kit is self-contained except for the established sealed J5 fixture; launch,
duplicate prevention, placement and failed-job collection use the proven workflow.

See [commands](../test/parallel_block_profile/ownership_274/README.md). Local validation
and final delivery hashes are recorded in the preparation record when complete.

## Local qualification and first evidence

Both serial builds and all three sequential four-rank J4 replays pass exactly:
optimized, checked and oracle numerical records, twelve step states and 62 VTK
members. The oracle retains 48 exact inverse-return proofs and 84 installations.
Each case exports 120 real plans, three for each of 40 Domains. Optimized and
checked structural reports agree exactly. Nine workflow and six contract tests
pass, including missing/corrupt data, no-Domain ranks, ordered partial writes,
same-level cross-owner blockers and explicit adapter handling.

For the first sampled J4 plan set, the conservative retained-Domain input scenario
has 71,056 remote bytes per layer, versus 23,040 for the current completed-tendency
export. All three sampled sets show zero detected same-level cross-owner source
dependencies and zero uninitialized source reads. The first set still has 25,434
source-input and 9,092 gradient-input values outside mapped block interiors.
These are structural observations, not runtime measurements or a J5 conclusion.
They justify testing the dependency cut before choosing a velocity-only migration;
the actual prototype must account for local field reuse and complete upstream costs.

The freshly extracted final kit verifies 48 sealed files, all source families,
exact reproduction of the diagnostic overlay, exclusion of Stage273, imports,
shell syntax and all fifteen delivered tests. It contains no macOS xattr headers.
See [the preparation record](parallel-blocks-stage274-preparation.json).
