# Stage275: prove the inputs before moving the numerical work

Stage274 supports local source/restriction/gradient execution under the sampled
block ownership, but conservative Domain-fed input packing at 256 ranks costs
3.40–3.68 times the completed tendency export. It does not establish that this
packing is necessary. Stage275 tests two concrete ways to reduce it before
implementing a production replacement.

## Implemented changes

A diagnostic overlay on accepted Stage266 changes only `src/multi_level.f90`.
The working Stage273 candidate remains intact and excluded from the sealed kit.
The original numerical kernels and production outputs remain unchanged.

1. **Read-minimal live execution.** The exact `velocity_qperp` expressions use
   19 distinct components of each declared 3-by-9 flux/PV argument. The probe
   copies only those components and active-target physics inputs into scratch
   storage; unused locations contain quiet NaNs. It executes the existing ordered
   direct/restrict kernel and compares every source bit to the normal call at
   the same numerical phase. A corresponding gradient check supplies only the
   active gradient stencil nodes. Checks cover every call, every layer and all
   local Domains; zero-Domain ranks still produce complete phase records.
2. **Resident scalar epoch check.** At selected common transaction ordinals and
   the first/last atmospheric layers, record the actual block interior state,
   mean and sum alongside mapped gradient inputs read from Domain workspaces.
   Join by transaction, layer, global block ID, mesh level and exact coordinate
   bits. Report missing values and mismatches independently for all six fields.
   This tests current-state freshness, not merely whether an array is allocated.
3. **Revised transfer budget.** Recompute per-field minimal input budgets from
   sampled real stencils, deduplicated by field/component, destination and level.
   Retain the earlier full-argument budget for comparison. Do not subtract the
   separate scalar reuse samples as if they proved all halo inputs reusable.

This is a diagnostic prototype of the input contract. It still uses Domain-fed
shadow inputs and allocates scratch arrays for validation. It is not the compact
block-owned production executor and supplies no speedup measurement. No new MPI
exchange or generic per-value numerical interpreter is introduced.

## Decision the returned evidence enables

If owner-interior state/mean/sum values match, credit only those sampled mapped
interiors as candidates for local reconstruction. The next executor still needs
validated destination halo/alias routes, complete epoch coverage and compact
plans. If they differ, identify the publication/trial-state transition before
attempting reuse. Do not introduce a bulk copy to hide a failed epoch contract.

Flux is consumed at each level's exact restriction phase. Bernoulli and Exner
are consumed after final restriction; PV and physics retain their own producer
contracts. These upstream fields are not proven block-resident by this probe.
A smaller read footprint alone does not remove Domain-bound producer work or
prove faster execution. It quantifies the remaining producer/consumer boundary.

After input validity is established, implement a complete compact replacement
that retains useful producer fusion and retires displaced Domain computation and
export. Gate full-step timing against Stage266 and supported-rank main before
promotion. Eventual efficiency at >=640 ranks for J5's 160 Domains remains the
strategic target; this single-node diagnostic does not qualify it.

## Campaign and review

Use the established four-build/eight-validation Stage274 structure, with serial
builds and main/Stage266/checked/oracle/science gates retained. Run 16-rank main,
Stage266, candidate, checked, oracle; 83/160 candidate; and 256 candidate oracle.
The reviewer independently reconstructs both budgets and epoch comparisons from
raw files and verifies all rank phase sequences. Missing/corrupt evidence fails
review. Genuine mismatches are a reportable design result, not a reason to stop
collecting the remaining cases.

Selected transaction ordinals are 1, 2, 3, 9, 17, 33, 41 and 48. They are common
RHS-call counters, independent of the rank-local topology generations used by
Stage274. Only first/last atmospheric layers are exported. No claim is made that
these samples cover every adaptation/remap transition or prove halo freshness.
Unmapped boundary/scaffold inputs remain explicit. Duplicate spatial keys fail
rather than being silently merged. Stage274 structural sample ordinals remain
rank-local and are not a synchronized mesh census.

Estimated execution is 20–30 minutes, excluding queue and collection. Stage274's
combined runner time was about 14:10; shadow execution and field I/O introduce
additional uncertainty. Accounting is collected automatically. Commands are in
`test/parallel_block_profile/inputs_275/README.md`.

## Local validation

Two serial builds and three sequential four-rank replays (optimized, checked,
oracle) pass exact Stage266 checkpoint/state/VTK gates. The oracle retains 48
exact inverse comparisons and 84 oracle installations. Each replay checks
172,800 source calls and 57,600 gradient calls over all 48 RHS transactions and
30 layers. The 16 exported phase/layer combinations contain 339,836 mapped input
records, including 12,146 crossing Domain/block owner ranks; none are missing
or different. Budgets and input-epoch reports agree across all three variants.
Raw plan generations can differ because the oracle adds publication work.

After component and destination deduplication, local Domain-fed input budgets
fall from 71,056 to 68,048; 66,312 to 63,240; and 56,320 to 53,680 bytes per
layer for the three rank-local sample groups (4.23–4.69%). Do not extrapolate
these reductions into a runtime gain or to J5 at 256 ranks. The returned cluster
evidence must establish its own budgets and input-freshness results.
