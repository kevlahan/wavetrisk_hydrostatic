# Stage231 — remove false state prerequisites before the owner migration

Stage230 is committed as `4ea64012` on `codex/performance-investigation`.
Main remains `6a245db7f1632b0205c2bebe95e5f0cd535c2582`.
The larger objective remains removal of unnecessary Domain numerical storage
and execution pathways, including the complete producer/consumer cycle. The
changes below do not complete that migration and do not claim improved scaling.

## Locally validated changes

The complete retained tendency assembler no longer receives a Domain solution
argument in production. Its native velocity delivery, scalar restriction,
recomposition and retained-stage path propagate an absent optional argument.
Independent oracle calls retain their explicit reference argument. Entry guards
reject absent state in non-native/reference modes. This establishes an actual
argument-free production boundary; it does not eliminate geometry-owner native
velocity production or its still-required transport to block owners.

The persistent block hydrostatic cache is also removed from production
readiness, allocation and eager refresh at restart, grid changes, provisional
RK transitions, final prognostic refresh and vertical remap. Its three arrays
(surface pressure, dynamic Exner and air temperature) are consumed by explicit
validation/statistics code, not the numerical tendency kernel. Production now
requires these arrays to be absent. A request to construct them while disabled
fails with `block hydrostatic cache requested without a numerical consumer`.
The migration validation suite temporarily enables its diagnostic cache and
restores the production policy on exit; dynamics oracles retain the cache.
Readiness is not faked: absent storage remains not ready.

The live Domain primitive pressure calculation is unchanged. The separate
`native_columns` data used by `compute_block_exner_difference_kernel` are also
unchanged: those values participate in the accepted velocity add/subtract
recomposition and its roundoff. Removing that computation algebraically is a
different numerical change, not a proven dead-state deletion.

The previous prototype approach is not being reinstated. There is no new
Domain-shaped solution/RK-base snapshot, halo-cache variant or broad profiling
campaign. A speedup threshold is not used to decide whether these false
prerequisites should survive. No cluster kit is prepared for this cleanup yet.

## Local validation

All seven final cases passed against the fixed Stage230 J4/J6 fixture. The
optimized production runs at 4 and 41 ranks, optimized all-oracle run at 4 ranks,
and checked production run at 4 ranks each completed 12 steps, two remaps and
two checkpoint initializations. Checkpoint fields, step states and all 62 VTK
payload members match their same-rank Stage230 production references exactly.
The 41-rank case exercises a rank without a Domain. Production reports zero
Domain tendency allocations on every rank; oracle runs retain reference storage.

Two injected test-only executables explicitly reject forbidden hydrostatic
cache creation and non-native tendency assembly without Domain reference state.
A separate one-step control enables adaptation/bootstrap validation with the
dynamics oracle off; it completes through the restored production absence
guards. Its first reviewer incorrectly required verbose-only messages; the
corrected review checks emitted selectors, normal completion and the storage
census. This control is not an additional full numerical-equivalence run.

`parallel-blocks-stage231-local-results.json` records source, binary, input and
log identities and the complete case results. Reproduce its evidence review:

```bash
python3 test/parallel_block_profile/review_local_231.py
python3 test/parallel_block_profile/audit_domain_dependencies_231.py
```

These local, oversubscribed tests establish correctness coverage only. No
performance gain, full Domain-state removal or cluster qualification is claimed.
The final solver changes remain uncommitted after the Stage230 commit.

## Live dependencies and migration contract

`parallel-blocks-stage231-dependencies.json` anchors the reviewed ledger to 22
source routines in 9 families. It is a manual producer/consumer classification,
with automated source locations/hashes; it is not a dynamic read trace.

1. **Primitives and diffusion:** `cal_surf_press`, `integrate_pressure_up`,
   `basic_operators`, `step1` and diffusion prepasses still read Domain `q` and
   `sol_mean`. They produce real pressure/geopotential, Bernoulli, PV, Exner,
   viscosity and flux data needed downstream. Move their full columns and
   stencil dependencies to block owners. The restricted Exner consumed by
   velocity is not interchangeable with a separately reconstructed column.
2. **Mass/velocity and restriction:** `native_mass(d)` and `native_velocity(d)`
   still execute on geometry owners. Compile their ordered stencils, aliases
   and restriction delivery for the block owners; produce retained tendencies
   directly, then remove scalar capture/native velocity return transfers when
   their last geometry-side consumers are gone.
3. **RK and forward/inverse transitions:** retain the existing block timestep
   checkpoint as the base state. Replace the live Domain `q1`/`sol` consumers
   and canonical alias/scaffold versions together before deleting those arrays.
   Do not recreate their full shape under a new name or treat a resident
   snapshot as a complete producer/consumer port.
4. **Physics, soil, adaptation and I/O:** move the supported physics/soil state
   and ordinary-step consumers to block ownership. Bound any necessary external
   serialization adapters explicitly. Geometry/connectivity may remain Domain
   described; persistent prognostic replication must not be hidden there.

Acceptance requires complete ordinary steps, topology changes, remaps and
restart/output; independent intermediate oracle comparisons; absent-allocation
and rejected-fallback checks; and ranks with no Domains but owned blocks.
The final production cycle must execute with the old Domain numerical operands
unavailable. Interface optionality alone is not that proof.

## MPI scaling follow-up

The current catalogue splits a Domain root into its four child subtrees when
all four are present, independently of MPI task count. Thus block count greater
than Domain count is not itself a new mode at the high-rank transition. The
important additional limitation is that for P>D, at least P-D ranks own no
Domains. They execute block kernels but cannot share work in `do d=1,size(grid)`
primitive, diffusion, geometry-native mass/velocity and RK loops.

For the existing Stage226 samples, only 40 of 160 ranks can execute the Domain
loops in J4; at most 160 of 256 can in J5. These are worker-count bounds for those
phases, **not** measured timestep fractions or predicted total efficiencies.
The existing strong-scaling screen was J4 0.830→0.651s (40→160 tasks) and J5
0.600→0.657s (160→256), with one production run/count. The observations are
consistent with this ownership limit plus transfer, synchronization, imbalance
and memory-bandwidth costs; they do not apportion those causes.

The block allocator balances a leaf-cell proxy, not measured computational
cost. `subtree_weight_Domain` sums `PATCH_SIZE**2` per live leaf patch; it does
not weight internal-level kernels, halo sizes, message counts or phase timings.
The ordered next-fit assignment begins with 5% overload tolerance and multiplies
that tolerance by 1.2 until every indivisible block fits. Thus 5% is an initial
target, not a guaranteed final bound. As tasks approach block count, unequal
subtree weights cannot be smoothed by assigning several blocks to a rank. At
the J4 160-task/160-block boundary, each rank has exactly one block: changing
the assignment cannot remove unequal subtree sizes. Splitting an oversized
subtree further would address that constraint, but not Domain-owner work.
Even ideal balance in this proxy would leave the Domain-owner phases above.
The existing optional catalogue report already prints minimum/mean/maximum
rank weights; this establishes proxy balance only, not timestep balance.

After migrating the joined cycle, revisit optimization with fixed workload,
block topology and controlled placement. Use the existing high-precision
ordinary-step medians and independent ABBA runs. If phase evidence is needed,
collect a bounded set of per-rank compute/wait durations and ownership/work
counts, buffered outside timed sections; distinguish critical-rank elapsed
cost from MPI wait inherited from upstream imbalance. Do not compare nk/bb as
one scaling curve, or assume deeper splitting fixes owner-limited execution.
