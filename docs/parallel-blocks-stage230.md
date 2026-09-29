# Stage230: remove the production Domain tendency storage family

The revised objective is to eliminate unnecessary Domain numerical storage and
move complete producer/consumer chains to block ownership. Architectural
progress no longer depends on exceeding a fixed timing threshold. Geometry and
connectivity may remain Domain-described where useful; numerical replication
and geometry-owner execution must not be mistaken for final block ownership.

This change starts from `5c457a5dc869ce18938dd24db5727053b8e8dc2a` on
`codex/performance-investigation`. Main is not checked out or edited. Its ref
on the September29 continuation is `6a245db7f1632b0205c2bebe95e5f0cd535c2582`,
including the user's separate commits made during the pause.

## Implemented deletion

The qualified restarted climate RK4 / compressible / unsplit / Simple-physics
path, without NCAR topography or SSO, no longer allocates the Domain `trend`
array. This covers all prognostic tendency components and vertical layers,
including field descriptors and numerical payloads. Initialization, refinement,
restart reconstruction and teardown respect its absence. Production timestep
guards reject a recreated array or a fallback to nonresident dynamics.

`prepare_block_multistage_tendency` passes no Domain tendency argument in
production. The shared primitive interface permits an absent argument only in
the native mass transaction; it rejects a request for the legacy path. Native
mass, velocity and temperature producers continue to feed the RK update.
Production RK closure no longer traverses the generic legacy RK tendency
consumer with skip flags. The untouched temperature scaffold is copied
explicitly; native publication supplies its evolved values.
Primitive `step1` calls also omit the unused tendency argument; only the oracle
restriction path within `basic_operators` can require a Domain tendency.

Legacy configurations, independent dynamics oracles and historical poison-test
builds retain their reference arrays. The legacy tendency entry point explicitly
rejects a call when storage is absent. This preserves an independent numerical
comparison while removing the production allocation itself.

VTK's vertical-velocity calculation was a hidden consumer of the same global
array. It now allocates only its own mass-divergence scratch for the output
transaction, initializes its shape without copying prognostic values, and frees
it immediately after the last divergence consumer. It cannot recreate `trend`.
This is an explicit remaining output compatibility boundary, not block-native
output.

## Remaining storage and ownership work

| Storage/path | Current consumer requiring it | Required next change |
| --- | --- | --- |
| Domain `sol` and RK `q1` | `prepare_block_multistage_tendency`, `basic_operators`, RK base/closure and provisional wavelet interfaces | Move the entire inverse → primitives → tendencies → RK → forward-transform cycle to block ownership; delete the stage field only when all these consumers have moved. |
| Domain `sol_mean`, pressure, Bernoulli, vorticity and diffusion scratch | Primitive and external-physics calculations | Give each producer and downstream consumer block-owned inputs and outputs; preserve column dependencies and halo semantics. |
| `native_mass(:)` / `native_velocity(:)` geometry-indexed workspaces | Native tendency/RK calculations currently indexed by local geometry Domain | Change execution ownership and indexing as part of the same chain; renaming these arrays does not remove Domain ownership. |
| Domain wavelet fields and canonical inverse workspace | Reconstruction, provisional/accepted wavelet transitions and adaptation | Make block-resident accepted/provisional fields authoritative through forward/inverse and topology changes, with explicit generation checks. |
| Domain SOL publication/writeback | The consumers above, external physics and I/O | Remove each transfer only after its last consumer has moved. Retain bounded adapters solely at explicitly supported external interfaces. |
| Legacy `trend` in other configurations | Legacy integrators, diagnostics and independent oracles | Retain as compatibility/reference code until those configurations are migrated or explicitly retired. |

The next architectural unit is the full block-owned dynamics cycle in the first
row, including the shared primitives it invokes. An isolated resident inverse
snapshot would leave the ownership problem intact. Stages218/219 establish
useful numerical checks, but their Domain-shaped snapshots and geometry-owner
execution are not a production replacement.

The horizontal-flux arrays cannot simply be removed alongside `trend`:
`block_tendency_compatibility_ml` calls `cal_Laplacian_scalars` before the native
tendency transaction, and that diffusion prepass still uses Domain flux scratch.
The existing poison tests deliberately poison it only after this prepass.
`step1` also reads Domain-addressed state for the scalar physics residual and
uses Domain geometry for pressure/vorticity/Bernoulli calculations. These are
concrete remaining consumers to migrate with the full cycle, not evidence of a
removable timing budget.

Keep the existing intermediate mass/temperature/velocity comparisons, adaptation
and remap oracles, output/restart checks, and zero-Domain receiver coverage.
Every deleted storage family needs an absent-allocation assertion and a rejected
fallback test. Do not add broad timing instrumentation or another halo-cache
variant. Test architecture and correctness first; use the established 16 → 83 →
ABBA protocol to investigate performance and any material regressions.

## Validation

The repeatable local driver is `test/parallel_block_profile/run_local_230.py`.
It checks all twelve J4 timesteps against the accepted Stage229 production
reference, both checkpoint reloads, both vertical remaps, an exact final
checkpoint and exact decompressed VTK output members at both save events.
It also requires the correct storage policy on both initializations, including
an MPI sum proving allocation is absent on every production rank. VTK comparison
requires a reference at the same MPI rank count; the driver checks this before
launching a run.

The September22 temporary builds, seed and raw logs were purged during the
pause. The September29 continuation regenerates an evolved J4/J6 seed with the
pinned Stage229 production reference, using the established spin-up input. It
then reruns the matched reference and candidate gates. This is a new local
fixture identity, not a substituted cluster fixture or a claim to have recovered
the earlier raw evidence. Build trees, fixture and raw outputs are kept under
`.stage230-validation/20260929/` in the investigation worktree, outside system
temporary storage. The directory is ignored by Git.

`test/parallel_block_profile/reject_fallback_230.py` builds an isolated test-only
main module that deliberately requests the legacy consumer on the first
production step. It requires a nonzero MPI exit and the explicit absent-storage
diagnostic, and verifies that production source is unchanged.

The first 41-rank output review compared against a four-rank reference and
correctly flagged different VTK bytes despite an exact checkpoint. A matched
41-rank run of the unchanged reference reproduced the output-layout difference:
the exporter deduplicates vertices locally and gathers those local lists. For
example, one layer has 5,552 exported points at four ranks and 6,090 at 41.
All 62 candidate VTK payloads match the 41-rank reference exactly. The baseline
checkpoint also matches across four and 41 ranks. The fix is a matched-rank VTK
control, with no tolerance relaxation or solver change.

All final local gates passed on September29; results and source/binary identities
are recorded in `parallel-blocks-stage230-local-results.json`.

| Gate | Result |
| --- | --- |
| Optimized production, 4 ranks | Twelve matching step records, exact checkpoint, 62 exact VTK payloads; storage count 0/4 on both initializations. |
| Independent dynamics/adaptation/remap oracles, 4 ranks | Same exact final comparisons; reference storage count 4/4 on both initializations. |
| Optimized production, 41 ranks / 40 Domains | Same exact comparisons against the matched 41-rank reference; storage count 0/41 on both initializations. |
| Checked production, 4 ranks | All twelve steps, reloads and remaps pass with bounds checks and floating-point traps; exact checkpoint and VTK payloads; storage count 0/4. |
| Deliberately injected legacy fallback | Nonzero MPI exit with `legacy tendency requested without Domain tendency storage`; production source unchanged. |
| Mismatched-rank VTK reference | Rejected before launching a solver run. |

The unchanged reference also has an exact checkpoint and matching step records
across four and 41 ranks. All 23 external-physics objects match between baseline
and candidate; the pinned baseline's tracked source is verified unchanged.
Candidate source identities cover 144 source/build-input files. Optimized and
checked builds use the validated Stage202/203 defaults with FP contraction off.

Local MPI timings are not evidence of a speedup. This stage is a complete
deletion of one production storage family, not completion of Domain-free
dynamics. J5 cluster qualification subsequently passed in job55089: all15 runs
completed, including167 tasks /160 Domains. The83-task ABBA average was
0.757403s baseline versus0.754893s candidate, within the3% regression guard
but with opposite-sign paired savings. No demonstrated speedup is claimed.
See `parallel-blocks-stage230-review.md` and `.json`. Changes remain uncommitted.
