# Stage234 handoff

Working branch: `codex/performance-investigation`; HEAD remains Stage230
`4ea640128a71a857f9b8c6b9411dc34fc0aeb641`. Main must remain untouched at
`6a245db7f1632b0205c2bebe95e5f0cd535c2582`. Stages231–234 are uncommitted.
The workspace also contains many unrelated historical untracked campaign
files: never stage everything.

The latest request was to implement the next stage to complete full tendency
migration. Stage234 makes actual additional owner connections, but that full
request remains incomplete. Read `parallel-blocks-stage234-owner-primitives.md`
and the source-bound `parallel-blocks-stage234-local-results.json`; do not
present this as Domain-free dynamics or a performance improvement.

New modules are `parallel_block_primitive_walk`, `parallel_block_resident_operands`,
`parallel_block_primitives`, `parallel_block_pentagons`, `parallel_block_gradients`.
The Stage232 pure arithmetic module is unchanged. Stage233 columns now expose
owner result addresses. The inverse registers the exact regular stencil closure
and exposes vector operands from the provisional epoch. `ops`, `multi_level`,
`parallel_block_mpi` and `time_integr` contain the live dispatch/retention hooks.

Important discovered ordering details:
- Pressure visits and primitive visits differ. Record pressure visits before
  each patch's regular traversal; do not count regular visits as pressure writes.
- Some extra boundary regular visits have no preceding pressure producer.
  Their B writes are unused. Per-layer definition checks cover the actual B
  restriction and gradient readers before omission is permitted.
- SW consumes its own preceding regular PV context. Pentagon corrections consume
  last regular/SW vorticity writes; do not read the later diffusion overwrite of
  Domain vorticity. Cross-pentagon dependencies are currently rejected explicitly.
- Bernoulli/Exner restrictions propagate producer addresses, so final gradients
  fetch the actual restricted producer values from their owners.
- The retained gradient map covers every compact node, including explicit zero
  definitions for inactive nodes. Active outputs are retained across layers on
  block owners and copied directly into block_velocity_remainder. The old native
  geometry-owner capture is rejected once those owner results are ready.

BLOCK_PRIMITIVES and BLOCK_COLUMNS remain off by default. The enabled route
handles stages after a provisional inverse. Stage-one/restart/final-transform/
post-physics reference lifetimes remain separate work. Current aliases still
use geometry-owned reference means; pressure still imports surface reference
values. Primitive compatibility outputs/returns remain for diffusion/source
consumers, and the explicit native velocity return remains for RK.

Do not resume broad instrumentation, local-halo cache variants or small timing
campaigns. The next coherent implementation is diffusion prepasses/exchanges
and their scalar/velocity physics consumers, then the mass and ordered velocity
source programs. Remove their old captures/returns when their last readers move.
Do not use a full Domain-shaped solution snapshot as a shortcut. Preserve exact
intermediate oracles and full adaptive/restart checks; then disable former
numerical inputs before deleting storage. Return to scaling and block splitting
after that dependency removal. Exclude bb04 from any future bb campaign.

Qualification uses `.stage234-validation/qualified-candidate`, with helpers in
`test/parallel_block_profile/primitives_234`. Earlier candidate/final-candidate
builds are prototypes. The earlier `full-checked-4` was interrupted deliberately
when direct retention was added, and is not a completed validation case.
No remote jobs, branch changes or commits were performed.

Final qualification: all ten checks passed and the review was generated. Both
full trajectories (4 checked, 41 production) matched checkpoints, step states
and 62 VTK members exactly through 12 steps, two remaps and restart. The complete
one-step oracle, feature-off control and six fault injections also passed.
All local validation processes have completed; no jobs remain running.
