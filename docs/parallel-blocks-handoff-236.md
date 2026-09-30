# Stage236 handoff — RK migration complete

The restarted climate RK4 route with Simple physics now closes every RK stage
from retained block state and registered, versioned aliases. `BLOCK_RK` remains
off by default and requires `BLOCK_TENDENCY_INITIAL=on/check` when enabled.

Authoritative worktree: `/Users/kevlahan/.codex/worktrees/7173/wavetrisk_hydrostatic`.
Branch: `codex/performance-investigation`, HEAD
`4ea640128a71a857f9b8c6b9411dc34fc0aeb641`. Main remains
`6a245db7f1632b0205c2bebe95e5f0cd535c2582`. Changes are uncommitted and unstaged;
no remote jobs were created. The automatically created 26d8 checkout was not
used as implementation source.

## Implementation

Mass tendencies stay on final block owners. Canonical velocity aliases use the
existing block checkpoint plus exact retained gradients in the original RK
rounding order. Sparse temperature results and versioned noncanonical aliases
feed the joined forward/inverse boundaries. Provisional low scaffold, fixed
coarse scaffold and final boundary epochs preserve their distinct lifetimes.

The production route no longer allocates q1 numerical arrays or native
mass/velocity tendency and RK return arrays. It bypasses Domain RK arithmetic,
interstage temperature publication, legacy stage boundary updates and provisional
inverse publication. Accidental use of the removed RK closure is rejected.
Independent oracle arithmetic and transform comparisons remain available.
No full Domain-shaped backup was introduced.

## Qualification

Final acceptance PASSED. `parallel-blocks-stage236-local-results.json` is generated
by `test/parallel_block_profile/rk_236/review_rk.py` and binds all 69 frozen build
inputs, 115 resolved sources, four binaries and the following evidence:

- Full 4- and 41-rank runs: 12 exact step states, two remaps, restart/checkpoint
  fields and all 62 VTK members match the same-rank references exactly.
- Rank40 has zero Domains and executes retained RK work on its owned blocks.
- Checked native and independent oracle runs pass. All four RK boundaries
  compare original registered aliases bit-for-bit; existing transform checks pass.
- Old registered physical RK inputs remain NaN across all four stages in the
  poison test, without a backup or restoration; the completed step is exact.
  q1 and native return-storage absence is checked throughout production runs.
- Missing mass output, stale stage, missing base checkpoint, corrupted velocity
  producer and forbidden Domain fallback each hit their intended rejection.
- The RK feature-off control passes and matches the reference step state.

Frozen sources/builds are in `.stage236-validation/qualified-candidate`. Only
`qualified-*`, `control-4`, `poison-4` and `reject-*` cases bound by the final report
are acceptance evidence. Development failures and the superseded
`pre-boundary-candidate` remain separate. Stage235 frozen sources and binaries
were reverified unchanged. The earlier preflight handoff is archived as
`parallel-blocks-stage236-preflight.md`.

## Remaining scope

Initial physical inputs and final publication for physics, adaptation, restart
and output remain explicit lifecycle boundaries. Upstream physics/soil, Domain
topology, adaptation/restart storage and wavelet/reference aliases need their own
migration before deleting those fields. This result completes the qualified RK
boundary, not all Domain storage removal. Block splitting/load balance is a
subsequent workstream.

There is no performance or scaling claim; keep the selector disabled by default.
No broad cluster timing campaign is justified by this correctness result alone.
Future bb partition work must continue to avoid bb04. See
`parallel-blocks-stage236-rk-migration.md` for implementation details.
