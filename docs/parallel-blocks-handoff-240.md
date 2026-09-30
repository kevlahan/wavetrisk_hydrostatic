# Stage240 — owner physics dependency scheduling

All ten qualification cases and `physics_240/review_physics.py` passed. The
source-bound results are recorded in `parallel-blocks-stage240-local-results.json`.

`PHYSICS_OWNER_BATCH=on/check` replaces the source-rank serial request/reply
loop with bulk owner import, dependency-wave execution, and bulk publication.
Its default is off; it requires the existing owner and velocity-plan selectors.

Each column carries six storage identities and their expected writer ordinals.
Version zero retains the entry value; positive versions receive the exact edge
component returned by the preceding writer. Only integer topology is replicated.
Numerical input packets, hydrostatic scratch, package work, and returned state
are held by execution owners during the physics epoch. Independent columns can
execute in the same wave. The tested trajectory uses 35 or 38 waves per step.

The production publication callback writes complete double-precision temperature,
velocity and soil results without reading Domain numerical operands. Entry export
and final publication still use Domain compatibility storage. Persistent owner
state across soil transforms, adaptation and I/O remains outstanding; this change
does not complete the migration or establish a performance improvement.

The check build replays the original source-column order and compares complete
single-precision package packets and double-precision state bit for bit. The queue
pass consumes the legacy hydrostatic scratch, so replay explicitly resets surface
pressure and geopotential before starting. Failed pre-fix snapshots are preserved
under `.stage240-validation/candidate-attempt2` and `oracle-4-attempt2`.

Completed checks:

- Independent four-rank replay: 35,812 exact column comparisons.
- Full four-rank and 41-rank runs: all 12 step states, two remaps and 62 VTK
  members exactly match Stage239 at each rank count.
- Each full run executes 460,608 owner columns. Remote counts are 14,796 and
  210,104 respectively; zero-Domain rank 40 executes 5,824 remote columns.

The serial control matched exactly. Wrong routes, missing edges, future versions
and corrupted owner results were rejected at their intended guards. Destructive
source-field poison required no backup/restoration and matched the first reference
step. Package-soil poison with forced late initialization matched two reference
steps. The reviewer also verified unchanged Stage239 evidence and main.

The local host has four usable cores. Solver cases run strictly sequentially with
library threads set to one. The 41-rank run is deliberately oversubscribed and is
correctness coverage only. No bbserv run has been launched.
