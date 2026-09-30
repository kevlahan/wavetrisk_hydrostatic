# Stage236 — retained RK closure

Qualification passed. `parallel-blocks-stage236-local-results.json` binds the
frozen sources and binaries to exact 4/41-rank trajectories, independent oracles,
checked storage absence, poison/fault tests and the feature-off control. See
`parallel-blocks-handoff-236.md` for the completed handoff.

`BLOCK_RK=on` (or `check`) requires `BLOCK_TENDENCY_INITIAL=on/check`. The default
remains `off`. The qualification target is restarted climate RK4 with Simple
physics and resident block dynamics, including a rank with no geometric Domains.

The RK driver now joins the retained owner tendency to the native forward and
inverse transform boundaries. Canonical scalar updates use the existing block
candidate. Canonical velocity aliases evaluate the original `base + h*gradient`
expression from the existing block checkpoint and retained owner gradient, to
preserve the original rounding order. They do not reconstruct that expression
from a rounded, recomposed velocity candidate.

Mass results are retained on their final block owners. Compiled metadata maps
canonical nodes to actual owner-program producers and requests only registered
alias results. Missing live producers remain errors. Temperature uses its
existing sparse owner result exchange. No full geometry-owner mass or gradient
return is needed by the production RK path.

Registered noncanonical aliases have an explicit timestep-base epoch and a
boundary-completed stage epoch. Fixed coarse scaffold and provisional low-level
scaffold have different lifetimes: the latter preserves the original untouched
`q1` initialization and previous provisional inverse values in a sparse keyed
registry. Initial adaptation runs outside this RK epoch. Final coarse velocity
restriction preserves the original evaluation order, followed by the native
signed boundary exchange. Final scaffold captures the completed boundary epoch.
Soil and surface levels retain their separate upstream lifecycle and correct
positions in writeback buffers; they are outside atmospheric RK arithmetic.

On the production route, `q1` retains only field metadata; its numerical arrays
are absent. The old native mass and velocity `tendency`/`rk` return arrays are
also absent. Interstage Domain RK arithmetic, temperature publication, boundary
updates, and provisional inverse publication are bypassed. The old RK closure
explicitly rejects accidental production fallback. The oracle route retains
independent Domain arithmetic and compares every original registered alias
bit-for-bit at each RK boundary, in addition to the existing transform checks.

The initial physical input boundary and final publication for physics,
adaptation, checkpointing and output remain explicit. Existing wavelet storage,
Domain topology, upstream physical fields and external lifecycle consumers are
not deleted. No Domain-shaped backup was introduced by the RK migration.

Qualification helpers are in `test/parallel_block_profile/rk_236/`. The final
review requires frozen source and binary identities, full exact 4/41-rank
adaptive/restart trajectories, checked native and independent oracle runs,
feature-off control, destructive old-primary poison without restoration, and
specific missing/stale/corrupt producer and forbidden-fallback rejections.
Local oversubscribed timings support no performance or scaling claim.
