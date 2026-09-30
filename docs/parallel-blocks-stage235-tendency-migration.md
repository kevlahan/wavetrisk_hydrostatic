# Stage 235 — complete climate tendency kernel migration

Implementation is in `codex/performance-investigation`. `main` is untouched.
Local qualification passed. The [source-bound results report](parallel-blocks-stage235-local-results.json)
binds the final source, four build variants, executable hashes and validation logs.
This stage makes no performance claim and does not enable the new route by default.

## What moved

The new distributed operation programs execute the remaining tendency arithmetic
on final block owners: scalar/velocity diffusion, scalar physics and advective
flux, mass flux restriction and divergence, rotational diffusion, ordered velocity
source restriction, and the final velocity source-to-gradient transition. They
join the previously migrated hydrostatic columns, regular/southwest/pentagon
primitives, temperature restriction/divergence, and retained RK tendency assembly.

The route now includes **stage one of every RK step**. Stage one binds a fresh
physical-alias epoch after the upstream boundary update. Canonical interiors are
read from resident `Block_Data`; the bridge imports only registered aliases and
scaffold from the upstream restart/adaptation/physics field. Interstages consume
the completed provisional inverse epoch. The saved final-inverse consumer epoch
is distinct and is never substituted for either input epoch.

| Producer | Numerical consumer | Current route |
| --- | --- | --- |
| Resident mass, temperature, velocity | Columns, primitives, diffusion | Deduplicated owner operand requests |
| Owner pressure/geopotential | Regular primitive | Owner column result |
| Owner primitive/PV/vorticity | Pentagon correction and velocity source | Explicit producer identities; compact result routes |
| Owner mass flux and diffusion | Velocity physics and source | Direct reads of current owner program results |
| Restricted owner Bernoulli/Exner and source | Final gradient | Direct producer requests; no Domain scratch operands |
| Final velocity gradient | Retained block assembler | Direct retained patch result |
| Scalar physics/flux | Existing block temperature restriction | Explicit compact scalar transport |
| Mass/velocity tendency | Remaining RK closure | Explicit native workspace return |

Geometry compilation preserves the original sequence of remote snapshots,
ordered local boundary copies, remote installs, selected-edge restriction, and
last-writer precedence. Compilation is reused until the grid generation changes;
numerical operands and output readiness are refreshed for each physical layer.

Missing producer definitions remain missing. Only slots proven outside **all**
write destinations are assigned the original scratch-initialization value zero.
A live request for an absent producer fails. Unused outputs may remain undefined;
they are not repaired by reading a Domain field. Independent check mode retains
the original arithmetic and compares intermediate results bit for bit.

In the complete production route, per-layer compatibility writes of pressure/geopotential,
PV, vorticity, kinetic energy, Bernoulli, Exner, and diffusion products are removed
from the tendency kernels. Surface initialization still exports its existing
compatibility values before the layer loop. Source values are no longer returned to native geometry
storage before the gradient. The mass/velocity final tendency return is retained
because the remaining RK closure still consumes it.
The RK return also spans untouched scaffold slots: these preserve the original
`zero_float` initialization (+0). They are separate from numerical graph inputs.
The complete runtime oracle exposed an earlier return path that poisoned those
unused slots; the corrected return preserves initialization without introducing
zero substitution for a missing live producer.


## Deliberate boundaries

This completes the locally validated **climate tendency kernel execution** migration. It does not delete every Domain object or all Domain storage.
Topology and geometry remain Domain-described. Stage-one aliases, alias reference
means and surface geometry still have explicit upstream inputs. Compact scalar
transport and RK closure returns remain. Restart, adaptation, post-step physics,
wavelet scaffold, diagnostics and output have separate storage lifetimes.

The initial input bridge does not recreate a full Domain-shaped solution snapshot.
Deleting the remaining upstream fields requires migrating those producers and
consumers first. Likewise, removing the RK return requires moving its closure
updates, not merely deleting the return while the old RK code still reads it.

The migrated diffusion physics currently selects the climate callback and the
primitive route requires compressible dynamics. The existing feature-off path is
retained. This is a correctness and ownership migration; local oversubscribed
MPI elapsed times are not evidence of a speedup or scalability improvement.
The generic distributed-operation executor remains experimental. An initial
four-rank full run hit the 2400-second watchdog during its final step. That run is
preserved as failed qualification evidence. A bounded CPU sample then found
2108 of 3596 main-thread samples in producer lookup on one rank. The original
hash translated contiguous producer IDs into overlapping ranges and long
linear-probing clusters. Mixing the identity bits before table reduction removes
that pathology without changing operation order or numerical arithmetic; all
acceptance tests use the corrected source. Products in the mixer fit signed
64-bit integers without relying on integer overflow.

The corrected local runs are still slower than the stored Stage230 reference.
These are not paired performance measurements and do not establish a cluster
speedup. Correct ownership is not sufficient for production promotion. The next
bounded performance work should examine remaining dispatch and plan/transport
costs, then batch or simplify execution while retaining the exact numerical and
no-fallback tests. A new broad cluster measurement campaign is not justified by
this correctness result alone.


## Build and qualification

Use a fresh build directory for each selector combination. All of the following
must be enabled for the complete route:

```
TEST_CASE=climate PARAM=param_J4 FP_CONTRACT=off STAGE202=on STAGE203=on
BLOCK_COLUMNS=on BLOCK_PRIMITIVES=on BLOCK_DIFFUSION=on
BLOCK_MASS_OWNER=on BLOCK_VELOCITY_OWNER=on BLOCK_TENDENCY_INITIAL=on
```

Use `check` for the first five owner selectors to compare intermediate numerical
results. `DEBUG=check` separately enables bounds, initialization and floating-point
checks. Runtime dynamics/adaptation/remap oracles remain independent switches.

The reproducible helpers are in `test/parallel_block_profile/tendency_235/`:

- `build_tendency.py`: frozen source and oracle/native/production/control builds.
- `run_tendency.py`: all-stage execution and exact adaptive/restart validation.
- `reject_invalid.py`: missing/wrong-layer producer rejection, corrupted source
  detection, and a successful run with removed scratch kept poisoned. It also
  poisons the compatibility primary fields after explicit operand binding, then
  restores them only for the downstream RK boundary.
- `review_tendency.py`: binds results to source, executable and log hashes;
  requires exact 4/41-rank full trajectories, independent oracles, a feature-off
  control and injected tests before writing the durable report.

The 41-rank fixture has 40 Domains and 160 blocks. Acceptance requires actual
owner computation on a rank with zero Domains. Full comparisons include all 12
steps, two remaps, final checkpoint fields, step states and 62 VTK members.
Cluster qualification and performance measurement follow local correctness; no
bbserv job has been submitted by this stage.

## Final local results

- Fresh oracle, native checked, optimized production and feature-off builds passed.
- Four and 41 ranks: 12 steps, two remaps, exact checkpoint fields and step states,
  and all 62 VTK members matched their same-rank Stage230 references.
- Rank 40 had zero Domains and owned blocks; it executed all 12 initial tendency
  epochs and contributed 48 compiled owner programs over the adaptive trajectory.
- Independent intermediate arithmetic and runtime dynamics/adaptation/remap oracles
  passed through all four RK stages of the one-step checked run.
- Removed Domain scratch and primary compatibility mass/temperature/velocity
  fields remained poisoned during the layer kernels; the step state still matched.
- Injected missing mass producer, wrong-layer diffusion read and corrupted velocity
  source each failed at the intended guard/oracle. The feature-off control passed.
- Source/executable/log integrity checks and `git diff --check` passed. `main`
  remains `6a245db7f1632b0205c2bebe95e5f0cd535c2582`. Changes are uncommitted.

The timed-out and deliberately interrupted pre-fix runs are preserved separately
and do not contribute to these passes. See the [handoff](parallel-blocks-handoff-235.md)
for explicit remaining storage boundaries and the next bounded work.
