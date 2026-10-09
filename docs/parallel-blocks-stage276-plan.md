# Stage276 — executable fused primitive prototype

Prepared against accepted Stage266 `ff740744608731a05d53771c32a3055309f3a586`
and legacy main `6a245db7f1632b0205c2bebe95e5f0cd535c2582`. The separate,
unpromoted Stage273 partition candidate stays out of the sealed sources.

## Why this change

Stage275 proved sampled scalar-owner interior residency, including 96,618
matches on 96 zero-Domain ranks at 256 ranks. It did not prove halo residency or
execute the upstream producer there. Flux/PV account for approximately 53% of
the remaining minimal Domain-fed input budget. Stage270 identified the shared
producer as the largest individual qualified candidate region (193.987 ms self,
195.821 ms inclusive); these figures do not establish recoverable time.

Stage276 executes part of that producer on actual block owners. It preserves
fused arithmetic and checks the prerequisites for replacing a complete slice,
rather than introducing the previous dependency-depth execution machinery.

## Bounded implementation

The isolated source overlay changes `multi_level.f90` and `ops.f90`:

- A fused kernel computes three mass fluxes (including the retained mass-physics
  residual), three edge PV values, kinetic energy and two centre-triangle
  vorticities. Arithmetic order follows the original `step1`.
- Mass state plus mean and velocity come from actual block storage. Seven compact
  node addresses per output node are cached at the native velocity-plan
  generation, covering patch, boundary and ghost records. Direct neighbour
  records avoid a dependency on later restriction-field metadata.
- Reference geometry, physics residuals, original operands and results travel
  in a bulk diagnostic exchange. These are validation packets, not a production
  communication design. No reference numerical input silently replaces a stale
  resident input.
- Optional `step1` output records capture results at their original assignment.
  This matters because later stencil visits overwrite shared vorticity entries
  with a different floating-point addition order. Comparing the initial fused
  output with the final overwritten entry is not a valid bitwise arithmetic
  test. A production replacement must preserve the authoritative final writes.
- The prototype stores results in private local buffers. It does not change
  tendencies or retire any production work. The exact checkpoint/VTK science
  gates and the existing dynamics/adaptation/remap oracles remain in force.

Scope is positive interiors of represented block patches, selected common RHS
ordinals 1,2,3,9,17,33,41,48 and layers 1/30, at every visited mesh level. Missing
Domain scaffold work is counted separately. Final boundary/alias publication,
flux restriction, B/Exner, other scalar fluxes, physics and production consumers
remain outside this replacement. No full-step speedup is claimed.

## Two goals, two kinds of evidence

The bbserv campaign uses four serial builds and thirteen runs on exclusive
nk/bb13:

| Runs | Purpose |
| --- | --- |
| 16-rank main and Stage266 | Established science controls |
| 16-rank optimized, checked and oracle prototype | Exact arithmetic, resident halo checks, science and oracle checks |
| 83/160-rank optimized prototype; 256-rank oracle prototype | Same input/output checks at different ownerships; actual fused work on ranks without Domains |
| Single main/Stage266 timings at 83 and 160 | Same-node legacy gap at supported rank counts |
| Single Stage266 timing at 256 | Intermediate baseline scaling from 160 to 256 |

There are no repeated identical timing configurations. The prototype reports
map preparation, reference capture/exchange, address preparation, resident
input gathering, fused computation/storage and validation separately. Diagnostic
phase timers are exploratory, not calibrated full-step performance measurements.
Geometry/reference traffic and Domain physics remain explicit costs/limitations.

The result is an executable and checked primitive implementation plus evidence
about input freshness and kernel distribution. It is **not yet** an accepted
Domain-path removal. Useful work on zero-Domain ranks is necessary but cannot
establish efficient scaling. J5 has 160 Domains; >=640 ranks remains a later
multinode qualification, with finer splitting/load balance still relevant.

## Decision after review

1. Require exact reference arithmetic and exact science. A failed resident-input
   hypothesis is reported by storage class and blocks promotion.
2. Use address/gather/compute costs and zero-Domain-rank output counts to decide
   whether this implementation is a viable upstream producer. Counts are not
   time-weighted balance or parallel efficiency.
3. Connect one complete producer/restriction/consumer slice, with cached geometry
   and bounded halo publication. Retire the displaced Domain computation and
   exports in that production change, preserving shared vorticity write order.
4. Gate that uninstrumented replacement immediately against Stage266 and main
   at P<=D and against Stage266 above D. Stop on a repeatable regression; do not
   accumulate more removals before resolving it.

Estimated cluster execution is 18–25 minutes after allocation, with uncertainty
for new prototype costs. Stage275's four builds/eight validations took 14:35;
Stage276 adds five timing controls. Queue and review transfer are extra. The
collector records Slurm start/end/elapsed automatically.

Reproducible source overlay, preparation and one-command launch/status/collection:
`test/parallel_block_profile/primitive_276/`. Delivery details and actual local
validation are recorded in `parallel-blocks-stage276-preparation.json`.

## Packaging correction: v2

The initial launch failed in preflight, before `sbatch`. Renaming Stage275 to
Stage276 by unrestricted text substitution had also altered `275` inside the
expected grid SHA256. The packaged grid retained the correct established bytes.
V2 restores the checksum and verifies it against the actual packaged grid during
preparation, sealing and preflight. Tests cover a corrupted expected checksum and
an altered external fixture. All three fixture hashes also match the successful
Stage275 campaign. Solver archive, source manifest, campaign and local MPI
validation are unchanged; the repair does not require new solver runs locally.
Use the v2 commands in the campaign README.
