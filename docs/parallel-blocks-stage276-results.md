# Stage276: fused primitive execution passes through 256 ranks

Reviewed 2026-10-09. Job **55493** completed all thirteen runs in **16:11**.
Independent review verified **3,304 sealed files**, reproduced the raw primitive
comparisons, science gates, placement and timing controls, and matched source
identity to the locally prepared v2 kit. Stage266 remains accepted; Stage273
remains separate and unpromoted. No production replacement is promoted here.

Archive: `rk276-review-55493.tar.gz`, SHA256
`4d4e0b8c7eeda53d56314c2cd15b867fff5de25bc10bbc5482a0e41236c9f30d`.
Machine-readable evidence is in `parallel-blocks-stage276-results.json`.

## Correctness and actual block execution

All thirteen runs pass science gates: twelve steps, two remaps and 62 VTK outputs
per run. Stage266-family exact checkpoint comparisons and supported main-relative
comparisons pass, including checked execution and retained 16/256-rank oracles.
Binary/checkpoint/VTK payloads are omitted from the review; recorded hashes and
numerical reports were verified, rather than rerunning those solvers locally.

Each of the six instrumented configurations checks the same **668,288 node
observations**, **10,692,608 inputs**, and **6,014,592 primitive outputs**. All
resident inputs, reference arithmetic and resident-derived outputs agree bitwise.
The input observations comprise 9,284,568 patch, 763,716 boundary and 644,324 ghost
values. This extends Stage275's owner-interior evidence to the actual sampled
producer stencil, using the existing refresh pathway.

Coverage is 46 common RHS/layer/level phases per rank: selected RHS ordinals
1,2,3,9,17,33,41,48 and layers 1/30 at visited mesh levels. Counts are repeated
observations, not unique mesh nodes. No mapped Domain nodes were skipped within
this scope. Unsampled phases/layers, other boundary outputs, final alias writes,
restriction and consumers remain outside this proof.

At 256 ranks, **all 96 ranks with no Domains execute the fused kernel**. They
compute 150,752 node observations and 1,356,768 primitive outputs, all exact.
Their share is **22.56%** of observed nodes versus 37.5% of ranks. This proves
useful arithmetic executes beyond Domain owners; it does not prove balanced or
efficient execution. Ownership comes from the actual Domain census, not rank IDs.

## Legacy gap and baseline scaling

These are uninstrumented main/Stage266 controls, with one run per configuration
and the established ordinary-step exclusions. They are not Stage276 replacement
timings. Main is the sealed main snapshot with established checkpoint I/O hardening.

| MPI ranks | Legacy main, ms/step | Stage266, ms/step | Stage266 gap |
|---:|---:|---:|---:|
| 83 | 543.124 | 755.113 | +211.989 ms / +39.03% |
| 160 | 368.380 | 548.911 | +180.531 ms / +49.01% |
| 256 | not run | 577.059 | +5.13% versus Stage266 at 160 |

The remembered 0.755-second result is reproduced at 83 ranks. The baseline still
slows down beyond 160 Domains: 160-to-256 speedup is 0.9512. These single-run
controls establish the current comparison, not uncertainty bounds or an exact
causal allocation of the performance gap.

## What the temporary timers support

For the optimized prototype runs, the sum across sampled phases of each phase's
maximum per-rank gather-plus-compute time is 2.914 ms at 83 ranks, 2.213 ms at
160, and 3.613 ms at 256. The corresponding reference-exchange aggregates are
139.222, 168.852 and 758.134 ms. Address-preparation aggregates are 2.279, 1.078
and 1.287 ms. Exact aggregates for all configurations are in the JSON report.

These are **uncalibrated diagnostic aggregates across the campaign's sampled
phases**, not milliseconds per production step. Maxima can come from different
ranks. The 256-rank run has oracles enabled; the 83/160 runs do not. Reference
exchange includes packing, allocations, collective synchronization and waiting,
so it is not a network-only measurement. No kernel scaling or recoverable-time
claim follows from these numbers.

They reinforce the design requirement to eliminate repeated reference packets
from production. Do not turn this validation Alltoall exchange into the new
producer transport or attempt to rescue it with shared-memory windows. The
direct arithmetic is correct; the complete production data path is untested.

## Recommended next substantial change

Proceed to one bounded, complete fused producer/restriction/consumer replacement
on accepted Stage266. Keep numerical fusion and direct compact addresses. Cache
geometry and addresses at their valid epochs; consume resident state/mean and
velocity through an explicitly validated refresh contract. Retain a bounded,
accounted physics adapter where physics remains Domain-owned.

Connect authoritative outputs to their real restriction and local consumers,
including required halo/alias publication and the original final shared-vorticity
write semantics. Retire the corresponding Domain computation and superseded
transfers in the same production path. Keep reference packets and duplicate
arithmetic exclusively in validation. Any retained Domain consumers and adapters
must be explicit; do not remove shared producer outputs they still need.

Validate all production calls/layers locally with optimized, checked and oracle
runs, then use the short cluster campaign: science/ownership gates through 256,
single uninstrumented main/Stage266/replacement controls at 83/160 and
Stage266/replacement at 256. Include preparation, refresh, physics transfer,
publication and restriction in total-step timing. Repeat a marginal or unexpected
timing result only as needed; resolve a regression before stacking another removal.

No reliable speedup percentage can yet be estimated. Stage270's roughly 194 ms
shared-producer region is not wholly removed by this primitive subset, and new
publication/consumer costs remain unmeasured. Efficient execution at >=640 ranks
still requires finer balance and multi-node qualification. Explicit MPI shared
memory remains a separate experiment after this production path is measured.

## Topology and campaign runtime

The archived allocation confirms bb13 has two AMD EPYC 9754 sockets, 128 cores
per socket, one online thread per core, two NUMA nodes (CPUs 0–127 and 128–255),
and 32 L3 cache groups of eight cores. Rank binding was independently verified.
This does not measure actual NUMA page residency or establish bb01–bb12 topology.

Slurm reports 2026-10-09 09:50:53–10:07:04, **00:16:11**, for four serial builds,
eight validations and five timing controls. Runner elapsed is 971.085 seconds.
The exact accounting and schedule are recorded in the runtime ledger. Queue and
collection/transfer are excluded. A comparable schedule can be budgeted at about
16–22 minutes; new production or multi-node work may take longer.
