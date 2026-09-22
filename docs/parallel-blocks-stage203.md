# Stage203: one inverse workspace and patch-wide scalar fields

Production default following [Stage229 requalification](parallel-blocks-stage229-review.md):
`STAGE202=on STAGE203=on FP_CONTRACT=off`. Job 54737 confirms a 45.715 ms (5.666%)
ordinary-step reduction at 83 tasks and passes checked/oracle correctness gates,
including 41 tasks/40 Domains. Always use fresh object and binary directories
when changing build settings; make does not track changes to compiler flags.

`STAGE203=off|fields|inverse|on` retains the reference and separate feature modes.
Use `STAGE202=on STAGE203=off` for the Stage226 comparison, or explicitly set
both off for the original path. Every non-off Stage203 mode requires Stage202 on.
The inverse workspace is incompatible with `WAVETRISK_STAGE192_SPARSE_OUTER=1`;
leave that experimental runtime option unset or 0, or build with `STAGE203=off`
(or `fields`) when investigating sparse outer installation. Existing guards remain.

The original Stage203 campaign, bb13 job 54695, passed correctness but missed its
former 100 ms threshold; see [that historical review](parallel-blocks-cluster-203-54695.md).
Its frozen comparison used Stage202-on commit `7c406f4d`. The implementation and
historical protocol below are retained; the later Stage229 acceptance supersedes
the original experimental status and threshold.

## Inverse execution

Scalar prolongation/reconstruction and inner-vector interpolation now execute
alongside outer-vector reconstruction in the existing inverse workspace. Integer
stencils are compiled at each topology generation. Geometry/masks are read at
execution, and arithmetic order and contraction-off policy are retained.

The inverse-loop gathers/scatters become local state transitions. Production
omits their intermediate block ghost refreshes. Required seed transfers, signed
alias communication and final publication remain. This changes the scalar/inner
compute owner from the block owner to the existing workspace owner; the campaign
records phase timings and per-rank work to expose any resulting imbalance.

Alias-overwritten interiors have a sparse canonical-value cache. Kernel writes
update that cache; phase transitions restore it and fixed scaffold values.
Outer installs update canonical values using the same conservative selected set
as Stage202, including alias effects. Final installation saves/restores just the
alias-overlap subset while exporting canonical interiors and fresh boundaries.
Scratch is retained per topology generation, not allocated per production install.
Consumer-visible boundary epochs remain separate from kernel dependency epochs.

In oracle mode the original block kernels still run independently. Before each
old gather replaces workspace results, the candidate's restored canonical values
are compared bitwise against those block results. Existing signed-alias and
Domain comparisons also remain. The production branch avoids that extra work.

## Scalar execution

New field batch APIs translate field slots once per patch span, preserving wire
order. They support compact production and full oracle storage. They are used for:

- compact received physics installation;
- dynamic restriction ghost packing/installation;
- direct same-rank patch-to-ghost field copies, bypassing wire scratch;
- selected field portions of full production ghost transport;
- direct-flux patch inputs/outputs, including the mass shortcut and separate
  shared horizontal geometry reads.

Restriction/divergence arithmetic retains its existing stencils; this stage
changes the surrounding field access and transport interfaces. Full geometry
transport stays deduplicated. General record APIs remain available. Poisoning,
coverage checks, physical/inactive classes and signed-zero bits are preserved.
Local full-rebuild ghosts still use the old buffer path; the repeated dynamic
refresh path takes direct copies. There is no claim that all representation
maintenance or all compatibility storage has been eliminated.

## Pass ledger and limits

| Operation | Stage202 | Stage203 production |
|---|---|---|
| Inverse-loop scalar/vector gather | transfer per synchronization | sparse canonical/scaffold restoration |
| Inverse-loop scatter/outer install | transfer per synchronization | canonical-cache update when required |
| Inverse-loop block ghost refresh | 4L+1 calls for L inverse levels | omitted; kernels use workspace aliases |
| Initial inverse seed | workspace import and boundary setup | retained, with additional stencil coverage |
| End-of-transform block install | incremental results already installed | complete canonical install at consumer boundary |
| Provisional final-boundary contract | final boundary refresh | retained; second install for this contract |
| Dynamic local scalar ghost | field pack then field install | direct selected field copy |
| Remote scalar ghost | per-sample record adapters | patch-wide selected field packing/install |

Thus this is not zero-copy inverse execution. Entry/exit transfers, sparse
canonical snapshots and extra stencil metadata must be included in the outcome.
A J5/J7 inverse transaction removes nine intermediate block ghost refresh calls;
this is a source count, not nine equal measured costs or a timestep timing claim.

Only observer builds define `STAGE203_LEDGER`. Their per-rank inverse records
report replaced gather/scatter calls, scalar outputs, inner cells, restored
canonical values and final installed interior values. Counters are whole-run,
including checkpoint steps. Bare builds omit counter updates and printing.
The report validates ledger rank/transaction coverage. This ledger verifies
mechanism; it does not convert copied bytes or MPI waiting into recoverable time.

## Tests and decision

Local J4/J6 four-rank, 12-step optimized, runtime-checked, independent-oracle and
isolated observer runs all pass exact Stage202 restart/workload comparisons and
post-restart remapping. The fixed legacy numerical screen also passes. See
`parallel-blocks-stage203-local-results.json` for completed evidence and exclusions.
Actual-module unit tests check batched fields against independent scalar accessors
at O0/O2, including compact/full storage, geometry classes, signed zero, untouched
fields and empty ranges. An actual-helper test covers overlapping canonical and
alias slots, outer canonical updates and preservation of distinct installation
epochs. Campaign tests preserve reversed pairing, fresh legacy
phase comparison, and rejection of marginal gains. Frozen-kit tests execute the
real patch/hash/submission sequence with compilation replaced by a stand-in.

Cluster: seven builds, 18 solver runs; fixed 83 ranks on exclusive bb13, 42/41
across sockets, identical physics objects, contraction off. Two checked modes,
three warmups, two reversed bare triples, three coarse observers (including a
fresh legacy control), and two reversed block core-PMU pairs. No memory-event,
stack-sampling or rank-scaling passes.

All numerical, adaptive-workload, restart/remap, placement and observer-quality
gates must pass. Require **at least 0.10 s/ordinary step in BOTH reversed pairs**
against Stage202. Report fresh coarse residuals versus legacy separately from
bare medians. A small gain is not a successful Stage203 outcome. Retain the
single review tarball workflow.
