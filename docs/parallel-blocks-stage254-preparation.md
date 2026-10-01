# Stage254: owner-program timing probe prepared

Stage253 job 55259 showed that eliminating two empty collectives is a
small benefit and leaves the large velocity-owner regression intact.
Stage254 measures where time remains inside the owner source path.
Starting from the exact sealed Stage251 c5 source, the probe inserts
timers into `parallel_block_program.f90` and a report call at the
existing parallel-block profile window in `parallel_block_mpi.f90`.
It does not change the numerical operation order.

The probe reports mean-rank time for complete owner execution,
resident fetch, internal reads, MPI within those reads, producer
evaluation/scans, final delivery reads and MPI within delivery. It
also reports read time by dependency depth. Complete execution
includes fetch, internal reads and producer loops; MPI time is a
subset of read time, so these categories are not additive. Reports
occur at the existing eight-step and four-step profile windows, with
no new collective communication in the hot path.

The J5 probe build passed locally with `make -j1`. The four-rank J4
probe build then completed a full restart replay with twelve exact
step states, two remaps, 62 byte-identical VTK members and the same
final checkpoint bytes as the prior J4 baseline (SHA-256
`345161aac6fc2a89a0241de04a56ea33847dd3293801e8f05734fbd61059f95c`).
Both owner timing windows were present. The local machine has four
usable cores; its timing values are not used as cluster evidence.

The sealed v1 kit SHA-256 is
`fd3d988679989c8cc2cd6125e5f52552ccc3decdd8732168a21415e24779bb21`.
It verifies the Stage251 result, source, baseline binary, input and
physical-core map; pins all 23 external physics objects; copies the
source; applies and hashes the two-file timing probe; and builds c5
serially. It runs unprofiled and profiled A-B-B-A batches at 83 ranks
on bb13, comparing the exact Stage251 c5 binary to the instrumented
c5 binary. This calibrates instrumentation cost while preserving the
full J5 numerical and placement gates. The archive unpack/hash,
Python and shell syntax, exact Stage251 state-file launcher layout,
local source transformation and packaged preflight path were checked.
Full cluster timing still requires the bbserv run.

No source promotion or further Domain removal follows automatically.
The no-edit commands are in
`test/parallel_block_profile/owner_breakdown_254/README.md`.
