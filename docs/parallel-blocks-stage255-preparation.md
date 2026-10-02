# Stage255: stable owner depth index prepared for bb13

The Stage254 partial review identified positive-depth reads and producer
traversal as the dominant measured owner costs. Stage255 tests whether
removing repeated scans of all owned operations helps directly and reduces
MPI waiting. The one-file candidate starts from the exact Stage251 c5
`parallel_block_program.f90` (SHA-256 `eff68d917b06eebea30aa81f12ec7bb1c8163de3c502e27b087025da4a66a8e8`).
At compilation, it builds stable buckets of operation indices for each
dependency depth. Execution visits the same operations in the same order
within every depth, but does not examine operations of other depths. The
output source SHA-256 is `a65f465f483e789f68f0f454bb91088fe81309d4cf51d62950a0d5bdd299f2dc`.
No Stage253 change is included.

The candidate built locally with `make -j1` and the c5 feature selectors.
The four-rank J4 restart replay passed all 12 printed step states, two
remaps, checkpoint SHA-256
`345161aac6fc2a89a0241de04a56ea33847dd3293801e8f05734fbd61059f95c`,
and byte-identical contents of both 31-member VTK archives against the
Stage254 local baseline. The local machine has only four usable cores, so
these tests establish functional behavior, not cluster performance.

The sealed bbserv kit reuses the proven Stage253 A-B-B-A runner structure,
checks Stage251 job 55256 inputs, baseline binary, physical-core map and
the 23 shared physics objects, builds only the changed c5 source, and runs
eight full 83-rank J5 trajectories. It includes raw run logs and per-rank
placement records in the review archive. The kit remains a trial until its
bb13 evidence is independently reviewed. Instructions are in
`test/parallel_block_profile/owner_depth_255/README.md`.
The v1 kit archive SHA-256 is
`74769893f6e0d2a3aeedd27534a9b8ad652d425470ecc3fbd187a9cd42464ab8`.
