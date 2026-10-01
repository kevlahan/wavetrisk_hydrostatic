# Stage252: measure where velocity ownership spends the additional time

Stage251 confirmed that enabling `BLOCK_VELOCITY_OWNER` adds 0.386426 s
per ordinary 83-rank J5 step in two A-B-B-A batches. The Stage252 kit
reuses the exact `c4` and `c5` executables from job 55256. Its launcher
checks the SHA-256 of the verified Stage251 result, both binaries, the J5
input and the physical-core map before submitting. There is no rebuild or
source change.

On `nk`/bb13, the kit runs an unprofiled `c4-c5-c5-c4` batch followed by
the same batch with the existing parallel-block profiler enabled. Each
full run must again complete twelve steps, two remaps and 62 VTK members,
match the established state and checkpoint hashes, and keep all ranks on
the verified distinct physical cores. The runner records independent
ordinary-step medians, and the profiled runs collect the built-in
`native velocity source` and `native velocity gradient` phase reports.
The unprofiled batch checks whether the Stage251 timing gap reproduces in
the new allocation; its comparison with the profiled batch bounds the
profiler's impact. Profile phase times will localize the jump to source
execution or gradient delivery. They cannot by themselves separate
useful arithmetic from MPI wait inside a phase.

The sealed kit passed Python and shell syntax checks, archive unpack/hash
verification, profile-report parser checks, fake-Slurm submission and
wrong-prior-job rejection. Cluster results are pending. The
[kit README](../test/parallel_block_profile/velocity_profile_252/README.md)
contains no-edit transfer, submission and collection commands. Domain
removal and promotion remain paused until a measured bottleneck is fixed
and compared against the same-source baseline.
