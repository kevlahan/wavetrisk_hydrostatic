# Stage256: direct depth-indexed c4/c5 comparison prepared

Stage255 reduced c5 ordinary-step time by 70.408 ms, but its native
velocity-source phase fell by only 27.963 ms per all-step. The one-file
depth index applies to every owner program and therefore also benefits
c4. It would be misleading to count the entire c5 speedup as closure of
the c4-to-c5 velocity-owner gap.

Stage256 builds c4 and c5 from the same exact Stage251 source with the
same depth-indexed `parallel_block_program.f90` (SHA-256
`a65f465f483e789f68f0f454bb91088fe81309d4cf51d62950a0d5bdd299f2dc`).
The builds differ only in the original Stage251 feature-ladder selector
`BLOCK_VELOCITY_OWNER`. The runner checks the 23 compiled physics objects,
build flags, input and physical-core map, then runs unprofiled and
profiled A-B-B-A pairs at 83 ranks on bb13. Exact J5 state, checkpoint,
VTK and placement gates are required on every run. The resulting
paired difference is the residual velocity-owner feature cost after the
common depth optimization; it is not yet a comparison with `main`.

Both feature variants built locally with `make -j1`. Four-rank J4 restart
replays of c4 and c5 matched all 12 printed states and two remaps, the
exact checkpoint SHA-256
`345161aac6fc2a89a0241de04a56ea33847dd3293801e8f05734fbd61059f95c`,
and the contents of both 31-member VTK archives. The local machine has
four usable cores, so these are numerical checks only. The sealed nine-file
v1 kit unpack-verified with SHA-256
`8732af0fb2f6c2770e7251c7178ae564f41d23d71a13c951fb17dc365e276007`.

No source is promoted or further Domain removal undertaken while the
remaining regression is being measured. The no-edit transfer and launch
commands are in `test/parallel_block_profile/depth_ladder_256/README.md`.
