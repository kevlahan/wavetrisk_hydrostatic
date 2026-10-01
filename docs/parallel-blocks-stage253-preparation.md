# Stage253: owner-program empty-collective trial prepared

Stage252 job 55258 reproduced the c4-to-c5 ordinary-step gap at
+0.405565 s (83 ranks) and attributed 89.2% of the profiled all-step
gap to native velocity-source execution. The first experiment changes
only `src/parallel_block_program.f90`: a program with no resident
producer skips the empty resident fetch, and all programs skip the
empty depth-zero read. The consumer-epoch and layer guards remain in
place. The compile step makes the resident decision collectively so all
ranks take the same execution path. Other owner communications remain.

The modified file builds serially with the J4 and J5 c5-related build
selectors on the local four-core machine. A full four-rank J4 restart
replay compared the previous baseline and optimized builds over 12
steps, two remaps and 62 VTK members. Step states, VTK bytes and final
checkpoint bytes were exact; the exact numerical checkpoint gate passed.
The common checkpoint SHA-256 was
`345161aac6fc2a89a0241de04a56ea33847dd3293801e8f05734fbd61059f95c`.
No local timing is a scaling result.

The sealed Stage253 v2 kit SHA-256 is
`8f4b453b3baf07042897abccaec54dc59a96361b498a3f443fe05bb782700648`.
It pins Stage251 job 55256, verifies all sealed candidate source files,
replaces only the program source in a private copy, builds c5 with
`make -j1`, and compares the optimized binary against the exact Stage251
c5 binary. Its eight J5 runs use unprofiled and profiled A-B-B-A batches
with state/checkpoint/VTK/physical-core gates. The kit has passed Python
and shell syntax, source-manifest checks through the local archived
Stage251 evidence, and tar unpack/hash verification. The local Stage251
review archive omits binaries, so full preflight can only finish on
bbserv. No cluster performance result or promotion exists yet.

The kit and no-edit bbserv commands are in
`test/parallel_block_profile/owner_collectives_253/README.md`.
