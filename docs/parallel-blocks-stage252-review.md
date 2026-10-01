# Stage252: the velocity-owner cost is in source execution

Job 55258 completed all eight full 83-rank J5 runs. The independent
reviewer verified all 1,349 archive files, the exact Stage251 `c4` and
`c5` binaries, 83 physical-core placements, twelve steps, two remaps and
62 VTK members per run. Every run preserved the established exact state
and checkpoint hashes. No new source was built in this campaign.

The unprofiled A-B-B-A reproduced the Stage251 velocity-owner jump:
`c4` averaged 1.060839 s/ordinary step and `c5` averaged 1.466404 s,
a +0.405565 s gap. Both paired differences were positive (+0.400929
and +0.410201 s). With profiling enabled, the gap was +0.378467 s;
its two paired differences were +0.378575 and +0.378359 s. The profiler
changes absolute timings somewhat, but the large regression remains.

Across the profiler's twelve-step windows, the mean-rank complete-timestep
gap was 0.384908 s/step. The native velocity-source phase accounts for
0.343379 s/step, or 89.2% of that gap; the native velocity-gradient phase
accounts for 0.009434 s/step, or 2.5%. The remaining 8.3% lies outside
those two phases. These phase totals include startup/remap/output steps,
whereas the ordinary-step figures exclude them; the phase attribution is
an all-step measurement, not a direct decomposition of the ordinary-step
median. Both profile windows gave similar source-phase differences.

Code inspection identifies two unconditional collective calls in every
owner-program execution even when the program has no resident producers
or a dependency-depth-zero read has no inputs. The velocity-owner program
uses external primitive/diffusion values and has no resident producers;
its compiled plans have dependency depths two or four. Skipping only these
provably empty collectives is the first narrow optimization to test.
Other communication at positive depths remains necessary and may account
for much of the source-phase cost. Promotion and further Domain removal
remain paused until the optimized executable passes full numerical gates
and a matched bb13 comparison.
