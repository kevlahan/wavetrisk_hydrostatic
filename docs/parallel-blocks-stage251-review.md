# Stage251: the velocity-owner transition is the largest J5 regression

The independent reviewer verified all 2,705 files in job 55256's return
archive, eight serial bb13 builds, 83 distinct physical-core placements,
sixteen full J5 runs and fifteen passing main-relative checkpoint numerical
gates. Every run completed twelve steps, two remaps and 62 VTK members.
All candidates had the established state and checkpoint hashes. These
results support timing comparisons across the dependency ladder; they do
not change the numerical acceptance policy or promote any source.

The single-run 83-rank screens rose from 0.539078 s/ordinary step for
read-only `main` to 0.771932 s with the common Stage202/203 block pathway
(`c0`), then to 1.645313 s after all six early owner features (`c6`).
The largest adjacent screen was `c4` to `c5`, which enables
`BLOCK_VELOCITY_OWNER`: +0.416384 s. Two independent A-B-B-A batches
confirmed that transition with four full runs per variant. `c4` averaged
1.062451 s and `c5` 1.448877 s, a +0.386426 s/step (+36.37%) increase.
All four adjacent paired differences were positive, from +0.370923 to
+0.396205 s. This one transition accounts for roughly 35% of the
single-run `c6`-to-main gap, though the denominator and the other
adjacent screens have not received the same confirmation.

The remaining screens also reveal substantial costs at the common
Stage202/203 pathway (+0.232854 s to main), primitives (+0.105882 s),
diffusion (+0.143140 s) and initial tendency (+0.164833 s). Those are
search leads from one full run per build, not confirmed feature effects.
The small column and mass-owner differences are likewise only screens.

`BLOCK_VELOCITY_OWNER` replaces the direct native velocity-source
calculation with an owner program. Its plans have dependency depths 2 or
4; each execution fetches resident inputs and performs a collective read
at every depth, then delivers source values to the gradient path. The
Stage251 logs show 35 new owner-program plans per rank, but contain no
phase timing. Communication latency is a plausible mechanism, not yet
measured attribution. The next experiment should profile `c4` and `c5`
with the same source and isolate velocity-source execution from gradient
delivery before changing arithmetic or ownership. Domain removal and
promotion stay paused.
