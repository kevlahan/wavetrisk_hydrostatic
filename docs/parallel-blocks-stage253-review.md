# Stage253: empty collectives are only a small part of the regression

Job 55259 completed the Stage253 v4 A-B-B-A comparison on bb13. The
independent reviewer verified all 1,350 archive files, the serial c5
build selectors, unchanged 23 compiled physics objects, both binary
hashes, eight full 83-rank J5 runs, physical-core placement, twelve
printed step states and two remaps per run. The cluster runner recorded
the established exact checkpoint hash and 62 VTK members in every run.
The review archive does not contain checkpoint or VTK payloads, so those
two output gates are supported by the runner's recorded hashes and
counts rather than independent byte comparison after transfer.

Unprofiled ordinary-step means were 1.455481 s for the verified
Stage251 c5 baseline and 1.440593 s for the one-file optimization:
a 0.014888 s improvement. The two paired differences were -0.027576
and -0.002200 s. Profiled ordinary-step means differed by -0.011217 s;
its pairs were +0.003342 and -0.025775 s. The small effect varies
between pairs, so this campaign does not establish a precise speedup.

Across all twelve profiled steps, the mean-rank native velocity-source
phase improved by 0.005854 s/step. The complete timestep improved by
0.015167 s/step, while gradient time was essentially unchanged. These
all-step phase totals and the ordinary-step medians use different step
sets. Compared with the +0.405565 s c4-to-c5 ordinary-step gap seen in
Stage252, the apparent 0.014888 s gain is only about 3.7%; subtracting
between separate jobs is illustrative, not a direct paired estimate of
the remaining gap.

The current owner program performs positive-depth reads and scans its
owned operations once for each dependency depth. The Stage253 baseline
log contains owner plans with maximum depths 2, 4, 6 and 8, and 4200
native velocity-source calls over twelve profiled steps. The next
measurement should divide that source phase into positive-depth
communication, producer evaluation/scans, external callbacks and
output delivery, with a matched uninstrumented run to calibrate timer
cost. That evidence is needed before another algorithmic change.

The empty-collective change remains an isolated trial on the
performance-investigation branch. It is not promoted to main; further
Domain removal remains paused while the much larger regression is
investigated.
