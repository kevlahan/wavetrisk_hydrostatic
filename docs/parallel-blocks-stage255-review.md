# Stage255: depth indexing helps, but most velocity-owner cost remains

Job 55261 completed all eight full J5 runs on bb13. The independent
reviewer verified all 1,350 archive files, the exact Stage251 c5 baseline
binary, the one-file depth-index candidate source and binary, unchanged
23 external physics objects, serial build selectors, 83 physical-core
placements, twelve raw step states and two remaps in every run. The cluster
runner recorded the established exact checkpoint hash and 62 VTK members
for all eight runs. Checkpoint and VTK payloads are absent from the review
archive, so their bytes are supported by the runner's recorded gates rather
than a second byte comparison after transfer.

In the unprofiled A-B-B-A batch, the baseline ordinary-step mean was
1.463083 s and depth-indexed c5 was 1.392675 s: a 0.070408 s (4.81%)
improvement. Both paired differences favored depth indexing (-0.065254
and -0.075561 s). The profiled batch likewise improved by 0.071333 s,
with both pairs favorable (-0.075368 and -0.067297 s).

Across twelve profiled steps, mean-rank complete-timestep time fell by
0.061042 s/step. Native velocity-source time fell by only 0.027963
s/step; native velocity-gradient time fell by 0.002563 s/step. The
remainder of the full-step saving lies in other phases. This is expected:
the one-file index applies to every owner program, including programs
already active in c4. Therefore the 70 ms c5 improvement cannot all be
credited against the c4-to-c5 velocity-owner regression. A matched
comparison of c4 and c5 with the same depth index is needed to measure
the residual feature gap.

Depth indexing is a sound performance improvement in this campaign but
is not yet promoted. The Stage253 empty-collective optimization was not
included. Further Domain removal remains paused while the remaining
velocity-owner penalty is measured and addressed.

Independent reproduction:

```bash
python3 test/parallel_block_profile/review_cluster_255.py \
  /Users/kevlahan/Desktop/rk255-review-55261.tar.gz
```
