# Stage254: owner-program timing breakdown (partial campaign)

The Stage254 bb13 job 55260 archive passed independent SHA-256 review of
all 1,013 archived files. Five of eight planned 83-rank J5 runs completed,
and the sixth solver run completed successfully but its Python postprocessor
raised `Incomplete owner timing windows`. The parser expected all native
velocity-source calls (2,880 and 1,320 per rank in the two windows), whereas
the owner timer counts only calls with the external-owner callback (2,160
and 990). The timing report itself has both complete windows and every
seven measured categories. The parser expectation is corrected in source.
No numerical solver or Slurm failure occurred.

The four unprofiled A-B-B-A runs passed the exact state and checkpoint
hashes, 12 steps, two remaps, 62 VTK members and 83 physical-core placement
gates. Their ordinary-step means were 1.444243 s for the sealed Stage251
c5 binary and 1.460236 s for the same source with the timing probe compiled
in and disabled: a +0.015993 s probe-build difference. The two paired
differences were -0.004118 and +0.036105 s, so the observed 16 ms is an
imprecise calibration, not a correction to subtract from future timings.
The first profiled baseline and instrumented solver runs preserved the same
raw 12-step state trace and placement. Their mean-rank native velocity-source
totals were 4.2547 s and 4.2332 s across 12 steps; this one pair shows no
large active-probe cost. The instrumented run's output gates passed inside
the runner before its parser exception, but its result row lacks the saved
checkpoint and VTK evidence. The archive omits checkpoint and VTK payloads,
so those bytes cannot be independently rechecked after transfer.

In the usable instrumented profile, mean-rank owner execution totaled
3.46941 s across the 12 steps (0.28912 s/step). Internal dependency reads
used 1.982633 s (0.16522 s/step), including 1.822852 s in MPI_Alltoallv
(0.15190 s/step). Producer loops used 1.331742 s (0.11098 s/step).
Final delivery reads used another 0.189441 s, almost all of it in MPI.
These are inclusive timers: MPI is contained in read time, and delivery
lies outside owner execution. The positive dependency depths 1–4 account
for 1.977031 s of the 1.982633 s internal-read total; depth zero is only
0.005598 s. Read time is highest at depths 2 and 4. The large MPI times
may include waiting for slow producers, so this profile alone cannot
attribute all of that time to communication overhead.

This establishes that the two empty collectives addressed in Stage253
were never the main cost. A depth-indexed producer traversal is a narrow
next candidate: it avoids rescanning every owned operation at each depth
while preserving operation order within each depth. A matched 83-rank
comparison is needed to measure whether it also reduces MPI waiting.
Further Domain removal and promotion remain paused during this investigation.

Independent reproduction:

```bash
python3 test/parallel_block_profile/review_cluster_254.py \
  /Users/kevlahan/Desktop/rk254-review-55260.tar.gz
```
