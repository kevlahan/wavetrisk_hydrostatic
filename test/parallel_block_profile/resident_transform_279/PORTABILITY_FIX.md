# Stage279 v2 — GCC13 checked-build repair

Job55505 failed after4:32 during the checked build, before any solver run.
GCC13.2 diagnosed selected.offset/lbound/ubound as maybe-uninitialized in
compile_fixed_inputs at the allocatable self-reassignment. The local GCC15
qualification did not flag this expression. Optimized cluster builds succeeded.

V2 explicitly counts, allocates and fills the selected indices in increasing
order. Empty selections remain valid. No numerical formula, index order, MPI
route, ownership policy or compiler/runtime check is changed. Only this solver
source file differs from v1. The runner now includes the last60 compiler-log
lines in Slurm output while retaining the full log and original failure status.

The exact packaged candidate passes clean J5 GCC13.2/MVAPICH4 checked and
optimized builds in an isolated bbserv directory. No solver run or Slurm
submission was performed during that compiler check. Three fresh local4-rank
J4/J6 optimized/checked/oracle full replays pass science and transform ownership
contracts;34 helper tests pass, including failure-diagnostic propagation.
All original v1 archives and failed-job evidence are preserved. Main is untouched.

These tests establish portability and preserve scientific qualification. They do
not establish speedup, legacy parity, memory parity or improved256 scaling.
The original49-run campaign and all performance/ownership gates are unchanged.
The new complete kit has not been uploaded or submitted. Only the single-source
compiler check was uploaded to an isolated temporary build directory.
