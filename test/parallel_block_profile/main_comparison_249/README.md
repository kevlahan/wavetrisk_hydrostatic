# Direct J5 comparison: current migration versus read-only main

This kit builds the qualified mirror-free migration and a read-only export of
`main` on the same exclusive `nk` node, bb13. The main export has only the
established Stage190 checkpoint I/O hardening needed on bbserv. Both builds use
the same compiler, FP contraction setting, J5 restart fixture and external
Simple-physics objects. Main is never checked out or modified.

The last direct comparison was Stage229, job 54737: at 83 ranks the then-current
block candidate averaged 0.758268 s per ordinary step versus 0.542883 s for
legacy, a 215.386 ms gap. Those runs predate the RK migration. Legacy observer
overhead failed calibration, so its phase decomposition is not an accepted
breakdown of that gap. This kit makes a new uninstrumented comparison.

It verifies the J5 trajectory at 16 and 83 ranks, including 12 steps, both
remaps, two VTK events, same-version exact output and the unchanged
cross-version checkpoint numerical gate. The candidate must also reproduce
the qualified v7 J5 state and checkpoint hashes. At 83 ranks it measures
separate full-run medians in A-B-B-A order, with at most one extra batch if the
3% guard is ambiguous. Ordinary steps are 2-7 and 10-12; steps within a run
are not independent samples. The result does not automatically promote code.

On the local machine, transfer the two new archives (the existing
`bbserv-candidate5.tar.gz` must remain in your bbserv home directory):

```bash
scp /Users/kevlahan/.codex/worktrees/7173/wavetrisk_hydrostatic/.stage249-validation/rk249-main-legacy.tar.gz /Users/kevlahan/.codex/worktrees/7173/wavetrisk_hydrostatic/.stage249-validation/rk249-main-performance-kit.tar.gz kevlahan@bbserv.math.mcmaster.ca:~/
```

On bbserv, extract and submit; no file editing or manual `sbatch` arguments:

```bash
tar -xzf "$HOME/rk249-main-performance-kit.tar.gz" -C "$HOME"
bash "$HOME/rk249-main-performance-kit/launch.sh"
```

The launcher checks both source archive hashes and the J5 input and seed
before submission. It prints the job ID, results path and log path. To inspect
the latest job:

```bash
source "$HOME/rk249-latest.state"
sacct -j "$job" -X -o JobID,JobName%25,State,ExitCode,Elapsed,NodeList
tail -n 80 "$HOME/rk249-${job}.log"
```

After the job ends, collect and transfer its verified review archive:

```bash
bash "$HOME/rk249-main-performance-kit/collect.sh"
```
