# Stage252: profile the confirmed J5 velocity-owner jump

Stage251 job 55256 confirmed that enabling `BLOCK_VELOCITY_OWNER` adds about
0.386 s per ordinary step at 83 ranks. This kit reuses the exact verified
`c4` and `c5` binaries and input from that job, so it does not rebuild the
solver. It runs one unprofiled and one profiled A-B-B-A batch on bb13. The
existing parallel-block profiler records native velocity-source and gradient
phase times, while full J5 state, checkpoint, VTK and placement gates protect
the comparison. The unprofiled batch calibrates profiling overhead in the
same allocation. No change is promoted automatically.

Keep the Stage251 job's result directory, source tree, kit and
`$HOME/rk251-job-55256.state` on bbserv. On the local machine:

```bash
scp /Users/kevlahan/.codex/worktrees/7173/wavetrisk_hydrostatic/.stage252-validation/rk252-velocity-profile-kit-v4.tar.gz kevlahan@bbserv.math.mcmaster.ca:~/
```

On bbserv, submit without editing files:

```bash
tar -xzf "$HOME/rk252-velocity-profile-kit-v4.tar.gz" -C "$HOME"
bash "$HOME/rk252-velocity-profile-kit-v4/launch.sh"
```

To inspect the job:

```bash
source "$HOME/rk252-latest.state"
sacct -j "$job" -X -o JobID,JobName%25,State,ExitCode,Elapsed,NodeList
tail -n 80 "$HOME/rk252-${job}.log"
```

After the job ends, verify and transfer its review archive to the local
Desktop:

```bash
bash "$HOME/rk252-velocity-profile-kit-v4/collect.sh"
```
