# Stage254: measure the remaining velocity-owner cost

The Stage253 empty-collective change recovered only about 15 ms of the
roughly 406 ms c4-to-c5 ordinary-step regression. This kit builds a timing
probe from the verified Stage251 c5 source and compares that build with
the exact Stage251 c5 binary on bb13. The probe records time in resident
fetch, dependency reads, MPI inside those reads, producer loops, and
final delivery. It reports by dependency depth at each existing profile
window; it makes no numerical algorithm change.

One unprofiled and one profiled A-B-B-A batch measure the probe's timing
cost. Every run must pass the established J5 state, checkpoint, VTK and
83-physical-core placement gates. The build uses `make -j1`. Nothing is
promoted automatically.

Keep Stage251 job 55256's result directory, source tree, kit and state
file on bbserv. On the local machine:

```bash
scp /Users/kevlahan/.codex/worktrees/7173/wavetrisk_hydrostatic/.stage254-validation/rk254-owner-breakdown-kit-v1.tar.gz kevlahan@bbserv.math.mcmaster.ca:~/
```

On bbserv, submit without editing files:

```bash
tar -xzf "$HOME/rk254-owner-breakdown-kit-v1.tar.gz" -C "$HOME"
bash "$HOME/rk254-owner-breakdown-kit-v1/launch.sh"
```

To inspect the submitted job:

```bash
source "$HOME/rk254-latest.state"
sacct -j "$job" -X -o JobID,JobName%25,State,ExitCode,Elapsed,NodeList
tail -n 80 "$HOME/rk254-${job}.log"
```

After the job ends, transfer its review archive to the local Desktop:

```bash
bash "$HOME/rk254-owner-breakdown-kit-v1/collect.sh"
```
