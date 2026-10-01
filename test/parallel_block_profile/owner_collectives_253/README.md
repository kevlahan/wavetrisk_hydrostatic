# Stage253: owner-program collective trial

Stage252 attributed 89.2% of the profiled c4-to-c5 timestep gap to native
velocity-source execution. This kit tests a narrow change: omit the resident
fetch when the compiled program has no resident producer, and omit the
zero-input depth-zero read. It copies the sealed Stage251 source, replaces
only `src/parallel_block_program.f90`, and builds c5 with `make -j1`.

The cluster job compares the verified Stage251 c5 binary with the optimized
binary in unprofiled and profiled A-B-B-A batches at 83 ranks on bb13. Every
run must satisfy the established J5 trajectory, checkpoint, VTK and core
placement gates. No source is promoted automatically. A local four-rank J4
replay already passed exact state, checkpoint and VTK comparisons; local
timing is not used to judge performance.

Keep the Stage251 job 55256 results and its state file on bbserv. On the
local machine, copy the sealed kit:

```bash
scp /Users/kevlahan/.codex/worktrees/7173/wavetrisk_hydrostatic/.stage253-validation/rk253-owner-collectives-kit-v2.tar.gz kevlahan@bbserv.math.mcmaster.ca:~/
```

On bbserv, submit without editing files:

```bash
tar -xzf "$HOME/rk253-owner-collectives-kit-v2.tar.gz" -C "$HOME"
bash "$HOME/rk253-owner-collectives-kit-v2/launch.sh"
```

Check the recorded job if needed:

```bash
source "$HOME/rk253-latest.state"
sacct -j "$job" -X -o JobID,JobName%25,State,ExitCode,Elapsed,NodeList
tail -n 80 "$HOME/rk253-${job}.log"
```

After the job ends, verify and transfer its review archive to the local
Desktop:

```bash
bash "$HOME/rk253-owner-collectives-kit-v2/collect.sh"
```
