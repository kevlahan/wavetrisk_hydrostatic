# Stage255: owner depth-index trial on bb13

The Stage254 probe found 0.111 s/step in mean-rank producer traversal and
0.152 s/step inside positive-depth MPI reads. Some of the MPI time may be
waiting for slower producers. This kit tests a stable depth index: compile
each owner's operation positions into depth buckets, then execute only the
bucket for the current dependency depth. The original order of operations
within each depth is preserved, as are all communication and numerical
calls. The kit transforms the exact sealed Stage251 c5 source; it does not
include the separate Stage253 empty-collective trial.

The cluster job builds the one-file candidate with `make -j1`, then compares
it with the verified Stage251 c5 binary in unprofiled and profiled A-B-B-A
batches at 83 physical cores on bb13. Every run must pass exact J5 state,
checkpoint, VTK and placement gates. No source is promoted automatically.
The local four-rank J4 replay passed 12 exact states, two remaps, the exact
checkpoint hash and all 62 byte-identical VTK members. Local timing is not
used to judge performance.

Keep Stage251 job 55256's results, source, kit and state file on bbserv.
From the local machine:

```bash
scp /Users/kevlahan/.codex/worktrees/7173/wavetrisk_hydrostatic/.stage255-validation/rk255-owner-depth-kit-v1.tar.gz kevlahan@bbserv.math.mcmaster.ca:~/
```

On bbserv, submit without editing files:

```bash
tar -xzf "$HOME/rk255-owner-depth-kit-v1.tar.gz" -C "$HOME"
bash "$HOME/rk255-owner-depth-kit-v1/launch.sh"
```

Check the recorded job:

```bash
source "$HOME/rk255-latest.state"
sacct -j "$job" -X -o JobID,JobName%25,State,ExitCode,Elapsed,NodeList
tail -n 80 "$HOME/rk255-${job}.log"
```

After it ends, verify and transfer the review archive to the local Desktop:

```bash
bash "$HOME/rk255-owner-depth-kit-v1/collect.sh"
```
