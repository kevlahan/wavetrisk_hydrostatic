# Stage256: residual velocity-owner gap with depth indexing

Stage255 improved c5 by about 70 ms per ordinary J5 step, but only about
28 ms per all-step came from native velocity source. The depth index also
helps programs active in c4. This kit builds both c4 and c5 from the same
sealed Stage251 source with the same one-file depth index, then compares
them directly. Their only feature-ladder difference is
`BLOCK_VELOCITY_OWNER`; all other selectors and the 83-core bb13 placement
match. This measures the remaining cost of enabling velocity ownership
after the common optimization.

Both builds use `make -j1`. The job runs unprofiled and profiled A-B-B-A
batches, each with two c4 and two c5 full J5 trajectories. Every run
must pass exact states, checkpoint, VTK and placement gates. The kit does
not promote any source. Keep Stage251 job 55256's source, results, kit and
state file on bbserv.

From the local machine:

```bash
scp /Users/kevlahan/.codex/worktrees/7173/wavetrisk_hydrostatic/.stage256-validation/rk256-depth-gap-kit-v1.tar.gz kevlahan@bbserv.math.mcmaster.ca:~/
```

On bbserv:

```bash
tar -xzf "$HOME/rk256-depth-gap-kit-v1.tar.gz" -C "$HOME"
bash "$HOME/rk256-depth-gap-kit-v1/launch.sh"
```

Check the job:

```bash
source "$HOME/rk256-latest.state"
sacct -j "$job" -X -o JobID,JobName%25,State,ExitCode,Elapsed,NodeList
tail -n 80 "$HOME/rk256-${job}.log"
```

After it ends, verify and transfer the archive to the local Desktop:

```bash
bash "$HOME/rk256-depth-gap-kit-v1/collect.sh"
```
