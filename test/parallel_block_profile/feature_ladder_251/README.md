# Stage251 J5 early-feature ladder on bb13

Stage250 established that suppressing the candidate's storage diagnostic
write does not recover its roughly 1.12 s/step gap to `main`. This kit tests
the six earlier owner features that were enabled in both Stage248 variants.
It uses the same sealed `bbserv-candidate5` source for all seven candidate
builds, with `MIRROR_FREE=off`, `BLOCK_RK=off` and later native features off.
Stage `c0` has all six early features off; `c1` through `c6` enable owner
columns, primitives, diffusion, mass owner, velocity owner and initial
tendency in dependency order. It also rebuilds the established read-only
`main` export as a numerical and timing anchor.

On `nk`/bb13, each binary gets a full 83-rank J5 validation run before any
timing comparison. Each must complete twelve steps, two remaps, both VTK
events, the same adaptive trajectory and the established main-to-candidate
numerical checkpoint gate. The kit screens each adjacent median once, then
confirms the largest positive jump with two independent A-B-B-A batches.
Each full run contributes one median of the nine ordinary steps; steps 1,
8 and 9 are excluded. Screens are for selecting the pair, not a final
feature-level estimate. No code is promoted automatically.

The `bbserv-candidate5.tar.gz` and `rk249-main-legacy.tar.gz` archives and
the Stage187 J5 fixture must still be in your bbserv home directory. On the
local machine, transfer only this small kit:

```bash
scp /Users/kevlahan/.codex/worktrees/7173/wavetrisk_hydrostatic/.stage251-validation/rk251-feature-ladder-kit.tar.gz kevlahan@bbserv.math.mcmaster.ca:~/
```

On bbserv, extract and submit without editing files:

```bash
tar -xzf "$HOME/rk251-feature-ladder-kit.tar.gz" -C "$HOME"
bash "$HOME/rk251-feature-ladder-kit/launch.sh"
```

To inspect the job:

```bash
source "$HOME/rk251-latest.state"
sacct -j "$job" -X -o JobID,JobName%25,State,ExitCode,Elapsed,NodeList
tail -n 80 "$HOME/rk251-${job}.log"
```

After the job ends, create and transfer its verified review archive:

```bash
bash "$HOME/rk251-feature-ladder-kit/collect.sh"
```
