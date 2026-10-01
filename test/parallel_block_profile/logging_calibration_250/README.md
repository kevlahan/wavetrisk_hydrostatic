# Stage250 J5 logging calibration on bb13

This kit measures how much the mirror-free candidate's per-rank storage
diagnostic contributes to the Stage249 main-relative regression. It extracts
the same sealed `bbserv-candidate5` and read-only `main` source archives used
in Stage249. Before building, it applies one sealed source patch that adds a
runtime switch to the diagnostic write in `src/domain.f90`. The Domain storage
assertion and numerical kernels remain in both modes. Loud and quiet candidate
runs use **the same binary**; only `WAVETRISK_QUIET_STORAGE_ABSENCE` changes.
The patched binary's loud mode must reproduce the prior qualified J5 state
and checkpoint hashes.

On `nk`/bb13, the kit validates full J5 83-rank runs for main, loud candidate
and quiet candidate. It checks twelve steps, two remaps, both VTK events,
physical-core placement, exact same-version output, and the established
main-to-candidate numerical checkpoint gate. It then runs two independent
A-B-B-A logging batches (loud-quiet-quiet-loud) and one A-B-B-A residual-gap
batch (main-quiet-quiet-main). Nine ordinary steps contribute one median per
full run; steps 1, 8 and 9 are excluded. The user-facing result is diagnostic,
with no automatic promotion or Domain removal.

The `bbserv-candidate5.tar.gz` and `rk249-main-legacy.tar.gz` archives and the
Stage187 J5 fixture must still be in your bbserv home directory. On the local
machine, transfer only the new small kit:

```bash
scp /Users/kevlahan/.codex/worktrees/7173/wavetrisk_hydrostatic/.stage250-validation/rk250-logging-performance-kit.tar.gz kevlahan@bbserv.math.mcmaster.ca:~/
```

On bbserv, extract and submit without editing any files:

```bash
tar -xzf "$HOME/rk250-logging-performance-kit.tar.gz" -C "$HOME"
bash "$HOME/rk250-logging-performance-kit/launch.sh"
```

The launcher checks the source archives, patch and J5 fixture before Slurm.
It prints the job ID, result path and log path. To inspect the job:

```bash
source "$HOME/rk250-latest.state"
sacct -j "$job" -X -o JobID,JobName%25,State,ExitCode,Elapsed,NodeList
tail -n 80 "$HOME/rk250-${job}.log"
```

After the job ends, create and transfer its verified review archive:

```bash
bash "$HOME/rk250-logging-performance-kit/collect.sh"
```
