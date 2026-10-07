Stage264 establishes the new starting point: exact Stage230 versus legacy main
on nk/bb13. It makes no solver changes. Sources and scripts are self-contained;
the only existing bbserv asset required is the established J5 fixture at
`$HOME/stage187-fixture-20260911-145307`. It requires no prior campaign kits,
results directories, source trees, or state files.

The job builds Stage230 optimized, legacy main optimized, and Stage230 checked
with serial make. One external Simple-physics build supplies the same 23 objects
to all solvers. Main uses the established checkpoint I/O hardening only. Each
solver retains its original build flags, with floating-point contraction off.

Six validation runs check J5 correctness at 16 ranks and Stage230 production and
oracle operation at 167 ranks, including ranks without geometric Domains. Eight
unprofiled 83-rank runs form two independent main/Stage230/Stage230/main batches.
Rank placement matches Stage230: 42/41 physical cores across bb13's sockets.
All runs check twelve states, both reloads, two remaps, checkpoint numerical
limits, 62 VTK members, and the appropriate Domain-tendency storage policy.
Main-relative checkpoint limits are unchanged; block comparisons remain exact.
Ordinary medians exclude steps 1, 8, and 9; validation timings are excluded.

On the local machine:

```bash
scp /Users/kevlahan/.codex/worktrees/26d8/wavetrisk_hydrostatic/.stage264-validation/rk264-stage230-reset-kit-v1.tar.gz kevlahan@bbserv.math.mcmaster.ca:~/
```

On bbserv, with the usual compiler/MPI modules available:

```bash
tar -xzf "$HOME/rk264-stage230-reset-kit-v1.tar.gz" -C "$HOME"
bash "$HOME/rk264-stage230-reset-kit-v1/launch.sh"
```

Check progress:

```bash
bash "$HOME/rk264-stage230-reset-kit-v1/status.sh"
```

After it ends, collect and transfer the review to murier Desktop:

```bash
bash "$HOME/rk264-stage230-reset-kit-v1/collect.sh"
```

The launcher refuses another submission until the previous review has been
collected. A fresh directory is created for every accepted submission. Collection
works after a failed build/run as well as after success, and queries active jobs
without requesting a possibly purged job ID from squeue. Repeating collection
verifies and retransfers the existing archive. Full solver outputs remain in the
results directory printed in the saved `$HOME/rk264-latest.json` state.

The review archive omits binary/checkpoint/VTK payloads. The cluster checks their
bytes; the local reviewer checks recorded hashes, numerical reports, physical
placement, storage policy, source/build identities, raw states, and recomputes
the timing result. Nothing is promoted automatically.
