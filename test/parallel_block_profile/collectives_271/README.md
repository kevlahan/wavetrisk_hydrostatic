Stage271 tests one isolated optimization and establishes production scaling
across the geometric-Domain boundary. Three consecutive payload-proof global
sums become one three-element integer reduction. Both byte-checksum scans,
the global checksum equality test, payload MPI, byte totals, error checks and
all numerical work remain. No audit timers are compiled into any solver.

The exact controls are accepted Stage266 ff740744 and legacy main 6a245db7
with its established checkpoint I/O hardening. Only
src/parallel_block_mpi.f90 differs in the candidate. No further Domain removal
or performance gain is claimed before the cluster comparison.

Stage270 measured 13.137 ms of mean-rank proof reduction time, with large
rank-to-rank variation. Fusion removes two collective calls per migration but
does not remove waiting for late ranks or the checksum scans. The full measured
cost is not an expected saving. The largest architectural priority remains the
Domain-bound producer/capture/restriction chain; this bounded change is tested
independently while obtaining the scaling baseline for that larger work.

One exclusive nk/bb13 allocation makes four serial builds using the same 23
physics objects: Stage266, main, optimized candidate, checked candidate.
It needs only the established fixture at
$HOME/stage187-fixture-20260911-145307 and the usual compiler/MPI modules.
No old campaign kit or result directory is required on bbserv.

There are twenty solver runs:

- Six validations: main, Stage266, candidate and checked candidate at 16 ranks;
  candidate oracles at 16 and 256 ranks. The 256-rank oracle covers 96 ranks
  without geometric Domains. All science, exact block checkpoints, 62-member
  VTK, restart/remap, storage, placement, 48 inverse proofs and 84 oracle
  installations remain required.
- At 83 ranks: main, Stage266, candidate, candidate, Stage266.
- At 160 ranks: the same order.
- At 256 ranks: Stage266, candidate, candidate, Stage266. Legacy main is not
  run beyond its geometric-Domain limit.

Each count has two independent runs per block version in A-B-B-A order;
legacy main has one anchor at each supported count. Small or mixed changes
remain inconclusive rather than triggering automatic repetition or promotion.
Exclude steps 1, 8 and 9. Report all paired gaps, run spreads, main gaps,
speedups and relative efficiencies from both 83- and 160-rank anchors.
Comparisons across counts use the same fixed-work fixture and source/binary;
physical-core placement spreads ranks across bb13's two sockets.

This is an intermediate single-node screen. Efficient execution at four or
more ranks per Domain means at least 640 ranks for this 160-Domain J5 fixture.
It still needs multi-node qualification and adequate block subdivision.
A favourable result at 256 does not establish the final objective.

Stage270 took 15 minutes 38 seconds for five builds and fifteen runs. Budget
roughly 20-25 minutes after allocation for this four-build/twenty-run campaign;
queue time is additional. The three-hour Slurm limit is a ceiling. Collection
now captures raw Slurm Start/End/ElapsedRaw/Elapsed automatically, alongside
runner UTC timestamps and monotonic elapsed time. If Slurm accounting is
unavailable, the archive records that explicitly and retains runner elapsed.

On murier:

```bash
scp /Users/kevlahan/.codex/worktrees/26d8/wavetrisk_hydrostatic/.stage271-validation/rk271-proof-scaling-kit-v1.tar.gz kevlahan@bbserv.math.mcmaster.ca:~/
```

On bbserv, extract and submit once:

```bash
tar -xzf "$HOME/rk271-proof-scaling-kit-v1.tar.gz" -C "$HOME"
bash "$HOME/rk271-proof-scaling-kit-v1/launch.sh"
```

Check status:

```bash
bash "$HOME/rk271-proof-scaling-kit-v1/status.sh"
```

After completion or failure, collect the review to murier Desktop:

```bash
bash "$HOME/rk271-proof-scaling-kit-v1/collect.sh"
```

The launcher records the unique output path and job ID in
$HOME/rk271-latest.json and rejects duplicates until the previous review is
collected. Status handles purged job IDs; collection preserves early failures.
Review archives omit source/binary/checkpoint/VTK payloads but preserve their
hashes, cluster comparison reports and raw solver timing/placement records.
The independent reviewer checks those reports; it cannot rerun omitted data.
