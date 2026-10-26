# Stage325 native verification v5

This package verifies saved job55569 captures. It launches zero solver, compiler
or performance runs. Job55574 and every predecessor input/package remain intact.
This is one separately authorized 1,800-second allocation, not an extension or
automatic replay of the exhausted earlier budget.

`launch.sh` first verifies all sealed files, runs the packaged offline test suite
with a bounded log, verifies the sealed files again, then submits one exclusive
bb13/nk allocation: one task,16CPUs,30minutes,16GiB process virtual-memory limits.
Only a CPU-bound Python verification step is launched through srun. The Python
supervisor enforces its own earlier deadline and aggregate memory/data limits.
There is no automatic submission retry, resumed run, budget extension or solver.

The new state is `$HOME/rk325-native-v5-latest.json`. Each launch uses a fresh
private work directory. `status.sh` requires confirmed scheduler/accounting
evidence; missing old squeue IDs fall back to sacct without hiding service errors.

After terminal state, `collect.sh` creates
`$HOME/rk325-native-v5-review-JOB.tar.gz`. The logical archive cap is32MiB.
It returns sealed code, complete manifests and bounded verification JSON/logs.
SQLite caches, captures, checkpoints, TSVs and full vectors remain on bbserv.
All returned files have an exact SHA256 inventory. Partial phase certificates
remain labeled by their staging path and never establish reusable completion.

A passing replay still requires independent evidence review. This campaign
establishes no elapsed solver speedup, performance qualification or promotion.
