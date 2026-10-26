# Stage325 native verification v5p1 — candidate awaiting reviewed approval

This package verifies saved job55569 captures. It launches zero solver, compiler
or performance runs. Job55575 failed in0schedulerseconds before the workflow
because its srun process did not inherit the batch-shell address-space limit.
The corrected runner installs and strictly verifies16GiB inside the actual step.
Job55575, its v5 kit/state/log and every predecessor input/package remain intact.
The earlier budget remains2656seconds consumed; the new1800-second budget was
not consumed by55575. Launching this corrected candidate requires explicit user
approval after reviewing the repair and its adversarial tests. No retry is automatic.

`launch.sh` first verifies all sealed files, runs the packaged offline test suite
with a bounded log, verifies the sealed files again, then submits one exclusive
bb13/nk allocation: one task,16CPUs,30minutes,16GiB process virtual-memory limits.
Only a CPU-bound Python verification step is launched through srun. The Python
supervisor enforces its own earlier deadline and aggregate memory/data limits.
There is no automatic submission retry, resumed run, budget extension or solver.

The new state is `$HOME/rk325-native-v5p1-latest.json`. Each launch uses a fresh
private work directory. `status.sh` requires confirmed scheduler/accounting
evidence; missing old squeue IDs fall back to sacct without hiding service errors.

After terminal state, `collect.sh` creates
`$HOME/rk325-native-v5p1-review-JOB.tar.gz`. The logical archive cap is32MiB.
It returns sealed code, complete manifests and bounded verification JSON/logs.
SQLite caches, captures, checkpoints, TSVs and full vectors remain on bbserv.
All returned files have an exact SHA256 inventory. Partial phase certificates
remain labeled by their staging path and never establish reusable completion.

A passing replay still requires independent evidence review. This campaign
establishes no elapsed solver speedup, performance qualification or promotion.
