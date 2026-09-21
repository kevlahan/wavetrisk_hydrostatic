# Stage226: transfer, build, run and retrieve

This campaign tests the production rank-limit integration and higher task/Domain
ratios. It uses independent Stage226 paths and preserves earlier evidence.
Use the usual bbserv compiler/MPI/Python
with NumPy environment. The J5 fixture remains
`$HOME/stage187-fixture-20260911-145307`; the kit includes its own J4 fixture.

## 1. On murier: transfer

```bash
scp "$HOME/Desktop/wavetrisk-stage226-cluster-kit.tar.gz" kevlahan@bbserv.math.mcmaster.ca:~/
```

## 2. On bbserv: extract and submit fresh builds

```bash
tar -xzf "$HOME/wavetrisk-stage226-cluster-kit.tar.gz" -C "$HOME" &&
bash "$HOME/wavetrisk-stage226-cluster-kit/bootstrap_226.sh"
```

Watch the build:

```bash
source "$HOME/stage226-state.sh"
tail -F "$HOME/wavetrisk-226-build-${build_job}.log"
```

Wait for `PASS: six fresh J4/J5 builds...`, then press Ctrl-C to stop following
the log. This does not cancel the job. Do not launch the campaign if the build
fails. Bootstrap refuses to overwrite an existing source tree.

## 3. On bbserv: submit and follow the campaign

```bash
bash "$HOME/wavetrisk-stage226-cluster-kit/launch_226.sh"
```

Then:

```bash
source "$HOME/stage226-state.sh"
tail -F "$HOME/wavetrisk-226-${job}.log"
```

The launcher verifies sources and all six binaries, selects their build node,
and records the actual job ID. A briefly missing log is normal while the job
starts. It refuses duplicate submission. The first solver failure stops later
gates; collect the archive for either success or failure.

The 21 sequential gates start with the known J4 40/41 boundary, then test
80, 120 and 160 tasks over 40 Domains (2:1, 3:1 and 4:1). J5 repeats 160/167
before using 256 tasks over 160 Domains. Checked observers verify useful work
on all ranks. Optimized production and oracle runs contain no campaign bypass
or lifecycle file output. Six fresh builds cover J4/J5 baseline, production
and checked observer. The production candidate has no Stage223 selector.

Reserve one exclusive nk node with two 128-core sockets. All ranks are pinned
to distinct physical cores; the 256-task case uses both complete sockets.
Success ends with `PASS: production integration and higher-ratio gates...`.
Checkpoint saving is intentionally enabled to validate both reloads/remaps;
this is a correctness campaign, not a speed comparison.

## 4. On bbserv: create, verify and transfer the review archive to murier

After the job ends, press Ctrl-C to stop tail and run:

```bash
bash "$HOME/wavetrisk-stage226-cluster-kit/collect_226.sh"
```

This command checks the active user queue, creates the archive only if absent,
verifies its complete file inventory, prints SHA256 and executes scp FROM
bbserv TO `kevlahan@murier.math.mcmaster.ca:~/Desktop/`. It does not depend on
Slurm accounting or query a potentially purged job ID. Interrupted transfers
can be retried with the same command; the verified archive is reused.

Attach `stage226-review-JOBID.tar.gz` from murier's Desktop for review.
For an explicit archive-and-transfer block after confirming the job ended:

```bash
(
set -e
source "$HOME/stage226-state.sh"
archive="$HOME/stage226-review-${job}.tar.gz"
if [ ! -f "$archive" ]; then
    python3 "$kit/collect_cluster_logs_226.py" --stage 226 --job "$job" --results "$out" --tree "$tree"
fi
python3 "$kit/collect_cluster_logs_226.py" --verify "$archive"
sha256sum "$archive"
scp "$archive" kevlahan@murier.math.mcmaster.ca:~/Desktop/
)
```

Prefer `collect_226.sh`, which additionally checks that the job has finished.
