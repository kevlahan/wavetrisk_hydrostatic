# Stage225: transfer, build, run and retrieve

This fresh campaign adds the locally validated empty-Domain oracle extent fix
to the Stage224 checkpoint correction and repaired Slurm placement launcher.
It uses independent Stage225 paths; keep Stage223/224 evidence and binaries. Use the usual bbserv compiler/MPI/Python
with NumPy environment. The J5 fixture remains
`$HOME/stage187-fixture-20260911-145307`; the kit includes its own J4 fixture.

## 1. On murier: transfer

```bash
scp "$HOME/Desktop/wavetrisk-stage225-cluster-kit.tar.gz" kevlahan@bbserv.math.mcmaster.ca:~/
```

## 2. On bbserv: extract and submit fresh builds

```bash
tar -xzf "$HOME/wavetrisk-stage225-cluster-kit.tar.gz" -C "$HOME" &&
bash "$HOME/wavetrisk-stage225-cluster-kit/bootstrap_225.sh"
```

Watch the build:

```bash
source "$HOME/stage225-state.sh"
tail -F "$HOME/wavetrisk-225-build-${build_job}.log"
```

Wait for `PASS: six fresh J4/J5 builds...`, then press Ctrl-C to stop following
the log. This does not cancel the job. Do not launch the campaign if the build
fails. Bootstrap refuses to overwrite an existing source tree.

## 3. On bbserv: submit and follow the campaign

```bash
bash "$HOME/wavetrisk-stage225-cluster-kit/launch_225.sh"
```

Then:

```bash
source "$HOME/stage225-state.sh"
tail -F "$HOME/wavetrisk-225-${job}.log"
```

The launcher verifies sources and all six binaries, selects their build node,
and records the actual job ID. A briefly missing log is normal while the job
starts. It refuses duplicate submission. The first solver failure stops later
gates; collect the archive for either success or failure.

The eleven gates start at J4 40/41 ranks, including the empty-root case, then
advance to J5 160/167. Reserve one exclusive nk node with two 128-core sockets.
All ranks are pinned to distinct physical cores across both sockets. Success
ends with `PASS: 41/40 and 167/160 capability gates...`.

## 4. On bbserv: create, verify and transfer the review archive to murier

After the job ends, press Ctrl-C to stop tail and run:

```bash
bash "$HOME/wavetrisk-stage225-cluster-kit/collect_225.sh"
```

This command checks the active user queue, creates the archive only if absent,
verifies its complete file inventory, prints SHA256 and executes scp FROM
bbserv TO `kevlahan@murier.math.mcmaster.ca:~/Desktop/`. It does not depend on
Slurm accounting or query a potentially purged job ID. Interrupted transfers
can be retried with the same command; the verified archive is reused.

Attach `stage225-review-JOBID.tar.gz` from murier's Desktop for review.
For an explicit archive-and-transfer block after confirming the job ended:

```bash
(
set -e
source "$HOME/stage225-state.sh"
archive="$HOME/stage225-review-${job}.tar.gz"
if [ ! -f "$archive" ]; then
    python3 "$kit/collect_cluster_logs_225.py" --stage 225 --job "$job" --results "$out" --tree "$tree"
fi
python3 "$kit/collect_cluster_logs_225.py" --verify "$archive"
sha256sum "$archive"
scp "$archive" kevlahan@murier.math.mcmaster.ca:~/Desktop/
)
```

Prefer `collect_225.sh`, which additionally checks that the job has finished.
