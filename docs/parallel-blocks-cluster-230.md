# Stage230: qualify removal of Domain tendency storage

The kit contains the committed Stage229 baseline (5c457a5d), the seven locally
validated Stage230 source edits, and frozen build/run/retrieval helpers. It uses
no Git bundle prerequisites and does not modify any Git branch, including main.

Three fresh J5 builds: baseline, optimized candidate and checked candidate.
Both versions retain STAGE202=on, STAGE203=on and FP contraction off, and share
identical external physics objects. The checked build uses DEBUG=check.
Build and campaign use the same exclusive nk host: two 128-core sockets, 256
physical cores. bb04 is excluded. Placement probes verify actual rank affinity;
16 ranks use one socket, 83 and 167 alternate sockets (42/41 and 84/83).

Correctness runs first: baseline, checked candidate and full candidate oracles
at 16 ranks; those three plus production candidate at 167 ranks / 160 Domains.
Each run must complete the 12-step restart/adaptation/remap trajectory. Exact
checkpoint and all 62 VTK payload comparisons use a baseline at the SAME rank
count. Cross-count checkpoint comparison also uses the established tolerance
policy. All-rank storage census must report 0/N twice for production/checked,
and N/N twice for oracle runs. Existing intermediate oracles remain enabled in
the oracle cases. The deliberately rejected legacy fallback was tested locally.

The established J5 fixture is $HOME/stage187-fixture-20260911-145307. It is copied
and hashed, with max_level=7, resume=3, time_end=.315, dt_write=.01 and CP_EVERY=1.
Output/reload events are retained to test the changed storage lifecycle. Timing
uses high-precision ordinary-step medians, excluding startup, saving/checkpoint
step and first post-reload step (steps 1, 8, 9); ordinary remaps remain included.
Timesteps within a run are not treated as independent repetitions.

Then one uninstrumented baseline/candidate pair at 16 ranks and one at 83 ranks.
A >=3% slowdown stops later performance runs for review. Otherwise confirmation
at 83 ranks is A-B-B-A (two independent runs/version). Repeat that batch only
if adjacent-pair decisions straddle the 3% regression threshold or the aggregate
is within 0.5 percentage points of it; at most two batches. The ordinary path
is 15 solver runs, maximum 19. No minimum speedup is required for this cleanup.
All numerical/storage gates apply to timing runs too. No broad instrumentation,
automatic promotion, further cache variant, or complete Domain-free claim.

## 1. On murier: transfer the kit

```bash
scp "$HOME/Desktop/wavetrisk-stage230-cluster-kit.tar.gz" \
    kevlahan@bbserv.math.mcmaster.ca:~/
```

## 2. On bbserv: extract and submit the build

Use the same compiler/MPI/Python environment as the successful recent campaigns
(Python with numpy, mpif90, git, make, numactl, lscpu and zstd available).

```bash
(
set -e
cd "$HOME"
test ! -e wavetrisk-stage230-cluster-kit || { echo 'Stage230 kit already exists; stop for inspection.'; exit 1; }
tar -xzf wavetrisk-stage230-cluster-kit.tar.gz
bash "$HOME/wavetrisk-stage230-cluster-kit/bootstrap_230.sh"
)
```

```bash
source "$HOME/stage230-state.sh"
tail -F "$HOME/wavetrisk-230-build-${build_job}.log"
```

Wait for `PASS: three fresh builds`. Press Ctrl-C to stop tail (not the job).
A log may take a moment to appear; tail -F keeps waiting.

## 3. On bbserv: submit and follow the campaign

```bash
bash "$HOME/wavetrisk-stage230-cluster-kit/launch_230.sh"
source "$HOME/stage230-state.sh"
tail -F "$HOME/wavetrisk-230-${job}.log"
```

The launcher checks all build/source identities and selects the build host for
native compilation compatibility. It refuses duplicate submissions. Wait for
`COMPLETE: ...`, or collect failure evidence if the job exits with an error.
Press Ctrl-C to stop tail before running the next command.

## 4. On bbserv: create, verify and transfer the review archive to murier

```bash
bash "$HOME/wavetrisk-stage230-cluster-kit/collect_230.sh"
```

That wrapper executes the following complete creation and scp sequence:

```bash
(
set -e
source "$HOME/stage230-state.sh"
active=$(squeue -h -u "$USER" -o '%A')
if printf '%s\n' "$active" | awk -v wanted="$job" '$1==wanted {found=1} END {exit !found}'; then
    echo "Job $job is still queued/running; collect after it ends."
    exit 1
fi
archive="$HOME/stage230-review-${job}.tar.gz"
if [ ! -f "$archive" ]; then
    python3 "$kit/collect_cluster_logs_230.py" \
        --stage 230 --job "$job" --results "$out" --tree "$tree"
fi
python3 "$kit/collect_cluster_logs_230.py" --verify "$archive"
sha256sum "$archive"
scp "$archive" kevlahan@murier.math.mcmaster.ca:~/Desktop/
)
```

Attach stage230-review-JOBID.tar.gz from murier's Desktop. Repeating collection
reuses and verifies the existing archive before retrying scp, so an interrupted
transfer does not require recreating it. Queue checks do not depend on retention
of old Slurm job IDs or sacct. Keep the full cluster result directory until review;
the small archive includes logs, hashes and exact-comparison reports rather than
large raw checkpoints/binaries/VTK files.
