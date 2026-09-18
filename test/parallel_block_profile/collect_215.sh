#!/bin/bash
# Atomic archive creation/reuse, verification, then transfer FROM bbserv to murier.
set -euo pipefail
source "$HOME/stage215-state.sh"
# Query by user instead of a purged job ID. Scheduler errors remain real errors.
active=$(squeue -h -u "$USER" -o '%A')
if printf '%s\n' "$active" | awk -v wanted="$job" '$1==wanted {found=1} END {exit !found}'; then
    echo "Job $job is still queued/running; collect after it ends." >&2
    exit 1
fi
archive="$HOME/stage215-review-${job}.tar.gz"
if [[ ! -f "$archive" ]]; then
    python3 "$kit/collect_cluster_logs_215.py" --stage 215 --job "$job" --results "$out" --tree "$tree"
fi
python3 "$kit/collect_cluster_logs_215.py" --verify "$archive"
sha256sum "$archive"
scp "$archive" kevlahan@murier.math.mcmaster.ca:~/Desktop/
