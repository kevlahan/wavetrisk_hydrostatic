#!/bin/bash
# Create/reuse and verify the archive, then transfer FROM bbserv to murier.
set -euo pipefail
source "$HOME/stage223-state.sh"
active=$(squeue -h -u "$USER" -o '%A')
if printf '%s\n' "$active" | awk -v wanted="$job" '$1==wanted {found=1} END {exit !found}'; then
    echo "Job $job is still queued/running; collect after it ends." >&2
    exit 1
fi
archive="$HOME/stage223-review-${job}.tar.gz"
if [[ ! -f "$archive" ]]; then
    python3 "$kit/collect_cluster_logs_223.py" --stage 223 --job "$job" --results "$out" --tree "$tree"
fi
python3 "$kit/collect_cluster_logs_223.py" --verify "$archive"
sha256sum "$archive"
scp "$archive" kevlahan@murier.math.mcmaster.ca:~/Desktop/
