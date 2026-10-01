#!/bin/bash
set -euo pipefail
state="$HOME/rk252-latest.state"
[[ -f "$state" ]] || { echo "No Stage252 submission recorded: $state" >&2; exit 1; }
source "$state"
active=$(squeue -h -u "$USER" -o '%A')
if printf '%s\n' "$active" | awk -v wanted="$job" '$1==wanted {found=1} END {exit !found}'; then
  echo "Job $job is still queued or running; collect after it ends." >&2
  exit 1
fi
archive="$HOME/rk252-review-${job}.tar.gz"
log="$HOME/rk252-${job}.log"
if [[ ! -f "$archive" ]]; then
  python3 "$kit/collect.py" "$out" "$archive" --slurm-log "$log" --state "$state"
fi
python3 "$kit/collect.py" --verify "$archive"
sha256sum "$archive"
scp "$archive" kevlahan@murier.math.mcmaster.ca:~/Desktop/
echo "Transferred review: $archive"
