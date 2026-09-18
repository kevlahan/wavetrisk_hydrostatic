#!/bin/bash
# Run with bash; do not source this script into the login shell.
set -euo pipefail
state="$HOME/stage215-state.sh"
source "$state"
if [[ -n "${job:-}" ]]; then
    echo "A campaign job is already recorded ($job); inspect it before resubmitting." >&2
    exit 1
fi
python3 - "$tree/restriction-215-builds.json" <<'PY'
import json,sys
r=json.load(open(sys.argv[1]));assert r['status']=='passed','Builds have not passed'
assert len(r['builds'])==3 and r['fp_contract']=='off'
PY
out="$HOME/stage215-nk-results-$(date +%Y%m%d-%H%M%S)"
cd "$HOME"
job=$(sbatch --parsable "$kit/cluster_215.sbatch" "$kit" "$tree" "$fixture" "$out")
printf 'out=%q\njob=%q\n' "$out" "$job" >> "$state"
echo "Test job: $job"
echo "Results: $out"
echo "Log: $HOME/wavetrisk-215-${job}.log"
