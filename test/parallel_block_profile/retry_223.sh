#!/bin/bash
# Run with bash, never source into the login shell.
set -euo pipefail
repair=$(cd -- "$(dirname -- "$0")" && pwd)
state="$HOME/stage223-state.sh"
source "$state"
# A retry must never overwrite a later campaign or hide evidence from solver runs.
if [[ -n "${job:-}" && "$job" != 54724 ]]; then
    echo "Campaign job $job is already recorded; no changes or duplicate submission."
    exit 0
fi
active=$(squeue -h -u "$USER" -o '%A')
if awk '$1 == 54724 {found=1} END {exit !found}' <<< "$active"; then
    echo 'Job 54724 is still active; wait until it finishes.' >&2
    exit 1
fi
if [[ "${job:-}" == 54724 ]]; then
    python3 - "$out/results.json" <<'PY'
import json,sys
r=json.load(open(sys.argv[1]))
assert r['stage']==223 and r['status']=='failed', 'Expected failed Stage223 result'
assert not r['runs'] and not r['schedule'], 'Solver runs exist; inspect before retrying'
assert 'probe' in r.get('error','').lower(), 'Expected placement-probe failure'
print('Confirmed failed placement probe; all old results will be preserved.')
PY
else
    [[ -f "$HOME/stage223-state.before-placement-repair.sh" ]] || { echo 'Missing repair state backup.' >&2; exit 1; }
fi
python3 "$repair/repair_223.py" "$kit"
if [[ "${job:-}" == 54724 ]]; then
    backup="$HOME/stage223-state.before-placement-repair.sh"
    if [[ -e "$backup" ]]; then
        cmp -s "$state" "$backup" || { echo 'Existing state backup differs; inspect before retrying.' >&2; exit 1; }
    else
        cp -p "$state" "$backup"
    fi
    temporary=$(mktemp "$HOME/stage223-state.repair.XXXXXX")
    printf 'kit=%q\ntree=%q\nfixture=%q\nbuild_job=%q\n' "$kit" "$tree" "$fixture" "$build_job" > "$temporary"
    mv "$temporary" "$state"
fi
bash "$kit/launch_223.sh"
