#!/bin/bash
set -euo pipefail
state="$HOME/stage230-state.sh"
source "$state"
[[ -z "${job:-}" ]] || { echo "Campaign job $job already recorded; inspect before resubmitting." >&2; exit 1; }
build_host=$(python3 - "$kit" "$tree" <<'PY'
import json,re,sys
from pathlib import Path
kit,tree=map(Path,sys.argv[1:]);sys.path.insert(0,str(kit))
from capability_223 import verify_kit,digest,verify_sources
from build_capability_230 import KEYS,source_key,source_dir
verify_kit(kit);m=json.loads((kit/'stage230-sources.json').read_text())
r=json.loads((tree/'capability-230-builds.json').read_text())
assert r['status']=='passed' and set(r['builds'])==set(KEYS) and r['fp_contract']=='off'
assert r['manifest_sha256']==digest(kit/'stage230-sources.json')
for kind in ('baseline','candidate','checked'):
    verify_sources(tree/source_dir(kind),m[source_key(kind)])
for b in r['builds'].values():assert digest(tree/b['binary'])==b['sha256']
host=r['build_host']['node'].split('.')[0]
assert re.fullmatch(r'[A-Za-z0-9_-]+',host) and host != 'bb04'
print(host)
PY
)
out="$HOME/stage230-nk-results-$(date +%Y%m%d-%H%M%S)"
cd "$HOME"
job=$(sbatch --parsable --nodelist="$build_host" "$kit/cluster_230.sbatch" "$kit" "$tree" "$fixture" "$out")
job=${job%%;*}
[[ "$job" =~ ^[0-9]+$ ]] || { echo 'Unexpected sbatch job ID' >&2; exit 1; }
printf 'out=%q\njob=%q\n' "$out" "$job" >> "$state"
echo "Campaign job: $job on $build_host"
echo "Results: $out"
echo "Log: $HOME/wavetrisk-230-${job}.log"
