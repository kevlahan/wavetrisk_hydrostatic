#!/bin/bash
set -euo pipefail
state="$HOME/stage224-state.sh"
source "$state"
[[ -z "${job:-}" ]] || { echo "Campaign job $job already recorded; inspect before resubmitting." >&2; exit 1; }
build_host=$(python3 - "$kit" "$tree" <<'PY'
import json,re,sys
from pathlib import Path
kit,tree=map(Path,sys.argv[1:]);sys.path.insert(0,str(kit))
from capability_223 import verify_kit,digest,verify_sources
from build_capability_224 import KEYS
verify_kit(kit);m=json.loads((kit/'stage224-sources.json').read_text())
r=json.loads((tree/'capability-224-builds.json').read_text())
assert r['status']=='passed' and set(r['builds'])==set(KEYS) and r['fp_contract']=='off'
assert r['manifest_sha256']==digest(kit/'stage224-sources.json')
for kind in ('reference','candidate'):verify_sources(tree/kind,m[kind+'_sources'])
for b in r['builds'].values():assert digest(tree/b['binary'])==b['sha256']
host=r['build_host']['node'].split('.')[0]
assert re.fullmatch(r'[A-Za-z0-9_-]+',host)
print(host)
PY
)
out="$HOME/stage224-nk-results-$(date +%Y%m%d-%H%M%S)"
cd "$HOME"
job=$(sbatch --parsable --nodelist="$build_host" "$kit/cluster_224.sbatch" "$kit" "$tree" "$fixture" "$out")
job=${job%%;*}
[[ "$job" =~ ^[0-9]+$ ]] || { echo 'Unexpected sbatch job ID' >&2; exit 1; }
printf 'out=%q\njob=%q\n' "$out" "$job" >> "$state"
echo "Campaign job: $job on $build_host"
echo "Results: $out"
echo "Log: $HOME/wavetrisk-224-${job}.log"
