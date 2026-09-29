#!/bin/bash
# Execute with bash; never source this script into the login shell.
set -euo pipefail
kit=$(cd "$(dirname "$0")" && pwd)
tree="$HOME/wt-stage230"
fixture="$HOME/stage187-fixture-20260911-145307"
state="$HOME/stage230-state.sh"
[[ ! -e "$tree" ]] || { echo "Existing source tree: $tree; nothing overwritten." >&2; exit 1; }
[[ -f "$fixture/simple.in" ]] || { echo "Missing J5 fixture: $fixture" >&2; exit 1; }
python3 - "$kit" "$tree" <<'PY'
import json,sys,tarfile,tempfile
from pathlib import Path
kit,tree=map(Path,sys.argv[1:]);sys.path.insert(0,str(kit))
from capability_223 import verify_kit,verify_sources
verify_kit(kit)
m=json.loads((kit/'stage230-sources.json').read_text())
with tempfile.TemporaryDirectory(dir=tree.parent,prefix='.stage230-bootstrap-') as tmp:
    root=Path(tmp)/'root';source=root/'reference';source.mkdir(parents=True)
    with tarfile.open(kit/'stage230-source.tar.gz') as t:t.extractall(source,filter='data')
    verify_sources(source,m['reference_sources'])
    root.rename(tree)
print('PASS: frozen Stage230 baseline source snapshot, no Git bundle prerequisites')
PY
if [[ -e "$state" ]]; then
    backup=$(mktemp "$HOME/stage230-state.previous.XXXXXX")
    cp -p "$state" "$backup"
    echo "Previous state preserved: $backup"
fi
printf 'kit=%q\ntree=%q\nfixture=%q\n' "$kit" "$tree" "$fixture" > "$state"
cd "$HOME"
build_job=$(sbatch --parsable "$kit/build_230.sbatch" "$kit" "$tree")
build_job=${build_job%%;*}
[[ "$build_job" =~ ^[0-9]+$ ]] || { echo 'Unexpected sbatch job ID' >&2; exit 1; }
printf 'build_job=%q\n' "$build_job" >> "$state"
echo "Build job: $build_job"
echo "Log: $HOME/wavetrisk-230-build-${build_job}.log"
