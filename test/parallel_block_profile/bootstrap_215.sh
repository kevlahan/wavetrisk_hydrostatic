#!/bin/bash
# Run as a script, never source into a login shell.
set -euo pipefail
assets=$(cd "$(dirname "$0")" && pwd)
kit="$assets"
tree="$HOME/wt-stage215"
fixture="$HOME/stage187-fixture-20260911-145307"
state="$HOME/stage215-state.sh"
if [[ -e "$tree" ]]; then
    echo "Source directory already exists: $tree; stopping without modifying it." >&2
    exit 1
fi
python3 - "$assets/stage196-source.tar.gz" "$kit" "$tree" <<'PY'
import hashlib,json,sys,tarfile,tempfile
from pathlib import Path
archive,kit,tree=map(Path,sys.argv[1:])
manifest=json.loads((kit/'stage215-sources.json').read_text())
assert manifest['baseline_commit']=='9bc3a19556dab79b3cba3944cf5df53517389da3'
with tempfile.TemporaryDirectory(dir=tree.parent,prefix='.stage215-source-') as tmp:
    root=Path(tmp)/'source';root.mkdir()
    with tarfile.open(archive) as t:t.extractall(root,filter='data')
    for name,expected in manifest['baseline_sources'].items():
        path=root/name
        data=str(path.readlink()).encode() if path.is_symlink() else path.read_bytes()
        if hashlib.sha256(data).hexdigest()!=expected:
            raise ValueError('Frozen source mismatch: '+name)
    root.rename(tree)
print('PASS: fresh source tree matches the frozen Stage196 inventory')
PY
if [[ -e "$state" ]]; then
    backup=$(mktemp "$HOME/stage215-state.previous.XXXXXX")
    cp -p "$state" "$backup"
    echo "Previous state preserved: $backup"
fi
printf 'kit=%q\ntree=%q\nfixture=%q\n' "$kit" "$tree" "$fixture" > "$state"
cd "$HOME"
build_job=$(sbatch --parsable "$kit/build_215.sbatch" "$kit" "$tree")
printf 'build_job=%q\n' "$build_job" >> "$state"
echo "Build job: $build_job"
echo "Log: $HOME/wavetrisk-215-build-${build_job}.log"
