#!/bin/bash
# Prepare a fresh source tree, submit once, and record the job and result paths.
set -euo pipefail
perf=$(cd "$(dirname "$0")" && pwd)
level=${1:-J4}
comparison=${2:-feature-off}
if (( $# > 2 )) || [[ "$level" != J4 && "$level" != J5 ]] || [[ "$comparison" != mirrored && "$comparison" != feature-off ]]; then
  echo 'Usage: bash launch.sh [J4|J5] [feature-off|mirrored]' >&2
  exit 2
fi
archive="$HOME/bbserv-candidate4.tar.gz"
[[ -f "$archive" ]] || { echo "Missing qualified source archive: $archive" >&2; exit 1; }
command -v sbatch >/dev/null || { echo 'sbatch is unavailable; run this on bbserv.' >&2; exit 1; }
printf '%s  %s\n' edc647fb8aff057444d0d4d9e93e2756eba1749c5fdcacff332e0404735b8eb0 "$archive" | sha256sum -c -
python3 - "$perf" <<'PY'
import hashlib,json,sys
from pathlib import Path
root=Path(sys.argv[1])
for name,expected in json.loads((root/'files-sha256.json').read_text()).items():
    actual=hashlib.sha256((root/name).read_bytes()).hexdigest()
    if actual!=expected:raise SystemExit('Performance kit file differs: '+name)
print('PASS: performance kit seal')
PY
if [[ "$level" == J5 ]]; then
  fixture="$HOME/stage187-fixture-20260911-145307"
  [[ -f "$fixture/simple.in" ]] || { echo "Missing J5 fixture: $fixture/simple.in" >&2; exit 1; }
fi
work=$(mktemp -d "$HOME/rk248-${level}-${comparison}-XXXXXXXX")
tar -xzf "$archive" -C "$work"
source_tree="$work/bbserv-candidate4"
[[ -f "$source_tree/files-sha256.json" ]] || { echo 'Source extraction failed.' >&2; exit 1; }
results="$work/results"
cd "$HOME"
if [[ "$level" == J5 ]]; then
  job=$(sbatch --parsable "$perf/performance.sbatch" "$perf" "$source_tree" "$results" "$level" "$comparison" "$fixture")
else
  job=$(sbatch --parsable "$perf/performance.sbatch" "$perf" "$source_tree" "$results" "$level" "$comparison")
fi
job=${job%%;*}
[[ "$job" =~ ^[0-9]+$ ]] || { echo "Unexpected sbatch response: $job" >&2; exit 1; }
state="$HOME/rk248-job-${job}.state"
printf 'job=%q\nkit=%q\nout=%q\nsource_tree=%q\nlevel=%q\ncomparison=%q\n' \
  "$job" "$perf" "$results" "$source_tree" "$level" "$comparison" > "$state"
cp "$state" "$HOME/rk248-latest.state"
echo "Submitted job: $job"
echo "Results: $results"
echo "Slurm log: $HOME/rk248-${job}.log"
echo "After the job ends: bash $perf/collect.sh"
