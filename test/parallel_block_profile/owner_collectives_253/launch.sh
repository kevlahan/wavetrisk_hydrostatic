#!/bin/bash
# Compare verified Stage251 c5 with a one-source-file optimized c5 build.
set -euo pipefail
if (( $# != 0 )); then
  echo 'Usage: bash launch.sh' >&2
  exit 2
fi
kit=$(cd "$(dirname "$0")" && pwd)
stage251="$HOME/rk251-feature-ladder-kit"
state="$HOME/rk251-job-55256.state"
[[ -f "$state" ]] || { echo "Stage251 submission state absent: $state" >&2; exit 1; }
source "$state"
[[ "$job" == 55256 ]] || { echo 'Wrong Stage251 job' >&2; exit 1; }
prior="$out"
# Stage251 records the extracted candidate kit root, not its source directory.
source_tree="$candidate/source"
command -v sbatch >/dev/null || { echo 'sbatch unavailable; run this on bbserv.' >&2; exit 1; }
python3 "$kit/preflight.py" "$kit" "$stage251" "$source_tree" "$prior"
work=$(mktemp -d "$HOME/rk253-collectives-XXXXXXXX")
results="$work/results"
cd "$HOME"
job=$(sbatch --parsable "$kit/performance.sbatch" "$kit" "$stage251" "$source_tree" "$prior" "$results")
job=${job%%;*}
[[ "$job" =~ ^[0-9]+$ ]] || { echo "Unexpected sbatch response: $job" >&2; exit 1; }
state="$HOME/rk253-job-${job}.state"
printf 'job=%q\nkit=%q\nout=%q\nprior=%q\n' "$job" "$kit" "$results" "$prior" > "$state"
cp "$state" "$HOME/rk253-latest.state"
echo "Submitted job: $job"
echo "Results: $results"
echo "Slurm log: $HOME/rk253-${job}.log"
echo "After the job ends: bash $kit/collect.sh"
