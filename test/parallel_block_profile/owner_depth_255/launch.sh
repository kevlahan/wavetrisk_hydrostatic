#!/bin/bash
# Compare verified Stage251 c5 with depth-indexed owner execution.
set -euo pipefail
if (( $# != 0 )); then
  echo 'Usage: bash launch.sh' >&2
  exit 2
fi
# The Stage251 state sets "kit"; keep this kit path under a distinct name.
perf=$(cd "$(dirname "$0")" && pwd)
stage251="$HOME/rk251-feature-ladder-kit"
state="$HOME/rk251-job-55256.state"
[[ -f "$state" ]] || { echo "Stage251 submission state absent: $state" >&2; exit 1; }
source "$state"
[[ "$job" == 55256 ]] || { echo 'Wrong Stage251 job' >&2; exit 1; }
prior="$out"
# Stage251 records the extracted candidate kit root, not its source directory.
source_tree="$candidate/source"
command -v sbatch >/dev/null || { echo 'sbatch unavailable; run this on bbserv.' >&2; exit 1; }
python3 "$perf/preflight.py" "$perf" "$stage251" "$source_tree" "$prior"
work=$(mktemp -d "$HOME/rk255-depth-XXXXXXXX")
results="$work/results"
cd "$HOME"
job=$(sbatch --parsable "$perf/performance.sbatch" "$perf" "$stage251" "$source_tree" "$prior" "$results")
job=${job%%;*}
[[ "$job" =~ ^[0-9]+$ ]] || { echo "Unexpected sbatch response: $job" >&2; exit 1; }
state="$HOME/rk255-job-${job}.state"
printf 'job=%q\nkit=%q\nout=%q\nprior=%q\n' "$job" "$perf" "$results" "$prior" > "$state"
cp "$state" "$HOME/rk255-latest.state"
echo "Submitted job: $job"
echo "Results: $results"
echo "Slurm log: $HOME/rk255-${job}.log"
echo "After the job ends: bash $perf/collect.sh"
