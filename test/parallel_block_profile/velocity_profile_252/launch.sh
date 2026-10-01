#!/bin/bash
# Reuse the exact verified Stage251 c4/c5 binaries for matched phase profiling.
set -euo pipefail
if (( $# != 0 )); then
  echo 'Usage: bash launch.sh' >&2
  exit 2
fi
perf=$(cd "$(dirname "$0")" && pwd)
stage251="$HOME/rk251-feature-ladder-kit"
state="$HOME/rk251-job-55256.state"
[[ -f "$state" ]] || { echo "Stage251 submission state absent: $state" >&2; exit 1; }
source "$state"
[[ "$job" == 55256 ]] || { echo 'Wrong Stage251 job' >&2; exit 1; }
prior="$out"
source_tree="$candidate"
command -v sbatch >/dev/null || { echo 'sbatch unavailable; run this on bbserv.' >&2; exit 1; }
python3 "$perf/preflight.py" "$perf" "$stage251" "$source_tree" "$prior"
work=$(mktemp -d "$HOME/rk252-profile-XXXXXXXX")
results="$work/results"
cd "$HOME"
job=$(sbatch --parsable "$perf/performance.sbatch" "$perf" "$stage251" "$source_tree" "$prior" "$results")
job=${job%%;*}
[[ "$job" =~ ^[0-9]+$ ]] || { echo "Unexpected sbatch response: $job" >&2; exit 1; }
state="$HOME/rk252-job-${job}.state"
printf 'job=%q\nkit=%q\nout=%q\nprior=%q\n' "$job" "$perf" "$results" "$prior" > "$state"
cp "$state" "$HOME/rk252-latest.state"
echo "Submitted job: $job"
echo "Results: $results"
echo "Slurm log: $HOME/rk252-${job}.log"
echo "After the job ends: bash $perf/collect.sh"
