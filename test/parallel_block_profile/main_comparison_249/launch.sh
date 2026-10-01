#!/bin/bash
# Validate both sealed sources and the fixture, then submit one fresh bb13 job.
set -euo pipefail
if (( $# != 0 )); then
  echo 'Usage: bash launch.sh' >&2
  exit 2
fi
perf=$(cd "$(dirname "$0")" && pwd)
candidate_archive="$HOME/bbserv-candidate5.tar.gz"
main_archive="$HOME/rk249-main-legacy.tar.gz"
fixture="$HOME/stage187-fixture-20260911-145307"
command -v sbatch >/dev/null || { echo 'sbatch unavailable; run this on bbserv.' >&2; exit 1; }
python3 "$perf/preflight.py" "$perf" "$candidate_archive" "$main_archive" "$fixture"
work=$(mktemp -d "$HOME/rk249-main-XXXXXXXX")
tar -xzf "$candidate_archive" -C "$work"
tar -xzf "$main_archive" -C "$work"
candidate="$work/bbserv-candidate5"
main="$work/rk249-main-legacy"
[[ -f "$candidate/files-sha256.json" && -f "$main/files-sha256.json" ]] || {
  echo 'Source extraction failed.' >&2
  exit 1
}
results="$work/results"
cd "$HOME"
job=$(sbatch --parsable "$perf/performance.sbatch" "$perf" "$candidate" "$main" "$results" "$fixture")
job=${job%%;*}
[[ "$job" =~ ^[0-9]+$ ]] || { echo "Unexpected sbatch response: $job" >&2; exit 1; }
state="$HOME/rk249-job-${job}.state"
printf 'job=%q\nkit=%q\nout=%q\ncandidate=%q\nmain=%q\nfixture=%q\n' \
  "$job" "$perf" "$results" "$candidate" "$main" "$fixture" > "$state"
cp "$state" "$HOME/rk249-latest.state"
echo "Submitted job: $job"
echo "Results: $results"
echo "Slurm log: $HOME/rk249-${job}.log"
echo "After the job ends: bash $perf/collect.sh"
