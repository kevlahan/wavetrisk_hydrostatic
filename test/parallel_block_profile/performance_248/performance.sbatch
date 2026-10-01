#!/bin/bash
#SBATCH --job-name=rk248-performance
#SBATCH --partition=nk
#SBATCH --nodes=1
#SBATCH --ntasks=256
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=04:00:00
#SBATCH --output=rk248-%j.log
set -euo pipefail
if (( $# < 5 || $# > 6 )); then
  echo 'Usage: sbatch performance.sbatch PERF_KIT CANDIDATE_KIT FRESH_RESULTS J4|J5 mirrored|feature-off [J5_FIXTURE]' >&2
  exit 2
fi
perf=$1
shift
python3 "$perf/run.py" "$@"
