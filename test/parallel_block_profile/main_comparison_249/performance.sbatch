#!/bin/bash
#SBATCH --job-name=rk249-main-direct
#SBATCH --partition=nk
#SBATCH --nodes=1
#SBATCH --ntasks=256
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=04:00:00
#SBATCH --output=rk249-%j.log
set -euo pipefail
if (( $# != 5 )); then
  echo 'Usage: sbatch performance.sbatch PERF_KIT CANDIDATE_KIT MAIN_TREE FRESH_RESULTS J5_FIXTURE' >&2
  exit 2
fi
perf=$1
shift
python3 "$perf/run.py" "$perf" "$@"
