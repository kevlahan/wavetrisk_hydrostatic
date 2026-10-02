#!/bin/bash
#SBATCH --job-name=rk255-owner-depth
#SBATCH --partition=nk
#SBATCH --nodes=1
#SBATCH --ntasks=256
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=02:00:00
#SBATCH --output=rk255-%j.log
set -euo pipefail
if (( $# != 5 )); then
  echo 'Usage: sbatch performance.sbatch KIT STAGE251_KIT CANDIDATE_SOURCE STAGE251_RESULTS FRESH_RESULTS' >&2
  exit 2
fi
perf=$1
shift
python3 "$perf/run.py" "$perf" "$@"
