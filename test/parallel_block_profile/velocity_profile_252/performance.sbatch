#!/bin/bash
#SBATCH --job-name=rk252-velocity-profile
#SBATCH --partition=nk
#SBATCH --nodes=1
#SBATCH --ntasks=256
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=02:00:00
#SBATCH --output=rk252-%j.log
set -euo pipefail
if (( $# != 5 )); then
  echo 'Usage: sbatch performance.sbatch PROFILE_KIT STAGE251_KIT CANDIDATE_SOURCE STAGE251_RESULTS FRESH_RESULTS' >&2
  exit 2
fi
perf=$1
shift
python3 "$perf/run.py" "$perf" "$@"
