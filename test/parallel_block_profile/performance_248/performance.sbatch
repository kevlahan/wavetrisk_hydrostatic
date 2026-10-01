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
if (( $# < 4 || $# > 5 )); then
  echo 'Usage: sbatch performance.sbatch CANDIDATE_KIT FRESH_RESULTS J4|J5 mirrored|feature-off [J5_FIXTURE]' >&2
  exit 2
fi
python3 "$(dirname "$0")/run.py" "$@"
