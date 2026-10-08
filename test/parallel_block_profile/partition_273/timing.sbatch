#!/bin/bash
#SBATCH --job-name=rk273-balanced-partition
#SBATCH --partition=nk
#SBATCH --nodes=1
#SBATCH --ntasks=256
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=03:00:00
set -euo pipefail
if (( $# != 3 )); then
  echo 'Usage: sbatch timing.sbatch KIT J5_FIXTURE RESULTS' >&2
  exit 2
fi
kit=$1
python3 "$kit/run.py" "$@"
