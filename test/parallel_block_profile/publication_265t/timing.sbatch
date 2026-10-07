#!/bin/bash
#SBATCH --job-name=rk265t-tiled-publication
#SBATCH --partition=nk
#SBATCH --nodes=1
#SBATCH --ntasks=256
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=03:00:00
set -euo pipefail
if (( $# != 3 && $# != 4 )); then
  echo 'Usage: sbatch timing.sbatch KIT J5_FIXTURE RESULTS [--repeat]' >&2
  exit 2
fi
kit=$1
python3 "$kit/run.py" "$@"
