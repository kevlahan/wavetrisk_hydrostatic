#!/bin/bash
#SBATCH --job-name=rk325-locality-qualification-v1
#SBATCH --partition=nk
#SBATCH --nodelist=bb13
#SBATCH --nodes=1
#SBATCH --ntasks=256
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:45:00
set -euo pipefail
if (( $# != 2 )); then
  echo 'Usage: sbatch qualification.sbatch KIT NEW_RESULTS' >&2
  exit 2
fi
python3 "$1/run325.py" "$@"
