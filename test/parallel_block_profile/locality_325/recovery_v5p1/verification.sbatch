#!/bin/bash
#SBATCH --job-name=rk325-native-v5p1
#SBATCH --partition=nk
#SBATCH --nodelist=bb13
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:30:00
set -euo pipefail
if (( $# != 2 )); then
  echo 'Usage: sbatch verification.sbatch KIT NEW_OUTPUT' >&2
  exit 2
fi
ulimit -S -v 16777216
ulimit -H -v 16777216
if [[ "$(ulimit -S -v)" != "16777216" || "$(ulimit -H -v)" != "16777216" ]]; then
  echo 'Require exact 16-GiB soft and hard virtual-memory limits' >&2
  exit 2
fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
exec srun --ntasks=1 --cpus-per-task=16 --cpu-bind=cores python3 -B "$1/native_runner325.py" supervise "$1" "$2"
