#!/bin/bash
#SBATCH --job-name=rk325-verification-recovery-v4
#SBATCH --partition=nk
#SBATCH --nodelist=bb13
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=0
#SBATCH --time=00:05:00
#SBATCH --signal=B:TERM@15
set -euo pipefail
if (( $# != 2 )); then
  echo 'Usage: sbatch verification.sbatch KIT NEW_OUTPUT' >&2
  exit 2
fi
ulimit -v 16777216
if [[ "$(ulimit -v)" != "16777216" ]]; then
  echo "16-GiB virtual-memory limit was not installed" >&2
  exit 2
fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
exec python3 -u "$1/recover_parallel325.py" "$@"
