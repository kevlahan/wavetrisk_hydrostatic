#!/bin/bash
#SBATCH --job-name=wavetrisk-223
#SBATCH --partition=nk
#SBATCH --nodes=1
#SBATCH --ntasks=256
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=04:00:00
#SBATCH --output=wavetrisk-223-%j.log
set -euo pipefail
[[ $# -eq 4 ]] || { echo 'Usage: sbatch cluster_223.sbatch KIT ROOT J5_FIXTURE OUTPUT' >&2; exit 2; }
# The driver records and clears inherited step restrictions before srun.
python3 "$1/cluster_223.py" "$1" "$2" "$3" "$4"
