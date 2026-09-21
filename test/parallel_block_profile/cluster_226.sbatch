#!/bin/bash
#SBATCH --job-name=wavetrisk-226
#SBATCH --partition=nk
#SBATCH --nodes=1
#SBATCH --ntasks=256
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=04:00:00
#SBATCH --output=wavetrisk-226-%j.log
set -euo pipefail
[[ $# -eq 4 ]] || { echo 'Usage: sbatch cluster_226.sbatch KIT ROOT J5_FIXTURE OUTPUT' >&2; exit 2; }
# The driver records and clears inherited step restrictions before srun.
python3 "$1/cluster_226.py" "$1" "$2" "$3" "$4"
