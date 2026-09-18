#!/bin/bash
#SBATCH --job-name=wavetrisk-215
#SBATCH --partition=nk
#SBATCH --nodelist=bb13
#SBATCH --nodes=1
#SBATCH --ntasks=256
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=02:00:00
#SBATCH --output=wavetrisk-215-%j.log
set -euo pipefail
if [[ $# -ne 4 ]]; then
    echo 'Usage: sbatch cluster_215.sbatch KIT_DIR SOURCE_TREE FIXTURE_DIR NEW_OUTPUT_DIR' >&2
    exit 2
fi
unset SLURM_CPUS_PER_TASK SLURM_EXACT SLURM_EXCLUSIVE
python3 "$1/submit_restriction_215.py" "$1" "$2" "$3" "$4"
