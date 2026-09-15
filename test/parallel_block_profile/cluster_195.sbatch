#!/bin/bash
#SBATCH --job-name=wavetrisk-195
#SBATCH --partition=nk
#SBATCH --nodes=1
#SBATCH --ntasks=83
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --output=wavetrisk-195-%j.log
set -euo pipefail
if [[ $# -ne 4 ]]; then
    echo 'Usage: sbatch cluster_195.sbatch KIT_DIR SOURCE_TREE FIXTURE_DIR NEW_OUTPUT_DIR' >&2
    exit 2
fi
python3 "$1/submit_restriction_195.py" "$1" "$2" "$3" "$4"
