#!/bin/bash
#SBATCH --job-name=wavetrisk-198
#SBATCH --partition=nk
#SBATCH --nodes=1
#SBATCH --ntasks=83
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --output=wavetrisk-198-%j.log
set -euo pipefail
if [[ $# -ne 4 ]]; then
    echo 'Usage: sbatch cluster_198.sbatch KIT_DIR SOURCE_TREE FIXTURE_DIR NEW_OUTPUT_DIR' >&2
    exit 2
fi
python3 "$1/submit_restriction_198.py" "$1" "$2" "$3" "$4"
