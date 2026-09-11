#!/bin/bash
#SBATCH --job-name=wavetrisk-188
#SBATCH --partition=bb
#SBATCH --ntasks=83
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --output=wavetrisk-188-%j.log
set -euo pipefail
if [[ $# -ne 4 ]]; then
    echo 'Usage: sbatch cluster_188.sbatch KIT_DIR SOURCE_TREE FIXTURE_DIR NEW_OUTPUT_DIR' >&2
    exit 2
fi
python3 "$1/submit_arithmetic_188.py" "$1" "$2" "$3" "$4"
