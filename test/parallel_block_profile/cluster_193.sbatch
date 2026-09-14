#!/bin/bash
#SBATCH --job-name=wavetrisk-193
#SBATCH --partition=nk
#SBATCH --nodes=1
#SBATCH --ntasks=83
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --output=wavetrisk-193-%j.log
set -euo pipefail
if [[ $# -ne 4 ]]; then
    echo 'Usage: sbatch cluster_193.sbatch KIT_DIR SOURCE_TREE FIXTURE_DIR NEW_OUTPUT_DIR' >&2
    exit 2
fi
python3 "$1/submit_restriction_193.py" "$1" "$2" "$3" "$4"
