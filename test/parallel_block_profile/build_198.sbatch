#!/bin/bash
#SBATCH --job-name=wavetrisk-198-build
#SBATCH --partition=nk
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --output=wavetrisk-198-build-%j.log
set -euo pipefail
if [[ $# -ne 2 ]]; then
    echo 'Usage: sbatch build_198.sbatch KIT_DIR FRESH_SOURCE_TREE' >&2
    exit 2
fi
python3 "$1/build_restriction_198.py" "$1" "$2"
