#!/bin/bash
#SBATCH --job-name=wavetrisk-223-build
#SBATCH --partition=nk
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=02:00:00
#SBATCH --output=wavetrisk-223-build-%j.log
set -euo pipefail
[[ $# -eq 2 ]] || { echo 'Usage: sbatch build_223.sbatch KIT ROOT' >&2; exit 2; }
python3 "$1/build_capability_223.py" "$1" "$2"
