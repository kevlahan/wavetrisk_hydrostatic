#!/bin/bash
#SBATCH --job-name=wavetrisk-226-build
#SBATCH --partition=nk
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=02:00:00
#SBATCH --output=wavetrisk-226-build-%j.log
set -euo pipefail
[[ $# -eq 2 ]] || { echo 'Usage: sbatch build_226.sbatch KIT ROOT' >&2; exit 2; }
python3 "$1/build_capability_226.py" "$1" "$2"
