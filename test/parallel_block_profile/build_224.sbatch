#!/bin/bash
#SBATCH --job-name=wavetrisk-224-build
#SBATCH --partition=nk
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=02:00:00
#SBATCH --output=wavetrisk-224-build-%j.log
set -euo pipefail
[[ $# -eq 2 ]] || { echo 'Usage: sbatch build_224.sbatch KIT ROOT' >&2; exit 2; }
python3 "$1/build_capability_224.py" "$1" "$2"
