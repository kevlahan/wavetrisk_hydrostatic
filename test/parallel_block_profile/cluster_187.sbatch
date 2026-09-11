#!/bin/bash
#SBATCH --job-name=wavetrisk-187
#SBATCH --partition=bb
#SBATCH --ntasks=83
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --output=wavetrisk-187-%j.log
set -euo pipefail

if [[ $# -ne 7 ]]; then
    echo 'Usage: sbatch cluster_187.sbatch KIT_DIR BASELINE_OPT CANDIDATE_OPT CANDIDATE_CHECKED FIXTURE_DIR NEW_OUTPUT_DIR BUILD_NOTES' >&2
    exit 2
fi
# Load the site's usual GNU Fortran/MPI/NetCDF/Python environment before sbatch,
# or add the corresponding module commands here. No compiles occur in this job.
python3 "$1/cluster_187.py" \
    --baseline "$2" --candidate "$3" --checked "$4" \
    --fixture "$5" --out "$6" --build-notes "$7" --expected-checkpoints 4
