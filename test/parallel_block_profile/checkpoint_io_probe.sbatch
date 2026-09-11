#!/bin/bash
#SBATCH --job-name=checkpoint-io-probe
#SBATCH --partition=bb
#SBATCH --ntasks=83
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --time=00:05:00
#SBATCH --output=checkpoint-io-probe-%j.log
set -euo pipefail
if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo 'Usage: sbatch checkpoint_io_probe.sbatch KIT_DIR [OLD_PROBE_DIRECTORY]' >&2
    exit 2
fi
kit_dir="$1"
probe_dir="$HOME/checkpoint-io-probe-$SLURM_JOB_ID"
mkdir "$probe_dir"
cp "$kit_dir/checkpoint_io_probe.f90" "$probe_dir/probe.f90"
cp "$kit_dir/verify_checkpoint_io_probe.py" "$probe_dir/verify.py"
cd "$probe_dir"
mpif90 --version
srun --ntasks=2 --ntasks-per-node=1 bash -c '
  hostname
  if command -v findmnt >/dev/null; then
    findmnt -T "$1" -o TARGET,SOURCE,FSTYPE,OPTIONS
  else
    df -T "$1"
  fi
' _ "$probe_dir"
if [[ $# -eq 2 ]]; then
    echo 'Rechecking the previous probe files before running the new probe:'
    python3 verify.py "$2" --ranks 83
fi
mpif90 -O2 -g -fcheck=all -ffree-line-length-132 probe.f90 -o probe
export OMP_NUM_THREADS=1
srun --ntasks=83 --cpus-per-task=1 --cpu-bind=cores --kill-on-bad-exit=1 ./probe
python3 verify.py "$probe_dir" --ranks 83
