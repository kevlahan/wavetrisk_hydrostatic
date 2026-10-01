#!/bin/bash
#SBATCH --job-name=rk248-diagnose
#SBATCH --partition=nk
#SBATCH --nodes=1
#SBATCH --ntasks=256
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=01:00:00
#SBATCH --output=rk248-diagnose-%j.log
# Reproduce the failed J4 feature-off baseline, then use a checked build if needed.
set -euo pipefail
if (( $# != 1 )) || [[ ! "$1" =~ ^[0-9]+$ ]]; then
  echo 'Usage: sbatch diagnose.sbatch FAILED_JOB_ID' >&2
  exit 2
fi
prior_job=$1
state="$HOME/rk248-job-${prior_job}.state"
[[ -f "$state" ]] || { echo "Missing submission state: $state" >&2; exit 1; }
source "$state"
[[ "$job" == "$prior_job" && "$level" == J4 && "$comparison" == feature-off ]] || {
  echo 'This diagnostic requires the failed J4 feature-off job state.' >&2
  exit 1
}
[[ -f "$out/results.json" && -f "$out/placement-16.json" && -f "$out/input/simple.in" ]] || {
  echo "Missing original result/input/placement files in $out" >&2
  exit 1
}
source "$source_tree/runtime.env"
export MV2_ENABLE_AFFINITY=0 OMPI_MCA_hwloc_base_binding_policy=none
unset SLURM_CPUS_PER_TASK SRUN_CPUS_PER_TASK SLURM_EXACT SLURM_EXCLUSIVE SLURM_TRES_PER_TASK SLURM_HINT SLURM_THREADS_PER_CORE
diag="$HOME/rk248-diagnostic-${prior_job}-${SLURM_JOB_ID}"
mkdir "$diag"
{
  echo "original_job=$prior_job"
  echo "diagnostic_job=$SLURM_JOB_ID"
  echo "host=$(hostname)"
  echo "stack_limit_before=$(ulimit -s)"
} > "$diag/summary.txt"
stack_raised=0
if ulimit -s unlimited 2> "$diag/stack-limit-error.log"; then
  stack_raised=1
else
  echo 'Could not raise the stack limit; proceeding directly to checked build.' | tee -a "$diag/summary.txt"
fi
echo "stack_limit_after=$(ulimit -s)" >> "$diag/summary.txt"
launcher="$source_tree/source/test/parallel_block_profile/rank_launch_223.py"
placement="$out/placement-16.json"
run_solver() {
  local name=$1 binary=$2 directory="$diag/$1" code
  mkdir "$directory"
  cp -a "$out/input/." "$directory/"
  cp "$binary" "$directory/climate"
  set +e
  (cd "$directory" && srun --nodes=1 --ntasks=16 --cpu-bind=none --kill-on-bad-exit=1 --time=00:10:00 \
    python3 "$launcher" "$placement" "$directory" ./climate simple.in) > "$directory/run.log" 2>&1
  code=$?
  set -e
  echo "$name exit_code=$code" | tee -a "$diag/summary.txt"
  tail -n 35 "$directory/run.log"
  return "$code"
}
if (( stack_raised )) && run_solver optimized "$source_tree/source/bin-baseline-J4/climate"; then
  echo 'Optimized binary passed with unlimited stack; inspect full numerical output before any timing claim.' | tee -a "$diag/summary.txt"
else
  echo 'Building checked baseline after the optimized retry did not pass or stack could not be raised.' | tee -a "$diag/summary.txt"
  set +e
  python3 - "$out/builds.json" "$source_tree/source" "$diag/build-check.log" <<'PY'
import json,subprocess,sys
from pathlib import Path
builds=Path(sys.argv[1]);source=Path(sys.argv[2]);log=Path(sys.argv[3])
command=json.loads(builds.read_text())['builds']['baseline']['command']
command=[('DEBUG=check' if x=='DEBUG=false' else
          'BUILD_DIR=build-diagnostic-J4' if x.startswith('BUILD_DIR=') else
          'BIN_DIR=bin-diagnostic-J4' if x.startswith('BIN_DIR=') else x) for x in command]
with log.open('w') as stream:
    stream.write(' '.join(command)+'\n');stream.flush()
    raise SystemExit(subprocess.run(command,cwd=source,stdout=stream,stderr=subprocess.STDOUT).returncode)
PY
  build_code=$?
  set -e
  echo "checked_build_exit_code=$build_code" | tee -a "$diag/summary.txt"
  if (( build_code == 0 )); then
    run_solver checked "$source_tree/source/bin-diagnostic-J4/climate" || true
  else
    tail -n 40 "$diag/build-check.log"
  fi
fi
python3 - "$diag" <<'PY'
import sys,tarfile
from pathlib import Path
root=Path(sys.argv[1]);archive=root.with_suffix('.tar.gz')
with tarfile.open(archive,'w:gz') as out:
    for path in sorted(root.rglob('*')):
        if path.is_file() and (path.suffix in ('.log','.json','.resources','.txt') or path.name.endswith('_log')):
            out.add(path,arcname=str(path.relative_to(root.parent)),recursive=False)
print('Diagnostic review:',archive)
PY
echo "Summary: $diag/summary.txt"
