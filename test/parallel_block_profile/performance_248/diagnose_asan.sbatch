#!/bin/bash
#SBATCH --job-name=rk248-asan
#SBATCH --partition=nk
#SBATCH --nodes=1
#SBATCH --ntasks=256
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=01:00:00
#SBATCH --output=rk248-asan-%j.log
# Locate the J4 feature-off startup SIGSEGV with a fresh ASan/UBSan build.
set -euo pipefail
if (( $# != 1 )) || [[ ! "$1" =~ ^[0-9]+$ ]]; then
  echo 'Usage: sbatch diagnose_asan.sbatch FAILED_JOB_ID' >&2
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
[[ -f "$out/builds.json" && -f "$out/placement-16.json" && -f "$out/input/simple.in" ]] || {
  echo "Missing original build/input/placement records in $out" >&2
  exit 1
}
for name in ${!WAVETRISK_@}; do unset "$name"; done
source "$source_tree/runtime.env"
export MV2_ENABLE_AFFINITY=0 OMPI_MCA_hwloc_base_binding_policy=none
export ASAN_OPTIONS=abort_on_error=1:detect_leaks=0:halt_on_error=1
export UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1
unset SLURM_CPUS_PER_TASK SRUN_CPUS_PER_TASK SLURM_EXACT SLURM_EXCLUSIVE SLURM_TRES_PER_TASK SLURM_HINT SLURM_THREADS_PER_CORE
diag="$HOME/rk248-asan-${prior_job}-${SLURM_JOB_ID}"
mkdir "$diag"
printf 'original_job=%s\ndiagnostic_job=%s\nhost=%s\nstack_limit=%s\n' \
  "$prior_job" "$SLURM_JOB_ID" "$(hostname)" "$(ulimit -s)" > "$diag/summary.txt"
set +e
python3 - "$out/builds.json" "$source_tree/source" "$diag/build.log" <<'PY'
import json,subprocess,sys
from pathlib import Path
builds=Path(sys.argv[1]);source=Path(sys.argv[2]);log=Path(sys.argv[3])
command=json.loads(builds.read_text())['builds']['baseline']['command']
command=[('DEBUG=asan' if x=='DEBUG=false' else
          'BUILD_DIR=build-asan-J4' if x.startswith('BUILD_DIR=') else
          'BIN_DIR=bin-asan-J4' if x.startswith('BIN_DIR=') else x) for x in command]
with log.open('w') as stream:
    stream.write(' '.join(command)+'\n');stream.flush()
    raise SystemExit(subprocess.run(command,cwd=source,stdout=stream,stderr=subprocess.STDOUT).returncode)
PY
build_code=$?
set -e
echo "asan_build_exit_code=$build_code" | tee -a "$diag/summary.txt"
if (( build_code == 0 )); then
  directory="$diag/run"
  mkdir "$directory"
  cp -a "$out/input/." "$directory/"
  cp "$source_tree/source/bin-asan-J4/climate" "$directory/climate"
  launcher="$source_tree/source/test/parallel_block_profile/rank_launch_223.py"
  set +e
  (cd "$directory" && srun --nodes=1 --ntasks=16 --cpu-bind=none --kill-on-bad-exit=1 --time=00:10:00 \
    python3 "$launcher" "$out/placement-16.json" "$directory" ./climate simple.in) > "$directory/run.log" 2>&1
  code=$?
  set -e
  echo "asan_run_exit_code=$code" | tee -a "$diag/summary.txt"
  tail -n 100 "$directory/run.log"
else
  tail -n 80 "$diag/build.log"
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
