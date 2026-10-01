#!/bin/bash
#SBATCH --job-name=rk248-gdb-v5
#SBATCH --partition=nk
#SBATCH --nodes=1
#SBATCH --ntasks=256
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:30:00
#SBATCH --output=rk248-gdb-v5-%j.log
# Capture the next startup fault in the existing optimized, symbol-bearing v5 binary.
set -euo pipefail
if (( $# != 1 )) || [[ ! "$1" =~ ^[0-9]+$ ]]; then
  echo 'Usage: sbatch diagnose_v5_gdb.sbatch FAILED_JOB_ID' >&2
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
[[ "$source_tree" == *'/bbserv-candidate4' ]] || {
  echo 'Expected the v5 corrected source tree.' >&2
  exit 1
}
[[ "${SLURM_JOB_PARTITION:-}" == nk && "$(hostname -s)" == bb13 ]] || {
  echo 'Require the nk partition on bb13.' >&2
  exit 1
}
command -v gdb >/dev/null || { echo 'gdb is unavailable on this node.' >&2; exit 1; }
binary="$source_tree/source/bin-baseline-J4/climate"
[[ -f "$binary" && -f "$out/builds.json" && -f "$out/placement-16.json" && -f "$out/input/simple.in" ]] || {
  echo 'Missing the failed job binary, build record, placement or input.' >&2
  exit 1
}
python3 - "$out/builds.json" "$binary" <<'PY'
import hashlib,json,sys
record=json.load(open(sys.argv[1]))['builds']['baseline']
with open(sys.argv[2],'rb') as stream:
    actual=hashlib.file_digest(stream,'sha256').hexdigest()
if actual!=record['sha256']:
    raise SystemExit('Failed-job baseline binary differs from the build record')
print('Verified failed-job baseline binary:',actual)
PY
for name in ${!WAVETRISK_@}; do unset "$name"; done
source "$source_tree/runtime.env"
export MV2_ENABLE_AFFINITY=0 OMPI_MCA_hwloc_base_binding_policy=none
unset SLURM_CPUS_PER_TASK SRUN_CPUS_PER_TASK SLURM_EXACT SLURM_EXCLUSIVE SLURM_TRES_PER_TASK SLURM_HINT SLURM_THREADS_PER_CORE
diag="$HOME/rk248-gdb-v5-${prior_job}-${SLURM_JOB_ID}"
mkdir "$diag" "$diag/run"
cp -a "$out/input/." "$diag/run/"
cp "$binary" "$diag/run/climate"
{
  echo "original_job=$prior_job"
  echo "debugger_job=$SLURM_JOB_ID"
  echo "host=$(hostname)"
  echo "gdb=$(command -v gdb)"
  echo "binary_sha256=$(sha256sum "$binary" | cut -d' ' -f1)"
  echo "stack_limit=$(ulimit -s)"
} > "$diag/summary.txt"
export DIAG_DIR="$diag"
export RUN_DIR="$diag/run"
export RANK_LAUNCHER="$source_tree/source/test/parallel_block_profile/rank_launch_223.py"
export PLACEMENT="$out/placement-16.json"
set +e
(
  cd "$RUN_DIR"
  srun --nodes=1 --ntasks=16 --cpu-bind=none --kill-on-bad-exit=0 --time=00:12:00 \
    bash -c 'python3 "$RANK_LAUNCHER" "$PLACEMENT" "$RUN_DIR" \
      gdb --nx --quiet --batch \
      -ex "set pagination off" \
      -ex "set confirm off" \
      -ex "set disable-randomization off" \
      -ex "handle SIGSEGV stop print nopass" \
      -ex run \
      -ex "thread apply all bt 40" \
      -ex "info registers rip rsp" \
      --args ./climate simple.in \
      > "$DIAG_DIR/gdb-rank-${SLURM_PROCID}.log" 2>&1'
) > "$diag/srun.log" 2>&1
code=$?
set -e
echo "srun_exit_code=$code" >> "$diag/summary.txt"
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
grep -h -m 2 -E 'Program received signal|runtime error|AddressSanitizer|^#0 ' "$diag"/gdb-rank-*.log | head -n 20 || true
