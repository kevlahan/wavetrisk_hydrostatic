#!/bin/sh
# Use under srun. Instrument selected application ranks, NOT the srun launcher.
# This helper never changes cluster permissions or silently falls back.
set -eu
: "${SLURM_PROCID:?Run this wrapper under srun}"
: "${WAVETRISK_SAMPLE_DIR:?Set a new shared output directory for this run}"
: "${WAVETRISK_SAMPLE_RANKS:=0}"
: "${WAVETRISK_SAMPLE_MODE:=record}"
if [ "$#" -lt 1 ]; then
  echo 'Usage: sh sample_rank.sh /absolute/path/to/climate simple.in' >&2
  exit 2
fi
case ",$WAVETRISK_SAMPLE_RANKS," in
  *",$SLURM_PROCID,"*) ;;
  *) exec "$@" ;;
esac
command -v perf >/dev/null 2>&1 || { echo 'perf is unavailable on this compute node' >&2; exit 2; }
mkdir -p "$WAVETRISK_SAMPLE_DIR"
sample_file="$WAVETRISK_SAMPLE_DIR/rank-$SLURM_PROCID.$WAVETRISK_SAMPLE_MODE"
if [ -e "$sample_file" ]; then
  echo "Refusing to overwrite $sample_file; use a new directory" >&2
  exit 2
fi
case "$WAVETRISK_SAMPLE_MODE" in
  record)
    exec perf record -e cpu-clock:u -F 99 --call-graph dwarf,4096 -o "$sample_file" -- "$@"
    ;;
  stat)
    exec perf stat -x ';' --no-big-num -e task-clock,cycles,instructions,cache-references,cache-misses \
      -o "$sample_file" -- "$@"
    ;;
  *) echo 'WAVETRISK_SAMPLE_MODE must be record or stat' >&2; exit 2 ;;
esac
