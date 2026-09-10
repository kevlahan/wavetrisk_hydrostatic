#!/bin/sh
set -eu
repo_dir=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
test_dir=$(mktemp -d "${TMPDIR:-/tmp}/wavetrisk-detail-profile.XXXXXX")
# Check source text, independently of the compiler's default line limit.
awk 'length($0)>132 && $0 !~ /^[[:space:]]*!/ {
  print FILENAME ":" FNR ": exceeds 132 columns"; bad=1
} END {exit bad}' \
  "$repo_dir/src/parallel_block_profile.f90" "$repo_dir/src/parallel_block.f90" \
  "$repo_dir/src/parallel_block_build.f90" "$repo_dir/src/parallel_block_inverse.f90" \
  "$repo_dir/src/parallel_block_mpi.f90" "$repo_dir/src/comm_mpi.f90" \
  "$repo_dir/src/multi_level.f90" "$repo_dir/src/time_integr.f90" "$repo_dir/src/adapt.f90" \
  "$repo_dir/src/remap.f90" \
  "$repo_dir/test/parallel_block_profile/test_profile.f90" \
  "$repo_dir/test/parallel_block_profile/legacy_profile_io.f90"
for optimization in 0 2; do
  "${FC:-gfortran}" -O"$optimization" -g -Wall -Wextra -Werror -fcheck=all -ffree-line-length-132 \
    -finit-real=snan -ffpe-trap=invalid,zero,overflow -J"$test_dir" -I"$test_dir" \
    "$repo_dir/src/parallel_block_profile.f90" \
    "$repo_dir/test/parallel_block_profile/test_profile.f90" -o "$test_dir/profile-O$optimization"
  "$test_dir/profile-O$optimization"
done
PYTHONDONTWRITEBYTECODE=1 python3 "$repo_dir/test/parallel_block_profile/test_analyze.py"
PYTHONDONTWRITEBYTECODE=1 python3 "$repo_dir/test/parallel_block_profile/test_experiment.py"
PYTHONDONTWRITEBYTECODE=1 python3 "$repo_dir/test/parallel_block_profile/test_legacy_profile.py"
sh -n "$repo_dir/test/parallel_block_profile/sample_rank.sh"
