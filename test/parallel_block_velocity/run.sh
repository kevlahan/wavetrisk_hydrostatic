#!/bin/sh
set -eu
repo_dir=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
test_dir=$(mktemp -d "${TMPDIR:-/tmp}/wavetrisk-velocity-kernels.XXXXXX")
printf 'Kernel test build directory: %s\n' "$test_dir"
for optimization in 0 2; do
  "${FC:-gfortran}" -O"$optimization" -g -Wall -Wextra -Werror -Wimplicit-interface \
    -fcheck=all -finit-real=snan -ffpe-trap=invalid,zero,overflow \
    -J"$test_dir" -I"$test_dir" \
    "$repo_dir/src/kind.f90" "$repo_dir/src/parallel_block_velocity.f90" \
    "$repo_dir/test/parallel_block_velocity/test_kernels.f90" \
    -o "$test_dir/kernels-O$optimization"
  "$test_dir/kernels-O$optimization"
done
