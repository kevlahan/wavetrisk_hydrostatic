#!/bin/sh
set -eu
repo_dir=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
test_dir=$(mktemp -d "${TMPDIR:-/tmp}/wavetrisk-inverse-routes.XXXXXX")
printf 'Inverse route test directory: %s\n' "$test_dir"
for optimization in 0 2; do
  "${FC:-gfortran}" -O"$optimization" -g -cpp -Wall -Wextra -Wno-unused-dummy-argument \
    -fcheck=all -finit-real=snan -ffpe-trap=invalid,zero,overflow \
    -J"$test_dir" -I"$test_dir" \
    "$repo_dir/src/kind.f90" "$repo_dir/src/param_J5.f90" \
    "$repo_dir/src/shared.f90" "$repo_dir/src/patch.f90" \
    "$repo_dir/src/parallel_block_profile.f90" "$repo_dir/src/parallel_block.f90" \
    "$repo_dir/test/parallel_block_inverse/test_routes.f90" -o "$test_dir/routes-O$optimization"
  "$test_dir/routes-O$optimization"
done
