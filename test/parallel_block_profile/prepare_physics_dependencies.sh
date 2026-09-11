#!/bin/bash
# Generate the physics package's ignored module-dependency include, without
# invoking its clean target (which also reformats source files).
set -euo pipefail
if [[ $# -eq 0 ]]; then
    echo "Usage: $0 SOURCE_TREE [SOURCE_TREE ...]" >&2
    exit 2
fi
for tree in "$@"; do
    (
        cd "$tree/src/physics/simple_physics/phyparam"
        scratch=$(mktemp -d "${TMPDIR:-/tmp}/wavetrisk-physics-deps.XXXXXX")
        trap 'rm -rf "$scratch"' EXIT
        for source in physics/*.F90; do
            unit=${source##*/}
            unit=${unit%.F90}
            awk -v unit="$unit" '{ print unit, NR, $0 }' "$source"
        done > "$scratch/all_code"
        bash bash/makedeps.sh "$scratch/all_code" > "$scratch/Makefile.inc"
        grep -q '^obj/iniphyparam_mod.o : obj/read_param_mod.o$' "$scratch/Makefile.inc"
        cp "$scratch/Makefile.inc" Makefile.inc
        echo "Prepared physics module dependencies in $PWD/Makefile.inc"
    )
done
