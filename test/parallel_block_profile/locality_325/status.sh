#!/bin/bash
set -euo pipefail
kit_dir="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
python3 "$kit_dir/control325.py" status "$kit_dir"
