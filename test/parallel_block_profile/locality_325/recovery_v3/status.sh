#!/bin/bash
set -euo pipefail
kit_dir="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
python3 "$kit_dir/control_recovery325.py" status "$kit_dir"
