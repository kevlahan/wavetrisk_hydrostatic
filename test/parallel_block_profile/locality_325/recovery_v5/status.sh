#!/bin/bash
set -euo pipefail
kit_dir="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
exec python3 -B "$kit_dir/control_native325.py" status "$kit_dir"
