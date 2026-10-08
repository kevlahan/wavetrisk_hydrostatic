#!/bin/bash
set -euo pipefail
kit=$(cd "$(dirname "$0")" && pwd)
python3 "$kit/recovery.py" "$kit"
