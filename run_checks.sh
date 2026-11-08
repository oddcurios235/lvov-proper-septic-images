#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p results
for script in verify_septic_trace.py verify_rep_filtration.py verify_deep_minors.py verify_deep_jacobians.py verify_deep_m3.py; do
    echo "[RUN] $script"
    PYTHONOPTIMIZE=0 python3 "code/$script"
done
