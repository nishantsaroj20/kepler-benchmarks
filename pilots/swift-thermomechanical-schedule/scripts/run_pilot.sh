#!/usr/bin/env bash
# Re-run the feasibility pilot from the pilot source directory.
set -euo pipefail
cd "$(dirname "$0")/../src"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
python3 run_pilot.py
