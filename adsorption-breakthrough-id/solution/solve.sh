#!/bin/bash
set -euo pipefail
# Resolve beside this script — Harbor may mount solution outside /workspace.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p /workspace/output
cp "${SCRIPT_DIR}/column_solver.py" /workspace/output/column_solver.py
chmod +x /workspace/output/column_solver.py
# Smoke-test on a visible sample so a broken reference fails during oracle.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
python3 /workspace/output/column_solver.py \
  --case /workspace/data/sample_cases/V1_base/case.json \
  --out /workspace/output/smoke_V1_result.json
