#!/bin/bash
set -euo pipefail
mkdir -p /logs/verifier
chmod 700 /logs/verifier
echo 0 > /logs/verifier/reward.txt
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
cd /tests
if python3 -m pytest -rA -p no:cacheprovider \
  --ctrf /logs/verifier/ctrf.json /tests/test_outputs.py; then
  echo 1 > /logs/verifier/reward.txt
else
  echo 0 > /logs/verifier/reward.txt
fi
