#!/bin/bash
set -u
# Pin BLAS/OpenMP threads for reproducible per-run timings.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

mkdir -p /logs/verifier
# Lock verifier log dir: only the root verifier process should write reward/CTRF.
chown -R root:root /logs/verifier 2>/dev/null || true
chmod 700 /logs/verifier 2>/dev/null || true
# Hide sealed hidden truth from non-root agent child processes.
if [ -d /tests/hidden ]; then
  chown -R root:root /tests/hidden 2>/dev/null || true
  chmod -R 700 /tests/hidden 2>/dev/null || true
fi

cd /tests
pytest test_evaluation.py --ctrf /logs/verifier/ctrf.json
TEST_EXIT_CODE=$?
if [ $TEST_EXIT_CODE -eq 0 ]; then
  echo "1" > /logs/verifier/reward.txt
else
  echo "0" > /logs/verifier/reward.txt
fi
# Freeze reward/CTRF against later agent tampering.
chmod 400 /logs/verifier/reward.txt 2>/dev/null || true
chmod 400 /logs/verifier/ctrf.json 2>/dev/null || true
exit 0
