#!/usr/bin/env bash
# Always write a literal 0 or 1 reward (swift-notch pattern).
set -uo pipefail
mkdir -p /logs/verifier
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

trap 'if [ ! -f /logs/verifier/reward.txt ]; then printf "0\n" > /logs/verifier/reward.txt.tmp && mv -f /logs/verifier/reward.txt.tmp /logs/verifier/reward.txt; fi' EXIT

printf '0\n' > /logs/verifier/reward.txt.tmp
mv -f /logs/verifier/reward.txt.tmp /logs/verifier/reward.txt

cd /tests
python3 -m pytest -rA -p no:cacheprovider --ctrf /logs/verifier/ctrf.json /tests/test_outputs.py
status=$?

if [ "$status" -eq 0 ]; then
  printf '1\n' > /logs/verifier/reward.txt.tmp
  mv -f /logs/verifier/reward.txt.tmp /logs/verifier/reward.txt
else
  printf '0\n' > /logs/verifier/reward.txt.tmp
  mv -f /logs/verifier/reward.txt.tmp /logs/verifier/reward.txt
fi
exit 0
