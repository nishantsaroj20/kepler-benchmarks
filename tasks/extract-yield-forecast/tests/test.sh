#!/bin/bash
mkdir -p /logs/verifier
if [ ! -f /logs/verifier/reward.txt ]; then
  echo 0 > /logs/verifier/reward.txt
fi
cd /tmp
if python3 -m pytest -p no:cacheprovider \
  --ctrf /logs/verifier/ctrf.json \
  /tests/test_outputs.py; then
  echo 1 > /logs/verifier/reward.txt
else
  echo 0 > /logs/verifier/reward.txt
fi
