#!/bin/bash
mkdir -p /logs/verifier
chmod 777 /logs/verifier || true
if [ "$(id -u)" -eq 0 ]; then
  chmod 700 /tests/truth || true
fi
cd /tmp
if python3 -m pytest -p no:cacheprovider \
  --ctrf /logs/verifier/ctrf.json \
  --junitxml /logs/verifier/junit.xml \
  /tests/test_outputs.py; then
  echo 1 > /logs/verifier/reward.txt
else
  echo 0 > /logs/verifier/reward.txt
fi
