#!/bin/bash
mkdir -p /logs/verifier /app
chmod 777 /logs/verifier || true
if [ "$(id -u)" -eq 0 ]; then
  chmod a+rx /app || true
  if [ -f /app/forecast.py ]; then
    chmod a+r /app/forecast.py || true
  fi
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
