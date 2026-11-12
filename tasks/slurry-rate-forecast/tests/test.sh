#!/bin/bash
set -u
mkdir -p /logs/verifier
on_exit() {
  if [ ! -f /logs/verifier/reward.txt ]; then
    printf '0\n' > /logs/verifier/reward.txt
  fi
}
trap on_exit EXIT
cd /tmp
python3 -m pytest -p no:cacheprovider \
  --ctrf /logs/verifier/ctrf.json \
  --junitxml /logs/verifier/junit.xml \
  /tests/test_outputs.py
status=$?
if [ "$status" -eq 0 ]; then
  printf '1\n' > /logs/verifier/reward.txt
else
  printf '0\n' > /logs/verifier/reward.txt
fi
exit 0
