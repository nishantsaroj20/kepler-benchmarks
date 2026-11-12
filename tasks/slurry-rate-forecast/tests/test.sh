#!/bin/bash
set -u
mkdir -p /logs/verifier
chmod 700 /logs/verifier
if [ -d /tests ]; then
  chmod 700 /tests
fi
reward=0
write_reward() {
  chmod 700 /logs/verifier
  printf '%s\n' "$reward" > /logs/verifier/reward.txt
  chmod 700 /logs/verifier/reward.txt || true
}
trap write_reward EXIT
cd /tmp
python3 -m pytest -p no:cacheprovider \
  --ctrf /logs/verifier/ctrf.json \
  --junitxml /logs/verifier/junit.xml \
  /tests/test_outputs.py
status=$?
if [ "$status" -eq 0 ]; then
  reward=1
fi
exit 0
