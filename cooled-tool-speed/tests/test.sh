#!/bin/bash
reward=0
mkdir -p /logs/verifier
trap 'printf "%s\n" "$reward" > /logs/verifier/reward.txt' EXIT
cd /tmp
if python -m pytest /tests/test_outputs.py -q -p no:cacheprovider --ctrf /logs/verifier/ctrf.json; then
  reward=1
fi
