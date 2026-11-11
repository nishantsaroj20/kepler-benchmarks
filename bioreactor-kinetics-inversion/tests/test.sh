#!/bin/bash
set -u
mkdir -p /logs/verifier
cd /tests
pytest test_evaluation.py --ctrf /logs/verifier/ctrf-report.json
TEST_EXIT_CODE=$?
if [ $TEST_EXIT_CODE -eq 0 ]; then
  echo "1" > /logs/verifier/reward.txt
else
  echo "0" > /logs/verifier/reward.txt
fi
exit 0
