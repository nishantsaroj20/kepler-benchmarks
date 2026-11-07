#!/usr/bin/env bash
# Each write site contains a literal binary reward; replace stale rewards atomically.
set -uo pipefail
mkdir -p /logs/verifier || exit 1
write_failure() {
  printf '0\n' > /logs/verifier/reward.txt.tmp &&
    mv -f /logs/verifier/reward.txt.tmp /logs/verifier/reward.txt
}
trap write_failure EXIT
write_failure || exit 1
if pytest /tests/test_outputs.py -q -p no:cacheprovider --ctrf /logs/verifier/ctrf.json --junitxml=/logs/verifier/junit.xml; then
  printf '1\n' > /logs/verifier/reward.txt.tmp &&
    mv -f /logs/verifier/reward.txt.tmp /logs/verifier/reward.txt || exit 1
  # Disable the fail-closed trap only after successful publication of literal 1.
  trap - EXIT
fi
exit 0
