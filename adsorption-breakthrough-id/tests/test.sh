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
CTRF_PATH=/logs/verifier/ctrf.json
JUNIT_PATH=/logs/verifier/junit.xml
rm -f "$CTRF_PATH" "$JUNIT_PATH"

# Prefer CTRF; also emit JUnit as a countable fallback.
pytest test_evaluation.py --ctrf "$CTRF_PATH" --junitxml "$JUNIT_PATH"
TEST_EXIT_CODE=$?

EXPECTED_TESTS=10
REWARD=0

python3 - <<'PY' "$CTRF_PATH" "$JUNIT_PATH" "$EXPECTED_TESTS" /logs/verifier/reward.txt
import json, sys, xml.etree.ElementTree as ET
from pathlib import Path

ctrf_path, junit_path, expected_s, reward_path = sys.argv[1:5]
expected = int(expected_s)
passed = failed = skipped = None
source = None

ctrf = Path(ctrf_path)
if ctrf.is_file() and ctrf.stat().st_size > 0:
    try:
        data = json.loads(ctrf.read_text())
        # pytest-ctrf variants nest under results.summary or summary.
        summary = None
        if isinstance(data, dict):
            if isinstance(data.get("results"), dict) and isinstance(data["results"].get("summary"), dict):
                summary = data["results"]["summary"]
            elif isinstance(data.get("summary"), dict):
                summary = data["summary"]
        if summary is not None:
            passed = int(summary.get("passed", summary.get("testsPass", 0)) or 0)
            failed = int(summary.get("failed", summary.get("testsFail", 0)) or 0)
            skipped = int(summary.get("skipped", summary.get("testsSkip", 0)) or 0)
            # Some schemas use "other" / pending; treat nonzero as skip-like.
            other = int(summary.get("other", summary.get("pending", 0)) or 0)
            skipped += other
            source = "ctrf"
    except Exception:
        passed = failed = skipped = None

if source is None:
    junit = Path(junit_path)
    if junit.is_file() and junit.stat().st_size > 0:
        try:
            root = ET.parse(junit).getroot()
            # May be <testsuites> or <testsuite>
            suites = []
            if root.tag == "testsuite":
                suites = [root]
            else:
                suites = list(root.iter("testsuite"))
            passed = failed = skipped = 0
            for suite in suites:
                tests = int(suite.attrib.get("tests", 0) or 0)
                fail = int(suite.attrib.get("failures", 0) or 0)
                err = int(suite.attrib.get("errors", 0) or 0)
                skip = int(suite.attrib.get("skipped", 0) or 0)
                failed += fail + err
                skipped += skip
                passed += max(0, tests - fail - err - skip)
            source = "junit"
        except Exception:
            passed = failed = skipped = None

reward = 0
if source is not None and passed is not None:
    if skipped == 0 and failed == 0 and passed == expected:
        reward = 1
    else:
        reward = 0
    print(
        f"reward_source={source} passed={passed} failed={failed} skipped={skipped} "
        f"expected={expected} reward={reward}"
    )
else:
    print("reward_source=unavailable; defaulting reward=0")
    reward = 0

Path(reward_path).write_text(f"{reward}\n")
raise SystemExit(0 if reward == 1 else 1)
PY
REWARD_PARSE_EXIT=$?

if [ $REWARD_PARSE_EXIT -ne 0 ]; then
  echo "0" > /logs/verifier/reward.txt
fi

# Freeze reward/CTRF against later agent tampering.
chmod 400 /logs/verifier/reward.txt 2>/dev/null || true
chmod 400 /logs/verifier/ctrf.json 2>/dev/null || true
chmod 400 /logs/verifier/junit.xml 2>/dev/null || true

# Harbor expects test.sh to exit 0 after writing reward.
exit 0
