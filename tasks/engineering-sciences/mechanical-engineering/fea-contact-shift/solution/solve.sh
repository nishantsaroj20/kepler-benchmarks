#!/bin/bash
set -euo pipefail
mkdir -p /app/solver /app/output
cp /solution/analyze.py /app/solver/analyze.py
chmod 644 /app/solver/analyze.py

# Smoke the practice packs so the oracle leaves a working analyzer in place.
for pack in /app/data/practice/practice_*; do
  out="/app/output/$(basename "$pack").json"
  python3 /app/solver/analyze.py --pack "$pack" --out "$out"
done
