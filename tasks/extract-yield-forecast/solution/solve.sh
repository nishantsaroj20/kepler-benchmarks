#!/bin/sh
set -eu
mkdir -p /app/output
cd "$(dirname "$0")"
python3 ./forecast.py --data /app/data --output /app/output/predictions.json
test -s /app/output/predictions.json
