#!/bin/sh
set -eu
mkdir -p /app /app/output
if [ -f /solution/forecast.py ]; then
  src=/solution/forecast.py
else
  src=$(CDPATH= cd -- "$(dirname "$0")" && pwd)/forecast.py
fi
cp "$src" /app/forecast.py
chmod a+r /app/forecast.py
test -s /app/forecast.py
