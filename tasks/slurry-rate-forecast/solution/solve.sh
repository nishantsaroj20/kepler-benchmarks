#!/bin/bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cp "$HERE/forecast.py" /app/forecast.py
chmod 644 /app/forecast.py
