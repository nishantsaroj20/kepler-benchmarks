"""Reference solution for the coupled on-chip bus."""

import json
from pathlib import Path

import bus_core

SPEC = json.loads(Path("/app/data/problem.json").read_text())
RESULT = Path("/app/output/bus_timing.json")


def main():
    answer = bus_core.evaluate(SPEC, 24, 2e-13)
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "far_end_cross_ps": answer["far_end_cross_ps"],
        "limit_clock_mhz": answer["limit_clock_mhz"],
    }
    RESULT.write_text(json.dumps(payload, indent=2) + "\n")


if __name__ == "__main__":
    main()
