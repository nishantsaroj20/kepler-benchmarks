"""Reference solution: medium fitted mesh for the cooled tool tip."""

import json
from pathlib import Path

import fem_core

SPEC = json.loads(Path("/app/data/problem.json").read_text())
RESULT = Path("/app/output/tool_result.json")


def main():
    answer = fem_core.evaluate(SPEC, "medium")
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "rake_temperature_c": answer["rake_temperature_c"],
        "edge_growth_um": answer["edge_growth_um"],
        "limit_speed_m_s": answer["limit_speed_m_s"],
    }
    RESULT.write_text(json.dumps(payload, indent=2) + "\n")


if __name__ == "__main__":
    main()
