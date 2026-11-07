"""Reference solution: medium-mesh axisymmetric fit, then a speed bisection."""

from __future__ import annotations

import json
from pathlib import Path

from fem_core import predict

MESH = {
    "nr_shaft": 8,
    "nr_sleeve": 8,
    "nz_contact": 16,
    "nz_side": 8,
    "penalty_scale": 5000.0,
}


def main():
    problem = json.loads(Path("/app/data/problem.json").read_text())
    result = predict(problem, MESH)
    destination = Path("/app/output/fit_result.json")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(result) + "\n")


if __name__ == "__main__":
    main()
