"""Write the fine-step truth and the coupling check."""

import json
import sys
from pathlib import Path

import bus_core

ROOT = Path(__file__).resolve().parents[1]
SPEC = json.loads((ROOT / "environment" / "data" / "problem.json").read_text())
TRUTH_PATH = ROOT / "tests" / "truth.json"
EVIDENCE_PATH = ROOT / "authoring" / "evidence" / "convergence.json"
CROSS_TOL = 0.03
CLOCK_TOL = 0.02
FINE = (48, 1e-13)
MEDIUM = (24, 2e-13)


def _max_relative(got, expected):
    return max(abs(a - b) / max(abs(b), 1.0) for a, b in zip(got, expected))


def main():
    fine = bus_core.evaluate(SPEC, *FINE)
    medium = bus_core.evaluate(SPEC, *MEDIUM)
    uncoupled = bus_core.evaluate(SPEC, *FINE, coupling_scale=0.0)
    cross_spread = _max_relative(medium["far_end_cross_ps"], fine["far_end_cross_ps"])
    clock_spread = abs(medium["limit_clock_mhz"] - fine["limit_clock_mhz"]) / abs(fine["limit_clock_mhz"])
    if cross_spread > CROSS_TOL or clock_spread > CLOCK_TOL:
        sys.exit(f"medium step is outside the band: {cross_spread}, {clock_spread}")
    misses = sum(
        abs(a - b) / max(abs(b), 1.0) > CROSS_TOL
        for a, b in zip(uncoupled["far_end_cross_ps"], fine["far_end_cross_ps"])
    )
    if misses < 4:
        sys.exit("uncoupled delays are inside the station band")
    truth = {
        "far_end_cross_ps": fine["far_end_cross_ps"],
        "limit_clock_mhz": fine["limit_clock_mhz"],
        "cross_rel_tol": CROSS_TOL,
        "clock_rel_tol": CLOCK_TOL,
        "uncoupled_far_end_cross_ps": uncoupled["far_end_cross_ps"],
    }
    evidence = {
        "medium_vs_fine_cross_rel": cross_spread,
        "medium_vs_fine_clock_rel": clock_spread,
        "uncoupled_cross_rel": _max_relative(uncoupled["far_end_cross_ps"], fine["far_end_cross_ps"]),
        "uncoupled_clock_rel": abs(uncoupled["limit_clock_mhz"] - fine["limit_clock_mhz"]) / abs(fine["limit_clock_mhz"]),
        "fine_far_end_cross_ps": fine["far_end_cross_ps"],
        "fine_limit_clock_mhz": fine["limit_clock_mhz"],
        "fine_segments": FINE[0],
        "fine_dt_ps": FINE[1] * 1e12,
        "medium_segments": MEDIUM[0],
        "medium_dt_ps": MEDIUM[1] * 1e12,
    }
    TRUTH_PATH.parent.mkdir(parents=True, exist_ok=True)
    EVIDENCE_PATH.parent.mkdir(parents=True, exist_ok=True)
    TRUTH_PATH.write_text(json.dumps(truth, indent=2) + "\n")
    EVIDENCE_PATH.write_text(json.dumps(evidence, indent=2) + "\n")
    print(json.dumps(evidence, indent=2))


if __name__ == "__main__":
    main()
