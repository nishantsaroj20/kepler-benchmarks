"""Write the fine-mesh truth and the mesh-convergence record."""

import json
import sys
from pathlib import Path

import fem_core

ROOT = Path(__file__).resolve().parents[1]
SPEC = json.loads((ROOT / "environment" / "data" / "problem.json").read_text())
TRUTH_PATH = ROOT / "tests" / "truth.json"
EVIDENCE_PATH = ROOT / "authoring" / "evidence" / "convergence.json"
TEMPERATURE_TOL = 0.03
GROWTH_TOL = 0.03
SPEED_TOL = 0.02


def _max_relative(got, expected):
    return max(abs(a - b) / max(abs(b), 1.0) for a, b in zip(got, expected))


def main():
    annulus_error = fem_core.annulus_temperature_error_c()
    if annulus_error > 0.5:
        sys.exit(f"annulus check failed: {annulus_error} C")
    fine = fem_core.evaluate(SPEC, "fine")
    medium = fem_core.evaluate(SPEC, "medium")
    constant = fem_core.evaluate(SPEC, "fine", constant_conductivity=True)
    solid = fem_core.evaluate(SPEC, "fine", with_hole=False)
    temperature_spread = _max_relative(medium["rake_temperature_c"], fine["rake_temperature_c"])
    growth_spread = abs(medium["edge_growth_um"] - fine["edge_growth_um"]) / abs(fine["edge_growth_um"])
    speed_spread = abs(medium["limit_speed_m_s"] - fine["limit_speed_m_s"]) / abs(fine["limit_speed_m_s"])
    if temperature_spread > TEMPERATURE_TOL or growth_spread > GROWTH_TOL or speed_spread > SPEED_TOL:
        sys.exit(
            f"medium mesh is outside the acceptance band: {temperature_spread}, {growth_spread}, {speed_spread}"
        )
    constant_misses = sum(
        abs(a - b) / max(abs(b), 1.0) > TEMPERATURE_TOL
        for a, b in zip(constant["rake_temperature_c"], fine["rake_temperature_c"])
    )
    if constant_misses < 4:
        sys.exit("constant conductivity is inside the station band")
    truth = {
        "rake_temperature_c": fine["rake_temperature_c"],
        "edge_growth_um": fine["edge_growth_um"],
        "limit_speed_m_s": fine["limit_speed_m_s"],
        "temperature_rel_tol": TEMPERATURE_TOL,
        "growth_rel_tol": GROWTH_TOL,
        "speed_rel_tol": SPEED_TOL,
        "constant_conductivity_rake_temperature_c": constant["rake_temperature_c"],
    }
    evidence = {
        "annulus_max_abs_error_c": annulus_error,
        "medium_vs_fine_rake_temperature_rel": temperature_spread,
        "medium_vs_fine_edge_growth_rel": growth_spread,
        "medium_vs_fine_limit_speed_rel": speed_spread,
        "constant_conductivity_rake_temperature_rel": _max_relative(
            constant["rake_temperature_c"], fine["rake_temperature_c"]
        ),
        "constant_conductivity_limit_speed_rel": abs(constant["limit_speed_m_s"] - fine["limit_speed_m_s"])
        / abs(fine["limit_speed_m_s"]),
        "solid_tool_rake_temperature_rel": _max_relative(solid["rake_temperature_c"], fine["rake_temperature_c"]),
        "solid_tool_edge_growth_rel": abs(solid["edge_growth_um"] - fine["edge_growth_um"]) / abs(fine["edge_growth_um"]),
        "fine_node_count": fine["node_count"],
        "medium_node_count": medium["node_count"],
        "fine_rake_temperature_c": fine["rake_temperature_c"],
        "fine_edge_growth_um": fine["edge_growth_um"],
        "fine_limit_speed_m_s": fine["limit_speed_m_s"],
    }
    TRUTH_PATH.parent.mkdir(parents=True, exist_ok=True)
    EVIDENCE_PATH.parent.mkdir(parents=True, exist_ok=True)
    TRUTH_PATH.write_text(json.dumps(truth, indent=2) + "\n")
    EVIDENCE_PATH.write_text(json.dumps(evidence, indent=2) + "\n")
    print(json.dumps(evidence, indent=2))


if __name__ == "__main__":
    main()
