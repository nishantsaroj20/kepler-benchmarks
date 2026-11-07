"""Build the sealed truth file from a fine mesh.

Run from the task directory:

    python3 authoring/generate_truth.py

The reference solution is not imported. A medium mesh is solved afterwards
only to record the spread used to set the acceptance bands.
"""

from __future__ import annotations

import json
from pathlib import Path

from fem_core import contact_solution, predict, station_pressures

ROOT = Path(__file__).resolve().parents[1]
PROBLEM_PATH = ROOT / "environment" / "data" / "problem.json"
TRUTH_PATH = ROOT / "tests" / "truth.json"
EVIDENCE_PATH = ROOT / "authoring" / "evidence" / "convergence.json"

FINE = {"nr_shaft": 16, "nr_sleeve": 16, "nz_contact": 32, "nz_side": 16, "penalty_scale": 5000.0}
MEDIUM = {"nr_shaft": 8, "nr_sleeve": 8, "nz_contact": 16, "nz_side": 8, "penalty_scale": 5000.0}


def plane_stress_pressure(problem):
    """Uniform long-cylinder pressure. This is the baseline the short hub must beat."""
    inner = problem["shaft_inner_radius_m"]
    radius = problem["shaft_outer_radius_m"]
    outer = problem["sleeve_outer_radius_m"]
    delta = 0.5 * problem["diametral_interference_m"]
    shaft = problem["shaft"]
    sleeve = problem["sleeve"]

    def compliance(young, poisson, r_in, r_out, loaded_inner):
        if loaded_inner:
            hoop = (r_out**2 + r_in**2) / (r_out**2 - r_in**2)
            radial = -1.0
            r = r_in
        else:
            hoop = -(r_out**2 + r_in**2) / (r_out**2 - r_in**2)
            radial = -1.0
            r = r_out
        return r / young * (hoop - poisson * radial)

    sleeve_u = compliance(
        sleeve["young_modulus_pa"], sleeve["poisson_ratio"], radius, outer, True
    )
    shaft_u = compliance(
        shaft["young_modulus_pa"], shaft["poisson_ratio"], inner, radius, False
    )
    return delta / (sleeve_u - shaft_u)


def plane_strain_pressure(young, poisson, inner, radius, outer, delta):
    def compliance(r_in, r_out, loaded_inner):
        if loaded_inner:
            hoop = (r_out**2 + r_in**2) / (r_out**2 - r_in**2)
            radial = -1.0
            r = r_in
        else:
            hoop = -(r_out**2 + r_in**2) / (r_out**2 - r_in**2)
            radial = -1.0
            r = r_out
        return r * (1.0 + poisson) / young * ((1.0 - poisson) * hoop - poisson * radial)

    return delta / (
        compliance(radius, outer, True) - compliance(inner, radius, False)
    )


def relative(got, expected):
    return abs(got - expected) / max(abs(expected), 1.0)


def main():
    problem = json.loads(PROBLEM_PATH.read_text())
    fine = predict(problem, FINE)
    medium = predict(problem, MEDIUM)
    lame = plane_stress_pressure(problem)
    station_error = [
        relative(lame, value) for value in fine["assembly_pressure_pa"]
    ]
    spread = {
        "pressure_rel_max": max(
            relative(a, b)
            for a, b in zip(medium["assembly_pressure_pa"], fine["assembly_pressure_pa"])
        ),
        "torque_rel": relative(
            medium["assembly_torque_capacity_nm"], fine["assembly_torque_capacity_nm"]
        ),
        "speed_rel": relative(medium["limit_speed_rad_s"], fine["limit_speed_rad_s"]),
    }
    too_far = (
        spread["pressure_rel_max"] > 0.01
        or spread["torque_rel"] > 0.01
        or spread["speed_rel"] > 0.01
    )
    if too_far:
        raise SystemExit(f"medium mesh is too far from the fine mesh: {spread}")
    if max(station_error) < 0.03:
        raise SystemExit("uniform Lamé pressure is inside the 3 percent station band")

    # Long, equal-material sleeve: the mid-face pressure should approach plane strain.
    check = json.loads(PROBLEM_PATH.read_text())
    check["sleeve_length_m"] = 0.080
    check["sleeve_z0_m"] = 0.020
    check["sleeve"] = json.loads(json.dumps(check["shaft"]))
    check_mesh = {
        "nr_shaft": 8,
        "nr_sleeve": 8,
        "nz_contact": 32,
        "nz_side": 4,
        "penalty_scale": 5000.0,
    }
    solved = contact_solution(
        check,
        0.0,
        check["assembly_temperature_c"],
        check["assembly_temperature_c"],
        check_mesh,
    )
    fem_mid = station_pressures(check, solved)[2]
    lame_long = plane_strain_pressure(
        check["shaft"]["young_modulus_pa"],
        check["shaft"]["poisson_ratio"],
        check["shaft_inner_radius_m"],
        check["shaft_outer_radius_m"],
        check["sleeve_outer_radius_m"],
        0.5 * check["diametral_interference_m"],
    )
    long_rel = relative(fem_mid, lame_long)
    if long_rel > 0.03:
        raise SystemExit(f"long-sleeve check is {long_rel:.3%} from plane strain")

    truth = {
        "assembly_pressure_pa": fine["assembly_pressure_pa"],
        "assembly_torque_capacity_nm": fine["assembly_torque_capacity_nm"],
        "limit_speed_rad_s": fine["limit_speed_rad_s"],
        "pressure_rel_tol": 0.03,
        "torque_rel_tol": 0.015,
        "speed_rel_tol": 0.02,
        "uniform_lame_pressure_pa": lame,
        "fine_mesh": FINE,
    }
    evidence = {
        "data": "synthetic",
        "generator": "authoring/generate_truth.py",
        "fine_mesh": FINE,
        "medium_mesh": MEDIUM,
        "medium_result": medium,
        "fine_result": fine,
        "mesh_spread": spread,
        "uniform_plane_stress_lame_pa": lame,
        "uniform_lame_station_relative_error": station_error,
        "long_equal_material_sleeve": {
            "fem_mid_pressure_pa": fem_mid,
            "plane_strain_lame_pa": lame_long,
            "relative_difference": long_rel,
        },
    }
    TRUTH_PATH.parent.mkdir(parents=True, exist_ok=True)
    TRUTH_PATH.write_text(json.dumps(truth, indent=2) + "\n")
    EVIDENCE_PATH.parent.mkdir(parents=True, exist_ok=True)
    EVIDENCE_PATH.write_text(json.dumps(evidence, indent=2) + "\n")
    print(
        json.dumps(
            {"spread": spread, "long_rel": long_rel, "limit": fine["limit_speed_rad_s"]},
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
