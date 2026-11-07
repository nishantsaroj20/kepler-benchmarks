"""Sealed checks for the short-hub shrink fit."""

import json
import math
from pathlib import Path

TRUTH = json.loads(Path("/tests/truth.json").read_text())
RESULT = Path("/app/output/fit_result.json")


def _relative_error(got, expected):
    return abs(got - expected) / max(abs(expected), 1.0)


def test_uniform_long_cylinder_pressure_is_outside_the_station_band():
    """A single Lamé pressure has to miss the short-hub stations."""
    lame = TRUTH["uniform_lame_pressure_pa"]
    misses = [
        _relative_error(lame, pressure) > TRUTH["pressure_rel_tol"]
        for pressure in TRUTH["assembly_pressure_pa"]
    ]
    assert sum(misses) >= 3


def test_fit_result_matches_the_axisymmetric_solution():
    """Assembly pressures, assembly torque, and hot limit speed match the fine mesh."""
    assert RESULT.is_file(), "fit_result.json is missing"
    result = json.loads(RESULT.read_text())
    pressures = result["assembly_pressure_pa"]
    assert len(pressures) == 5
    for got, expected in zip(pressures, TRUTH["assembly_pressure_pa"]):
        assert math.isfinite(got)
        assert _relative_error(got, expected) <= TRUTH["pressure_rel_tol"]
    torque = result["assembly_torque_capacity_nm"]
    speed = result["limit_speed_rad_s"]
    assert math.isfinite(torque) and math.isfinite(speed)
    assert _relative_error(torque, TRUTH["assembly_torque_capacity_nm"]) <= TRUTH["torque_rel_tol"]
    assert _relative_error(speed, TRUTH["limit_speed_rad_s"]) <= TRUTH["speed_rel_tol"]
