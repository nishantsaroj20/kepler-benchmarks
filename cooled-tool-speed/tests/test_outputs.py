"""Sealed checks for the cooled cutting-tool tip."""

import json
import math
from pathlib import Path

TRUTH = json.loads(Path("/tests/truth.json").read_text())
RESULT = Path("/app/output/tool_result.json")


def _relative_error(got, expected):
    return abs(got - expected) / max(abs(expected), 1.0)


def test_constant_conductivity_misses_the_rake_stations():
    """A temperature-independent conductivity has to miss the hot rake."""
    misses = [
        _relative_error(got, expected) > TRUTH["temperature_rel_tol"]
        for got, expected in zip(
            TRUTH["constant_conductivity_rake_temperature_c"],
            TRUTH["rake_temperature_c"],
            strict=True,
        )
    ]
    assert sum(misses) >= 4


def test_tool_result_matches_the_fitted_mesh():
    """Rake temperatures, edge growth, and limit speed match the fine mesh."""
    assert RESULT.is_file(), "tool_result.json is missing"
    result = json.loads(RESULT.read_text())
    temperatures = result["rake_temperature_c"]
    assert len(temperatures) == 5
    for got, expected in zip(temperatures, TRUTH["rake_temperature_c"], strict=True):
        assert math.isfinite(got)
        assert _relative_error(got, expected) <= TRUTH["temperature_rel_tol"]
    growth = result["edge_growth_um"]
    speed = result["limit_speed_m_s"]
    assert math.isfinite(growth) and math.isfinite(speed)
    assert _relative_error(growth, TRUTH["edge_growth_um"]) <= TRUTH["growth_rel_tol"]
    assert _relative_error(speed, TRUTH["limit_speed_m_s"]) <= TRUTH["speed_rel_tol"]
