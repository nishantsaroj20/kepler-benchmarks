"""Sealed checks for the coupled on-chip bus."""

import json
import math
from pathlib import Path

TRUTH = json.loads(Path("/tests/truth.json").read_text())
RESULT = Path("/app/output/bus_timing.json")


def _relative_error(got, expected):
    return abs(got - expected) / max(abs(expected), 1.0)


def test_uncoupled_wires_miss_the_crossing_band():
    """Dropping the neighbor capacitance has to miss the far-end times."""
    misses = [
        _relative_error(got, expected) > TRUTH["cross_rel_tol"]
        for got, expected in zip(
            TRUTH["uncoupled_far_end_cross_ps"],
            TRUTH["far_end_cross_ps"],
            strict=True,
        )
    ]
    assert sum(misses) >= 4


def test_bus_timing_matches_the_coupled_transient():
    """Far-end crossings and the limit clock match the fine transient."""
    assert RESULT.is_file(), "bus_timing.json is missing"
    result = json.loads(RESULT.read_text())
    crossings = result["far_end_cross_ps"]
    assert len(crossings) == 5
    for got, expected in zip(crossings, TRUTH["far_end_cross_ps"], strict=True):
        assert math.isfinite(got)
        assert _relative_error(got, expected) <= TRUTH["cross_rel_tol"]
    clock = result["limit_clock_mhz"]
    assert math.isfinite(clock)
    assert _relative_error(clock, TRUTH["limit_clock_mhz"]) <= TRUTH["clock_rel_tol"]
