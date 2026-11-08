"""Outcome-only checks. No participant code, numerical model, or input imports."""
import itertools
import json
import math
import os
from pathlib import Path
import pytest

ART = Path(os.environ.get("LDI_ARTIFACT", "/app/output/damage.json"))
GT = Path(__file__).with_name("ground_truth.json")
IDS = {"plate_01"}
BOUNDS = {
    "xc": (0.06, 0.44),
    "yc": (0.055, 0.345),
    "a": (0.014, 0.04),
    "b": (0.008, 0.02),
    "theta_deg": (-90.0, 90.0),
    "d": (0.4, 0.95),
}


def unique(pairs):
    found = {}
    for key, value in pairs:
        if key in found:
            raise ValueError("Duplicate JSON key: " + key)
        found[key] = value
    return found


def finite(x):
    assert type(x) in (int, float) and math.isfinite(x), (
        "Every physical value must be a finite JSON number, not a boolean or string"
    )
    return float(x)


def document():
    assert ART.is_file() and not ART.is_symlink(), "Missing regular output file (symlinks are not allowed)"
    assert 0 < ART.stat().st_size <= 1048576, "Output must be nonempty and no larger than 1 MiB"
    with ART.open() as handle:
        return json.load(
            handle,
            object_pairs_hook=unique,
            parse_constant=lambda token: (_ for _ in ()).throw(ValueError("Nonfinite JSON token " + token)),
        )


def entries():
    out = document()
    assert isinstance(out, dict) and isinstance(out.get("plates"), list), "Expected a root object with plates array"
    assert len(out["plates"]) == 1, "Exactly one plate is required"
    by = {}
    for entry in out["plates"]:
        assert isinstance(entry, dict), "Plate entries must be objects"
        pid = entry.get("plate_id")
        assert isinstance(pid, str) and pid in IDS and pid not in by, "Unknown or duplicated plate_id"
        by[pid] = entry
    assert set(by) == IDS
    return by


def test_regular_bounded_unique_json():
    document()


def test_schema_bounds_and_forecast_shape():
    for pid, entry in entries().items():
        assert isinstance(entry.get("zones"), list) and len(entry["zones"]) == 3, pid + " needs exactly three zones"
        for zone in entry["zones"]:
            assert isinstance(zone, dict)
            for key, (lo, hi) in BOUNDS.items():
                assert key in zone, pid + " missing " + key
                value = finite(zone[key])
                assert lo - 1e-12 <= value <= hi + 1e-12, pid + " outside inclusive bounds for " + key
            assert finite(zone["a"]) > finite(zone["b"]) - 1e-12, pid + " must keep the major semi-axis longer"
        assert isinstance(entry.get("fixture"), dict)
        for key in ("x0_multiplier", "yL_multiplier"):
            assert key in entry["fixture"]
            value = finite(entry["fixture"][key])
            assert 0.6 - 1e-12 <= value <= 1.6 + 1e-12
        assert isinstance(entry.get("branches"), dict) and set(entry["branches"]) == {f"record_{i:02d}" for i in range(1, 5)}
        assert all(finite(value) in (-1.0, 1.0) for value in entry["branches"].values()), (
            "Branch signs must be numeric -1 or +1, not booleans or strings"
        )
        assert "snap_load_N" in entry and finite(entry["snap_load_N"]) > 0, "Snap load must be finite and positive"
        forecast = entry.get("open_branch_forecast_Hz")
        assert isinstance(forecast, list) and len(forecast) == 2, pid + " needs two open-branch forecast rows"
        for row in forecast:
            assert isinstance(row, dict) and "force_N" in row and "frequencies_Hz" in row
            finite(row["force_N"])
            freqs = row["frequencies_Hz"]
            assert isinstance(freqs, list) and len(freqs) == 4, pid + " needs four forecast frequencies"
            assert all(finite(value) > 0 for value in freqs)


def angle_gap(left, right):
    return abs((left - right + 90.0) % 180.0 - 90.0)


def zone_matches(got, target, tol):
    return (
        abs(finite(got["xc"]) - target["xc"]) <= tol["zone_centre_m"] + 1e-12
        and abs(finite(got["yc"]) - target["yc"]) <= tol["zone_centre_m"] + 1e-12
        and abs(finite(got["a"]) - target["a"]) <= tol["zone_a_m"] + 1e-12
        and abs(finite(got["b"]) - target["b"]) <= tol["zone_b_m"] + 1e-12
        and angle_gap(finite(got["theta_deg"]), target["theta_deg"]) <= tol["zone_theta_deg"] + 1e-12
        and abs(finite(got["d"]) - target["d"]) <= tol["zone_d"] + 1e-12
    )


def test_preparations_snap_and_open_branch_forecast():
    by = entries()
    truth = json.loads(GT.read_text())
    tol = truth["tolerances"]
    for plate in truth["plates"]:
        entry = by[plate["plate_id"]]
        reported = entry["zones"]
        authored = plate["zones"]
        paired = False
        for order in itertools.permutations(range(len(authored))):
            if all(zone_matches(reported[i], authored[order[i]], tol) for i in range(len(authored))):
                paired = True
                break
        assert paired, plate["plate_id"] + " zones do not match the physical patches"
        for key in ("x0_multiplier", "yL_multiplier"):
            assert abs(finite(entry["fixture"][key]) - plate["fixture"][key]) <= tol["fixture_absolute"] + 1e-12, (
                plate["plate_id"] + " spring multiplier " + key + " is outside its limit"
            )
        assert all(
            abs(finite(entry["branches"][key]) - plate["branches"][key]) <= 1e-12
            for key in plate["branches"]
        ), plate["plate_id"] + " has incorrect prepared equilibrium signs"
        assert abs(finite(entry["snap_load_N"]) / plate["snap_load_N"] - 1) <= tol["snap_relative"] + 1e-12, (
            plate["plate_id"] + " has an inaccurate closed-branch snap load"
        )
        reported = entry["open_branch_forecast_Hz"]
        for got, target in zip(reported, plate["open_branch_forecast_Hz"]):
            assert abs(finite(got["force_N"]) - target["force_N"]) <= 1e-12, "Forecast forces must be the requested values"
            for left, right in zip(got["frequencies_Hz"], target["frequencies_Hz"]):
                assert abs(finite(left) / right - 1) <= tol["forecast_relative"] + 1e-12, (
                    "Open-branch frequency forecast is outside its relative limit"
                )
