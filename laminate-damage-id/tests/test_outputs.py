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
        assert isinstance(entry.get("zones"), list) and len(entry["zones"]) == 4, pid + " needs exactly four zones"
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
            values = [finite(value) for value in freqs]
            assert all(value > 0 for value in values)
            assert values == sorted(values), "Forecast frequencies must be sorted from low to high"


def angle_gap(left, right):
    return abs((left - right + 90.0) % 180.0 - 90.0)


ZONE_LIMITS = {
    "xc": "zone_centre_m",
    "yc": "zone_centre_m",
    "a": "zone_a_m",
    "b": "zone_b_m",
    "theta_deg": "zone_theta_deg",
    "d": "zone_d",
}


def zone_error(got, target, key):
    if key == "theta_deg":
        return angle_gap(finite(got[key]), target[key])
    return abs(finite(got[key]) - target[key])


def scored_plates():
    by_id = entries()
    truth = json.loads(GT.read_text())
    for target in truth["plates"]:
        yield by_id[target["plate_id"]], target, truth["tolerances"]


def test_zone_pairing():
    for entry, plate, tol in scored_plates():
        reported = entry["zones"]
        authored = plate["zones"]
        limits = {key: tol[name] + 1e-12 for key, name in ZONE_LIMITS.items()}
        # Minimize the worst normalized error over complete one-to-one pairings.
        # The same all-attribute limits apply to every zone in the chosen pairing.
        worst = None
        for order in itertools.permutations(range(len(authored))):
            candidate = max(
                (zone_error(zone, authored[j], key) / limit, i, key)
                for i, (zone, j) in enumerate(zip(reported, order))
                for key, limit in limits.items()
            )
            if worst is None or candidate < worst:
                worst = candidate
        ratio, index, key = worst
        assert ratio <= 1, (
            f"{plate['plate_id']} has no valid one-to-one zone pairing. "
            f"In the closest pairing, reported zone {index + 1} exceeds the {key} limit "
            f"({ratio:.3f} times the allowed error)."
        )


def test_fixture_and_preparations():
    for entry, plate, tol in scored_plates():
        for key in ("x0_multiplier", "yL_multiplier"):
            error = abs(finite(entry["fixture"][key]) - plate["fixture"][key])
            assert error <= tol["fixture_absolute"] + 1e-12, (
                f"{plate['plate_id']} {key} error {error:.6g} exceeds {tol['fixture_absolute']}"
            )
        assert all(
            abs(finite(entry["branches"][key]) - plate["branches"][key]) <= 1e-12
            for key in plate["branches"]
        ), plate["plate_id"] + " has incorrect prepared equilibrium signs"


def test_snap_load():
    for entry, plate, tol in scored_plates():
        error = abs(finite(entry["snap_load_N"]) / plate["snap_load_N"] - 1)
        assert error <= tol["snap_relative"] + 1e-12, (
            f"{plate['plate_id']} closed-branch snap relative error {error:.6g} "
            f"exceeds {tol['snap_relative']}"
        )


def test_open_branch_forecast():
    for entry, plate, tol in scored_plates():
        reported = entry["open_branch_forecast_Hz"]
        for got, target in zip(reported, plate["open_branch_forecast_Hz"]):
            assert abs(finite(got["force_N"]) - target["force_N"]) <= 1e-12, "Forecast forces must be the requested values"
            for mode, (left, right) in enumerate(zip(got["frequencies_Hz"], target["frequencies_Hz"]), 1):
                error = abs(finite(left) / right - 1)
                assert error <= tol["forecast_relative"] + 1e-12, (
                    f"Open-branch forecast at {got['force_N']} N, mode {mode}: "
                    f"relative error {error:.6g} exceeds {tol['forecast_relative']}"
                )
