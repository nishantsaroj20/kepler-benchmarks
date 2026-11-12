"""Compare the forecast file with the latent extraction truth."""

from __future__ import annotations

import json
import math
from pathlib import Path

PREDICTIONS = Path("/app/output/predictions.json")
TRUTH_PATH = Path("/tests/truth/expected.json")
REL_TOL = 0.02
ABS_RECOVERY = 0.012
ABS_ORGANIC = 0.60
MAX_BYTES = 2_000_000


def _pairs(pairs):
    obj = {}
    for key, value in pairs:
        if key in obj:
            raise ValueError(f"duplicate key {key}")
        obj[key] = value
    return obj


def _reject_constant(token: str):
    raise ValueError(token)


def load_document(path: Path):
    assert path.is_file()
    assert not path.is_symlink()
    resolved = path.resolve()
    assert resolved.parent == path.parent.resolve()
    raw = path.read_bytes()
    assert 0 < len(raw) <= MAX_BYTES
    return json.loads(raw.decode("utf-8"), object_pairs_hook=_pairs, parse_constant=_reject_constant)


def finite_number(value):
    assert not isinstance(value, bool)
    assert isinstance(value, (int, float))
    number = float(value)
    assert math.isfinite(number)
    return number


def test_forecast_matches_truth():
    document = load_document(PREDICTIONS)
    truth = json.loads(TRUTH_PATH.read_text(encoding="utf-8"))
    assert set(document.keys()) == {"charges"}
    rows = document["charges"]
    truth_rows = truth["charges"]
    assert isinstance(rows, list) and len(rows) == len(truth_rows)
    missed = []
    for got, truth_row in zip(rows, truth_rows):
        assert isinstance(got, dict)
        assert set(got.keys()) == {"case_id", "recovery", "organic_mol_m3"}
        assert got["case_id"] == truth_row["case_id"]
        recovery = finite_number(got["recovery"])
        organic = finite_number(got["organic_mol_m3"])
        true_recovery = float(truth_row["recovery"])
        true_organic = float(truth_row["organic_mol_m3"])
        if abs(recovery - true_recovery) > max(REL_TOL * abs(true_recovery), ABS_RECOVERY):
            missed.append(truth_row["case_id"])
            continue
        if abs(organic - true_organic) > max(REL_TOL * abs(true_organic), ABS_ORGANIC):
            missed.append(truth_row["case_id"])
    assert not missed, f"{len(missed)} charges outside the stated limits, including {missed[:8]}"
