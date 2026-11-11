"""Check the slurry-rate predictions against the latent campaign truth."""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

OUT = Path("/app/output/predictions.json")
TRUTH = json.loads(Path("/tests/truth/expected.json").read_text())
REL_TOL = 0.040
ABS_TOL = 0.0006
FRAC_TOL = 0.030
SUM_TOL = 0.0001
MAX_BYTES = 65536
NUMERIC_KEYS = (
    "rate_mol_m3_s",
    "resistance_fraction_gas_liquid",
    "resistance_fraction_liquid_solid",
    "resistance_fraction_kinetic",
)
REQUIRED_KEYS = ("case_id",) + NUMERIC_KEYS
CASE_IDS = [row["case_id"] for row in TRUTH["cases"]]
TRUTH_BY_ID = {row["case_id"]: row for row in TRUTH["cases"]}


def _reject_constant(token: str):
    raise ValueError(token)


def _pairs(pairs):
    obj = {}
    for key, value in pairs:
        if key in obj:
            raise ValueError(f"duplicate key {key}")
        obj[key] = value
    return obj


def load_document():
    assert not OUT.is_symlink()
    assert OUT.is_file()
    raw = OUT.read_bytes()
    assert 0 < len(raw) <= MAX_BYTES
    return json.loads(raw.decode("utf-8"), object_pairs_hook=_pairs, parse_constant=_reject_constant)


def case_map(document):
    assert isinstance(document, dict)
    assert set(document.keys()) == {"cases"}
    rows = document["cases"]
    assert isinstance(rows, list)
    assert rows, "cases is empty"
    seen = []
    mapped = {}
    for row in rows:
        assert isinstance(row, dict)
        assert set(row.keys()) == set(REQUIRED_KEYS)
        assert len(row.keys()) == len(REQUIRED_KEYS)
        case_id = row["case_id"]
        assert isinstance(case_id, str)
        assert case_id not in mapped
        seen.append(case_id)
        mapped[case_id] = row
    assert seen == CASE_IDS or set(seen) == set(CASE_IDS)
    assert set(seen) == set(CASE_IDS)
    assert len(seen) == len(CASE_IDS)
    return mapped


def finite_number(value):
    assert not isinstance(value, bool)
    assert isinstance(value, (int, float))
    number = float(value)
    assert math.isfinite(number)
    return number


def test_predictions_file_is_regular_and_bounded():
    assert not OUT.is_symlink()
    assert OUT.is_file()
    size = OUT.stat().st_size
    assert 0 < size <= MAX_BYTES


def test_document_parses_without_nonfinite_tokens():
    document = load_document()
    assert isinstance(document, dict)


def test_case_ids_match_exactly():
    mapped = case_map(load_document())
    assert set(mapped) == set(CASE_IDS)
    assert len(mapped) == len(CASE_IDS)


@pytest.mark.parametrize("case_id", CASE_IDS)
def test_rate(case_id: str):
    row = case_map(load_document())[case_id]
    predicted = finite_number(row["rate_mol_m3_s"])
    truth = float(TRUTH_BY_ID[case_id]["rate_mol_m3_s"])
    tolerance = max(REL_TOL * abs(truth), ABS_TOL)
    assert abs(predicted - truth) <= tolerance


@pytest.mark.parametrize("case_id", CASE_IDS)
def test_fractions(case_id: str):
    row = case_map(load_document())[case_id]
    values = []
    for key in NUMERIC_KEYS[1:]:
        predicted = finite_number(row[key])
        truth = float(TRUTH_BY_ID[case_id][key])
        assert abs(predicted - truth) <= FRAC_TOL
        values.append(predicted)
    assert abs(sum(values) - 1.0) <= SUM_TOL
