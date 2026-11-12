#!/usr/bin/env python3
"""Verifier: re-run analyzer on shifted held-out packs; exercise practice + schema."""

from __future__ import annotations

import json
import math
import os
from pathlib import Path

import pytest

from harness import PROGRAM, PER_RUN_TIMEOUT_SEC, load_output_json, run_analyzer
from hidden_family import materialize_heldout
from truth_solver import truth_for_pack

N_PACKS = 8
PRACTICE_ROOT = Path("/app/data/practice")

REQUIRED_KEYS = [
    "schema_version",
    "p_peak",
    "fn_total",
    "ft_total",
    "contact_length",
    "q_peak",
    "stick_fraction",
    "slip_length",
    "x_q_centroid",
    "cum_slip_l1",
]

REL_KEYS = {
    "p_peak": 0.05,
    "fn_total": 0.02,
    "ft_total": 0.02,
    "contact_length": 0.10,
    "q_peak": 0.10,
    "slip_length": 0.15,
    "cum_slip_l1": 0.30,
}
ABS_KEYS = {
    "stick_fraction": 0.08,
    "x_q_centroid": 1.0e-4,
}


def _is_json_int(x) -> bool:
    return isinstance(x, int) and not isinstance(x, bool)


def _finite_number(x) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(float(x))


def _validate_schema(data: dict) -> None:
    assert isinstance(data, dict)
    assert set(data.keys()) == set(REQUIRED_KEYS), f"key set mismatch: {sorted(data.keys())}"
    assert _is_json_int(data["schema_version"]) and data["schema_version"] == 1
    for key in REQUIRED_KEYS:
        if key == "schema_version":
            continue
        assert _finite_number(data[key]), f"non-finite or bad type for {key}: {data[key]!r}"


def _compare(pred: dict, truth: dict, name: str) -> None:
    for key, tol in REL_KEYS.items():
        p, t = float(pred[key]), float(truth[key])
        if abs(t) < 1e-30 and abs(p) < 1e-30:
            continue
        # relative; for near-zero truth with nonzero pred, fail
        den = max(abs(t), 1e-30)
        err = abs(p - t) / den
        assert err <= tol, f"{name} {key}: rel err {err} > {tol} (pred={p}, truth={t})"
    for key, tol in ABS_KEYS.items():
        p, t = float(pred[key]), float(truth[key])
        err = abs(p - t)
        assert err <= tol, f"{name} {key}: abs err {err} > {tol} (pred={p}, truth={t})"


def _run_and_score(pack: Path, out_json: Path) -> dict:
    result = run_analyzer(pack, out_json)
    assert not result["timeout"], f"timeout on {pack.name}"
    assert result["status"] == 0, f"nonzero exit on {pack.name}: {result['status']}"
    assert result["seconds"] <= PER_RUN_TIMEOUT_SEC
    pred = load_output_json(out_json)
    _validate_schema(pred)
    return pred


def test_program_present():
    assert PROGRAM.is_file(), f"required artifact missing: {PROGRAM}"
    assert not PROGRAM.is_symlink(), "program must not be a symlink"
    text = PROGRAM.read_text(encoding="utf-8", errors="ignore")
    assert len(text) > 50


def test_practice_packs_schema_and_science():
    """Instruction requires the program to work on every public practice pack."""
    assert PRACTICE_ROOT.is_dir(), "practice data missing from verifier image"
    packs = sorted(p for p in PRACTICE_ROOT.iterdir() if p.is_dir())
    assert len(packs) >= 1
    import tempfile

    with tempfile.TemporaryDirectory(prefix="practice_run_") as td:
        root = Path(td)
        for pack in packs:
            out = root / f"{pack.name}.json"
            pred = _run_and_score(pack, out)
            truth = truth_for_pack(pack)
            _compare(pred, truth, pack.name)


def test_heldout_shifted_packs():
    assert PROGRAM.is_file()
    import tempfile

    with tempfile.TemporaryDirectory(prefix="held_packs_") as td:
        root = Path(td)
        packs = materialize_heldout(PROGRAM, root / "packs", n_packs=N_PACKS)
        assert len(packs) == N_PACKS
        # ensure at least one no-contact and frictional case among drawn family
        saw_contact = False
        saw_open = False
        for pack in packs:
            out = root / f"{pack.name}.json"
            pred = _run_and_score(pack, out)
            truth = truth_for_pack(pack)
            _compare(pred, truth, pack.name)
            if truth["contact_length"] > 0:
                saw_contact = True
            else:
                saw_open = True
                # no-contact contract
                assert pred["p_peak"] == 0.0
                assert pred["fn_total"] == 0.0
                assert pred["ft_total"] == 0.0
                assert pred["q_peak"] == 0.0
                assert pred["cum_slip_l1"] == 0.0
        assert saw_contact and saw_open


def test_schema_rejects_extras_and_bool_schema_version():
    bad = {
        "schema_version": True,  # bool must not coerce
        "p_peak": 1.0,
        "fn_total": 1.0,
        "ft_total": 1.0,
        "contact_length": 1.0,
        "q_peak": 1.0,
        "stick_fraction": 0.5,
        "slip_length": 0.1,
        "x_q_centroid": 0.0,
        "cum_slip_l1": 0.0,
    }
    with pytest.raises(AssertionError):
        _validate_schema(bad)
    bad2 = dict(bad)
    bad2["schema_version"] = 1
    bad2["extra"] = 0
    with pytest.raises(AssertionError):
        _validate_schema(bad2)


def test_baselines_fail():
    import tempfile

    with tempfile.TemporaryDirectory(prefix="baseline_") as td:
        root = Path(td)
        packs = materialize_heldout(PROGRAM, root / "packs", n_packs=1)
        # pick a contacting pack if possible by drawing until contact
        truth = truth_for_pack(packs[0])
        if truth["contact_length"] == 0:
            packs = materialize_heldout(PROGRAM, root / "packs2", n_packs=8)
            packs = [p for p in packs if truth_for_pack(p)["contact_length"] > 0]
            assert packs
            truth = truth_for_pack(packs[0])
        zeros = {k: (1 if k == "schema_version" else 0.0) for k in REQUIRED_KEYS}
        with pytest.raises(AssertionError):
            _compare(zeros, truth, "zeros")
        # copy final loads style: cannot invent frictional field from δ alone without solve
        bogus = dict(zeros)
        bogus.update(
            {
                "p_peak": 1e6,
                "fn_total": 1.0,
                "ft_total": 1.0,
                "contact_length": 1.0,
                "q_peak": 1.0,
            }
        )
        with pytest.raises(AssertionError):
            _compare(bogus, truth, "bogus")
