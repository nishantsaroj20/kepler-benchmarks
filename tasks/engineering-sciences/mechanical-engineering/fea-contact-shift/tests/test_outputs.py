#!/usr/bin/env python3
"""Verifier: re-run the submitted analyzer on shifted held-out packs."""

from __future__ import annotations

import json
import math
import os
import resource
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import pytest

from hidden_family import materialize_heldout
from truth_solver import truth_for_pack

PROGRAM = Path("/app/solver/analyze.py")
N_PACKS = 8
# Per-run wall budget ~3x observed oracle (~0.1 s) with generous CI margin.
PER_RUN_TIMEOUT_SEC = 8.0
OMP_ENV = {
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
}

REL_KEYS = {
    "p_peak": 0.12,
    "p_rms": 0.12,
    "fn_total": 0.02,
    "ft_total": 0.08,
    "contact_length": 0.12,
}
ABS_KEYS = {
    "mobilization": 0.05,
    "x_p_centroid": 5.0e-4,
}


def _finite(x) -> bool:
    try:
        return math.isfinite(float(x))
    except (TypeError, ValueError):
        return False


def _run_program(pack: Path, out_json: Path) -> None:
    if not PROGRAM.is_file():
        pytest.fail(f"missing program artifact: {PROGRAM}")

    def _limit():
        # soft CPU/AS limits; wall enforced by subprocess timeout
        resource.setrlimit(resource.RLIMIT_CPU, (int(PER_RUN_TIMEOUT_SEC) + 2, int(PER_RUN_TIMEOUT_SEC) + 2))
        mem = 2048 * 1024 * 1024
        resource.setrlimit(resource.RLIMIT_AS, (mem, mem))
        resource.setrlimit(resource.RLIMIT_NOFILE, (256, 256))
        resource.setrlimit(resource.RLIMIT_NPROC, (64, 64))

    env = os.environ.copy()
    env.update(OMP_ENV)
    cmd = [
        sys.executable,
        str(PROGRAM),
        "--pack",
        str(pack),
        "--out",
        str(out_json),
    ]
    try:
        subprocess.run(
            cmd,
            check=True,
            timeout=PER_RUN_TIMEOUT_SEC,
            env=env,
            preexec_fn=_limit,
            cwd="/tmp",
            stdin=subprocess.DEVNULL,
            capture_output=True,
            text=True,
        )
    except subprocess.TimeoutExpired as exc:
        pytest.fail(f"analyzer exceeded {PER_RUN_TIMEOUT_SEC}s on {pack.name}: {exc}")
    except subprocess.CalledProcessError as exc:
        pytest.fail(
            f"analyzer failed on {pack.name} with code {exc.returncode}; stderr={exc.stderr[-500:]}"
        )


def _load_out(path: Path) -> dict:
    assert path.is_file(), f"missing output JSON: {path}"
    data = json.loads(path.read_text())
    assert isinstance(data, dict), "output must be a JSON object"
    required = [
        "schema_version",
        "p_peak",
        "p_rms",
        "fn_total",
        "ft_total",
        "contact_length",
        "x_p_centroid",
        "mobilization",
    ]
    for key in required:
        assert key in data, f"missing key {key}"
    assert int(data["schema_version"]) == 1
    for key in required:
        if key == "schema_version":
            continue
        assert _finite(data[key]), f"non-finite {key}"
    return data


def _compare(pred: dict, truth: dict, pack: Path) -> None:
    for key, tol in REL_KEYS.items():
        p, t = float(pred[key]), float(truth[key])
        den = max(abs(t), 1e-30)
        err = abs(p - t) / den
        assert err <= tol, f"{pack.name} {key}: rel err {err} > {tol} (pred={p}, truth={t})"
    for key, tol in ABS_KEYS.items():
        p, t = float(pred[key]), float(truth[key])
        err = abs(p - t)
        assert err <= tol, f"{pack.name} {key}: abs err {err} > {tol} (pred={p}, truth={t})"


def test_program_present():
    assert PROGRAM.is_file(), f"required artifact missing: {PROGRAM}"
    text = PROGRAM.read_text(encoding="utf-8", errors="ignore")
    assert len(text) > 50, "program file looks empty"


def test_heldout_shifted_packs():
    assert PROGRAM.is_file(), f"required artifact missing: {PROGRAM}"
    with tempfile.TemporaryDirectory(prefix="held_packs_") as td:
        root = Path(td)
        packs = materialize_heldout(PROGRAM, root / "packs", n_packs=N_PACKS)
        assert len(packs) == N_PACKS
        for pack in packs:
            truth = truth_for_pack(pack)
            out_json = root / f"{pack.name}.json"
            t0 = time.time()
            _run_program(pack, out_json)
            elapsed = time.time() - t0
            assert elapsed <= PER_RUN_TIMEOUT_SEC
            pred = _load_out(out_json)
            _compare(pred, truth, pack)


def test_baselines_would_fail_on_contract():
    """Sanity: trivial constant / zero outputs violate finite-contact expectations on a live pack."""
    with tempfile.TemporaryDirectory(prefix="baseline_") as td:
        root = Path(td)
        packs = materialize_heldout(PROGRAM if PROGRAM.is_file() else Path(__file__), root / "packs", n_packs=1)
        truth = truth_for_pack(packs[0])
        zeros = {
            "schema_version": 1,
            "p_peak": 0.0,
            "p_rms": 0.0,
            "fn_total": 0.0,
            "ft_total": 0.0,
            "contact_length": 0.0,
            "x_p_centroid": 0.0,
            "mobilization": 0.0,
        }
        with pytest.raises(AssertionError):
            _compare(zeros, truth, packs[0])
        meanish = dict(zeros)
        meanish.update(
            {
                "p_peak": truth["p_peak"] * 0.0 + 1e6,
                "fn_total": 1.0,
                "contact_length": 1.0,
            }
        )
        with pytest.raises(AssertionError):
            _compare(meanish, truth, packs[0])
