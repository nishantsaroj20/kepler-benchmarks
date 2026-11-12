#!/usr/bin/env python3
"""Verifier: run the submitted column_solver.py on held-out hidden cases."""

from __future__ import annotations

import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import solve_ivp

PROGRAM = Path("/workspace/output/column_solver.py")
HIDDEN_ROOT = Path("/tests/hidden")
MANIFEST = HIDDEN_ROOT / "manifest.json"
RESULT_DIR = Path("/tmp/adsorption_hidden_results")

_REPO = Path(__file__).resolve().parents[1]
if not PROGRAM.exists():
    local_prog = _REPO / "output" / "column_solver.py"
    if local_prog.exists():
        PROGRAM = local_prog
if not MANIFEST.exists():
    local_hidden = _REPO / "tests" / "hidden"
    if (local_hidden / "manifest.json").exists():
        HIDDEN_ROOT = local_hidden
        MANIFEST = local_hidden / "manifest.json"

REQUIRED_KEYS = ("qm", "K", "k_LDF", "D_ax", "L_min_m")


def _thread_env():
    env = os.environ.copy()
    for key in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "BLIS_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
    ):
        env[key] = "1"
    return env


def q_star(c, qm, K):
    return qm * K * c / (1.0 + K * c)


def max_outlet_under_truth(case: dict, params: dict, L: float, n_nodes: int = 64) -> float:
    """Independent MOL check used for L_min feasibility (not the agent's code)."""
    duty = case["design_duty"]
    eps = float(case["bed_voidage"])
    rho = float(case["particle_density_kg_m3"])
    c_in = float(duty["c_in_mol_m3"])
    u = float(duty["interstitial_velocity_m_s"])
    t_service = float(duty["service_time_s"])
    qm = float(params["qm"])
    K = float(params["K"])
    k_ldf = float(params["k_LDF"])
    D_ax = float(params["D_ax"])
    dz = L / n_nodes
    factor = ((1.0 - eps) / eps) * rho
    y0 = np.zeros(2 * n_nodes)

    def rhs(_t, y):
        c = y[:n_nodes]
        q = y[n_nodes:]
        dqdt = k_ldf * (q_star(c, qm, K) - q)
        dc = np.empty(n_nodes)
        dc[0] = u * (c[0] - c_in) / dz
        dc[1:] = u * (c[1:] - c[:-1]) / dz
        c_im1 = np.empty(n_nodes)
        c_ip1 = np.empty(n_nodes)
        c_im1[0] = c_in
        c_im1[1:] = c[:-1]
        c_ip1[:-1] = c[1:]
        c_ip1[-1] = c[-1]
        d2 = (c_ip1 - 2.0 * c + c_im1) / (dz * dz)
        return np.concatenate([-dc + D_ax * d2 - factor * dqdt, dqdt])

    sol = solve_ivp(
        rhs,
        (0.0, t_service),
        y0,
        method="BDF",
        t_eval=np.linspace(0.0, t_service, 201),
        rtol=1e-6,
        atol=1e-9,
    )
    assert sol.success, sol.message
    return float(np.max(sol.y[n_nodes - 1, :]))


@pytest.fixture(scope="module")
def manifest():
    assert MANIFEST.is_file(), f"missing manifest {MANIFEST}"
    data = json.loads(MANIFEST.read_text())
    assert isinstance(data, dict)
    return data


@pytest.fixture(scope="module")
def program_path():
    assert PROGRAM.is_file(), f"missing program artifact {PROGRAM}"
    assert not PROGRAM.is_symlink(), "program must not be a symlink"
    text = PROGRAM.read_text(errors="ignore")
    assert text.strip(), "program file is empty"
    return PROGRAM


@pytest.fixture(scope="module")
def hidden_results(program_path, manifest):
    timeout_s = float(manifest["per_run_timeout_s"])
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    results = {}
    env = _thread_env()
    for entry in manifest["hidden_cases"]:
        case_id = entry["case_id"]
        case_path = HIDDEN_ROOT / case_id / "case.json"
        out_path = RESULT_DIR / f"{case_id}.json"
        if out_path.exists():
            out_path.unlink()
        assert case_path.is_file(), f"missing hidden case {case_path}"
        cmd = [
            sys.executable,
            str(program_path),
            "--case",
            str(case_path),
            "--out",
            str(out_path),
        ]
        t0 = time.monotonic()
        try:
            proc = subprocess.run(
                cmd,
                env=env,
                capture_output=True,
                text=True,
                timeout=timeout_s,
            )
        except subprocess.TimeoutExpired as exc:
            pytest.fail(
                f"case {case_id}: program exceeded {timeout_s}s per-run limit "
                f"(timeout fails the task); stdout={exc.stdout!r} stderr={exc.stderr!r}"
            )
        elapsed = time.monotonic() - t0
        if proc.returncode != 0:
            pytest.fail(
                f"case {case_id}: program exited {proc.returncode} in {elapsed:.2f}s\n"
                f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
            )
        assert out_path.is_file(), f"case {case_id}: no output written"
        assert not out_path.is_symlink(), f"case {case_id}: symlink output rejected"
        payload = json.loads(out_path.read_text())
        assert isinstance(payload, dict), f"case {case_id}: output must be a JSON object"
        case = json.loads(case_path.read_text())
        results[case_id] = {
            "payload": payload,
            "elapsed_s": elapsed,
            "truth": entry,
            "case": case,
        }
    return results


def _finite_positive(value, name):
    assert isinstance(value, (int, float)) and not isinstance(value, bool), (
        f"{name} must be numeric"
    )
    value = float(value)
    assert math.isfinite(value), f"{name} must be finite"
    assert value > 0.0, f"{name} must be positive"
    return value


def test_program_exists(program_path):
    assert program_path.name.endswith(".py")


def test_all_hidden_cases_run(hidden_results, manifest):
    expected = {e["case_id"] for e in manifest["hidden_cases"]}
    assert set(hidden_results) == expected


@pytest.mark.parametrize(
    "case_id",
    [
        "H1_new_material",
        "H2_dispersion_regime",
        "H3_sensor_faults",
        "H4_high_cin",
        "H5_near_linear",
    ],
)
def test_hidden_case_accuracy(case_id, hidden_results, manifest):
    assert case_id in hidden_results, f"missing result for {case_id}"
    blob = hidden_results[case_id]
    payload = blob["payload"]
    truth = blob["truth"]
    case = blob["case"]
    tol = truth["tolerances"]
    for key in REQUIRED_KEYS:
        assert key in payload, f"{case_id}: missing key {key}"
        _finite_positive(payload[key], f"{case_id}.{key}")

    params = truth["parameters"]
    for key, rel in (
        ("qm", tol["qm_rel"]),
        ("K", tol["K_rel"]),
        ("k_LDF", tol["k_LDF_rel"]),
        ("D_ax", tol["D_ax_rel"]),
    ):
        est = float(payload[key])
        ref = float(params[key])
        assert abs(est - ref) / abs(ref) <= rel + 1e-12, (
            f"{case_id}: {key} rel err too large est={est} ref={ref} tol={rel}"
        )

    L_est = float(payload["L_min_m"])
    L_ref = float(truth["L_min_m"])
    # Anti-trivial: reject beds that ignore sizing and just pick a huge length.
    assert L_est <= float(tol["L_min_oversize_ceil"]) * L_ref + 1e-12, (
        f"{case_id}: L_min trivially oversized est={L_est} ref={L_ref}"
    )
    # Physical feasibility under the true model (standard-method friendly).
    c_in = float(case["design_duty"]["c_in_mol_m3"])
    c_limit = float(case["design_duty"]["outlet_mole_fraction_limit"]) * c_in
    slack = float(tol["outlet_feasibility_slack"])
    max_out = max_outlet_under_truth(case, params, L_est, n_nodes=64)
    assert max_out <= c_limit * slack + 1e-12, (
        f"{case_id}: true-model outlet {max_out} exceeds limit {c_limit} at L={L_est}"
    )
