#!/usr/bin/env python3
"""Harbor verifier for adsorption-breakthrough-id."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import solve_ivp

ARTIFACT = Path("/workspace/output/adsorption_identification.json")
GT_PATH = Path("/tests/ground_truth_params.json")
SPECS_PATH = Path("/workspace/data/bed_specs.json")

# Local fallbacks for host-side validation without Harbor mounts.
_REPO = Path(__file__).resolve().parents[1]
if not ARTIFACT.exists():
    ARTIFACT = _REPO / "output" / "adsorption_identification.json"
if not GT_PATH.exists():
    GT_PATH = Path(__file__).resolve().parent / "ground_truth_params.json"
if not SPECS_PATH.exists():
    SPECS_PATH = _REPO / "environment" / "data" / "bed_specs.json"

REQUIRED_KEYS = ("qm", "K", "k_LDF", "D_ax", "L_min_m")


def _load_json(path: Path) -> dict:
    assert path.is_file(), f"Missing artifact or resource: {path}"
    assert not path.is_symlink(), f"Symlink rejected: {path}"
    data = json.loads(path.read_text())
    assert isinstance(data, dict), f"Expected JSON object in {path}"
    return data


def _finite_number(value, name: str) -> float:
    assert isinstance(value, (int, float)) and not isinstance(value, bool), (
        f"{name} must be a number"
    )
    value = float(value)
    assert math.isfinite(value), f"{name} must be finite"
    return value


def q_star(c, qm, K):
    return qm * K * c / (1.0 + K * c)


def simulate_outlet_max(
    *,
    L,
    u,
    c_in,
    t_service,
    qm,
    K,
    k_ldf,
    D_ax,
    eps,
    rho_p,
    n_nodes=56,
):
    dz = L / n_nodes
    factor = ((1.0 - eps) / eps) * rho_p
    y0 = np.zeros(2 * n_nodes)

    def rhs(_t, y):
        c = y[:n_nodes]
        q = y[n_nodes:]
        dqdt = k_ldf * (q_star(c, qm, K) - q)
        dc_conv = np.empty(n_nodes)
        dc_conv[0] = u * (c[0] - c_in) / dz
        dc_conv[1:] = u * (c[1:] - c[:-1]) / dz
        c_im1 = np.empty(n_nodes)
        c_ip1 = np.empty(n_nodes)
        c_im1[0] = c_in
        c_im1[1:] = c[:-1]
        c_ip1[:-1] = c[1:]
        c_ip1[-1] = c[-1]
        d2c = (c_ip1 - 2.0 * c + c_im1) / (dz * dz)
        dcdt = -dc_conv + D_ax * d2c - factor * dqdt
        return np.concatenate([dcdt, dqdt])

    t_eval = np.linspace(0.0, t_service, 201)
    sol = solve_ivp(
        rhs,
        (0.0, t_service),
        y0,
        method="BDF",
        t_eval=t_eval,
        rtol=1e-6,
        atol=1e-9,
    )
    assert sol.success, sol.message
    return float(np.max(sol.y[n_nodes - 1, :]))


@pytest.fixture(scope="module")
def submission():
    return _load_json(ARTIFACT)


@pytest.fixture(scope="module")
def ground_truth():
    return _load_json(GT_PATH)


@pytest.fixture(scope="module")
def specs():
    return _load_json(SPECS_PATH)


def test_artifact_schema(submission):
    for key in REQUIRED_KEYS:
        assert key in submission, f"Missing key {key}"
    for key in REQUIRED_KEYS:
        val = _finite_number(submission[key], key)
        assert val > 0.0, f"{key} must be positive"


def test_parameter_accuracy(submission, ground_truth):
    truth = ground_truth["parameters"]
    tol = ground_truth["tolerances"]
    checks = [
        ("qm", tol["qm_rel"]),
        ("K", tol["K_rel"]),
        ("k_LDF", tol["k_LDF_rel"]),
        ("D_ax", tol["D_ax_rel"]),
    ]
    for key, rel_tol in checks:
        est = _finite_number(submission[key], key)
        ref = float(truth[key])
        assert abs(est - ref) / abs(ref) <= rel_tol + 1e-12, (
            f"{key} relative error too large: est={est}, ref={ref}, tol={rel_tol}"
        )


def test_design_length_accuracy(submission, ground_truth):
    est = _finite_number(submission["L_min_m"], "L_min_m")
    ref = float(ground_truth["design"]["L_min_m"])
    rel = float(ground_truth["tolerances"]["L_min_rel"])
    assert abs(est - ref) / abs(ref) <= rel + 1e-12, (
        f"L_min_m relative error too large: est={est}, ref={ref}, tol={rel}"
    )


def test_design_feasibility_under_truth(submission, ground_truth, specs):
    """Proposed L_min must keep true-model outlet under the duty limit."""
    L = _finite_number(submission["L_min_m"], "L_min_m")
    L_ref = float(ground_truth["design"]["L_min_m"])
    # Reject undersized beds even if relative-error rounding is favorable.
    assert L + 1e-12 >= 0.97 * L_ref, (
        f"L_min_m={L} is undersized vs reference {L_ref}"
    )
    duty = specs["design_duty"]
    truth = ground_truth["parameters"]
    c_in = float(duty["c_in_mol_m3"])
    c_limit = float(duty["outlet_mole_fraction_limit"]) * c_in
    max_out = simulate_outlet_max(
        L=L,
        u=float(duty["superficial_velocity_m_s"]),
        c_in=c_in,
        t_service=float(duty["service_time_s"]),
        qm=float(truth["qm"]),
        K=float(truth["K"]),
        k_ldf=float(truth["k_LDF"]),
        D_ax=float(truth["D_ax"]),
        eps=float(specs["bed_voidage"]),
        rho_p=float(specs["particle_density_kg_m3"]),
        n_nodes=56,
    )
    assert max_out <= c_limit * 1.03 + 1e-12, (
        f"True-model outlet {max_out} exceeds limit {c_limit} at proposed L"
    )
