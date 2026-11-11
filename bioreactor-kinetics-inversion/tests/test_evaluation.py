#!/usr/bin/env python3
"""Verifier for bioreactor-kinetics-inversion."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import solve_ivp

WORKSPACE = Path("/workspace")
OUTPUT = WORKSPACE / "output"
DATA = WORKSPACE / "data"
TESTS = Path("/tests")

if not (TESTS / "ground_truth_params.json").exists():
    ROOT = Path(__file__).resolve().parents[1]
    TESTS = ROOT / "tests"
    WORKSPACE = ROOT
    OUTPUT = ROOT / "output"
    DATA = ROOT / "environment" / "data"

GT_PATH = TESTS / "ground_truth_params.json"
PARAMS_PATH = OUTPUT / "estimated_parameters.json"
FEED_PATH = OUTPUT / "optimal_feed.json"
SPECS_PATH = DATA / "reactor_specs.json"

REQUIRED_PARAM_KEYS = [
    "mu_max",
    "K_S",
    "K_I",
    "P_crit",
    "t_switch",
    "alpha",
    "beta",
]


def _load_json(path: Path) -> dict:
    with open(path, encoding="utf-8") as f:
        return json.load(f)


@pytest.fixture(scope="module")
def gt() -> dict:
    return _load_json(GT_PATH)


@pytest.fixture(scope="module")
def specs() -> dict:
    return _load_json(SPECS_PATH)


def test_output_files_exist():
    assert PARAMS_PATH.is_file(), f"missing {PARAMS_PATH}"
    assert FEED_PATH.is_file(), f"missing {FEED_PATH}"


def test_estimated_parameters_schema_and_finite():
    data = _load_json(PARAMS_PATH)
    assert isinstance(data, dict)
    for key in REQUIRED_PARAM_KEYS:
        assert key in data, f"missing key {key}"
        val = data[key]
        assert isinstance(val, (int, float)), f"{key} must be numeric"
        assert math.isfinite(float(val)), f"{key} must be finite"
        assert float(val) >= 0.0, f"{key} must be non-negative"


def test_parameter_accuracy(gt):
    est = _load_json(PARAMS_PATH)
    truth = gt["parameters"]
    tol = gt["tolerances"]
    rel_primary = float(tol["relative_error_primary"])
    rel_haldane = float(tol["relative_error_haldane"])
    t_abs = float(tol["t_switch_abs_h"])
    primary_keys = tol.get("primary_keys", ["mu_max", "P_crit", "alpha", "beta"])
    haldane_keys = tol.get("haldane_keys", ["K_S", "K_I"])

    for key in primary_keys:
        tv = float(truth[key])
        ev = float(est[key])
        rel = abs(ev - tv) / max(abs(tv), 1e-12)
        assert rel <= rel_primary, f"{key} relative error {rel:.4f} exceeds {rel_primary}"

    for key in haldane_keys:
        tv = float(truth[key])
        ev = float(est[key])
        rel = abs(ev - tv) / max(abs(tv), 1e-12)
        assert rel <= rel_haldane, f"{key} relative error {rel:.4f} exceeds {rel_haldane}"

    assert abs(float(est["t_switch"]) - float(truth["t_switch"])) <= t_abs


def test_optimal_feed_schema(gt, specs):
    feed = _load_json(FEED_PATH)
    n = int(gt["tolerances"]["n_feed_intervals"])
    assert "feed_rates_L_per_h" in feed
    rates = feed["feed_rates_L_per_h"]
    assert isinstance(rates, list)
    assert len(rates) == n, f"feed profile must have length {n}"
    assert all(isinstance(x, (int, float)) for x in rates)
    assert all(math.isfinite(float(x)) for x in rates)
    assert all(float(x) >= 0.0 for x in rates)
    fmax = float(specs["max_pump_rate"])
    assert all(float(x) <= fmax + 1e-9 for x in rates)

    if "intervals" in feed:
        intervals = feed["intervals"]
        assert isinstance(intervals, list)
        assert len(intervals) == n
        for i, item in enumerate(intervals):
            assert int(item["index"]) == i
            assert math.isfinite(float(item["F"]))
            assert float(item["F"]) >= 0.0


def _mu_qp(S, P, t, p, specs):
    mu_max = p["mu_max"] * (p["alpha"] if t >= p["t_switch"] else 1.0)
    S = max(S, 0.0)
    P = max(P, 0.0)
    denom = p["K_S"] + S + (S * S) / p["K_I"]
    gamma = float(p.get("gamma", specs["gamma"]))
    inhib = max(1.0 - P / p["P_crit"], 0.0) ** gamma
    mu = mu_max * S / denom * inhib
    qp = specs["Y_px"] * mu if t < p["t_switch"] else p["beta"]
    return mu, qp


def _simulate_true(rates, gt, specs):
    p = gt["parameters"]
    n = len(rates)
    tf = float(specs["t_final"])
    edges = np.linspace(0.0, tf, n + 1)
    y0 = np.array(
        [
            specs["run03_initial_state"]["X0"],
            specs["run03_initial_state"]["S0"],
            specs["run03_initial_state"]["P0"],
            specs["run03_initial_state"]["CL0"],
            specs["run03_initial_state"]["V0"],
        ],
        dtype=float,
    )

    def feed_fn(t):
        idx = int(np.searchsorted(edges, t, side="right") - 1)
        idx = min(max(idx, 0), n - 1)
        return float(rates[idx])

    def rhs(t, y):
        X = max(float(y[0]), 0.0)
        S = max(float(y[1]), 0.0)
        P = max(float(y[2]), 0.0)
        CL = max(float(y[3]), 0.0)
        V = float(y[4])
        F = max(feed_fn(t), 0.0)
        V_eff = max(V, 1e-9)
        mu, qp = _mu_qp(S, P, t, p, specs)
        q_o2 = specs["Y_xo"] * mu + specs["m_o"]
        return [
            mu * X - (F / V_eff) * X,
            -(1.0 / specs["Y_xs"]) * mu * X
            - specs["m_s"] * X
            + (F / V_eff) * (specs["S_in"] - S),
            qp * X - (F / V_eff) * P,
            specs["kLa"] * (specs["C_star"] - CL) - q_o2 * X - (F / V_eff) * CL,
            F,
        ]

    t_eval = np.linspace(0.0, tf, 401)
    sol = solve_ivp(
        rhs,
        (0.0, tf),
        y0,
        t_eval=t_eval,
        method="Radau",
        rtol=1e-8,
        atol=1e-10,
    )
    assert sol.success, sol.message
    return sol.t, sol.y.T


def test_feed_constraints_and_performance(gt, specs):
    feed = _load_json(FEED_PATH)
    rates = [float(x) for x in feed["feed_rates_L_per_h"]]
    _t, Y = _simulate_true(rates, gt, specs)
    S = Y[:, 1]
    CL = Y[:, 3]
    V = Y[:, 4]
    P = Y[:, 2]

    S_lim = float(gt["tolerances"]["S_constraint_g_per_L"])
    S_slack = float(gt["tolerances"]["S_constraint_slack_g_per_L"])
    cl_frac = float(gt["tolerances"]["CL_fraction_of_Cstar"])
    CL_lim = cl_frac * float(specs["C_star"])
    Vmax = float(specs["V_max"])
    min_mass = float(gt["tolerances"]["min_product_mass_g"])

    assert float(np.max(S)) <= S_lim + S_slack
    assert float(np.min(CL)) >= CL_lim - 1e-4
    assert float(np.max(V)) <= Vmax + 1e-3
    product_mass = float(V[-1] * P[-1])
    assert product_mass + 1e-9 >= min_mass, (
        f"product mass {product_mass:.4f} below required {min_mass:.4f}"
    )
