#!/usr/bin/env python3
"""Independent reimplementation used for verification_explanation evidence.

This script does NOT import solution/solver.py. It re-simulates the true model
from tests/ground_truth_params.json, checks the reference feed feasibility, and
evaluates a naive baseline (constant mid-range feed + CSV-mean fake parameters)
to confirm it fails verifier thresholds.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[2]
GT = json.loads((ROOT / "tests" / "ground_truth_params.json").read_text())
SPECS = json.loads((ROOT / "environment" / "data" / "reactor_specs.json").read_text())
P = GT["parameters"]


def mu_qp(S, Pconc, t):
    mu_max = P["mu_max"] * (P["alpha"] if t >= P["t_switch"] else 1.0)
    S = max(S, 0.0)
    Pconc = max(Pconc, 0.0)
    denom = P["K_S"] + S + (S * S) / P["K_I"]
    inhib = max(1.0 - Pconc / P["P_crit"], 0.0) ** float(SPECS["gamma"])
    mu = mu_max * S / denom * inhib
    qp = SPECS["Y_px"] * mu if t < P["t_switch"] else P["beta"]
    return mu, qp


def simulate(rates):
    n = len(rates)
    tf = float(SPECS["t_final"])
    edges = np.linspace(0.0, tf, n + 1)
    y0 = np.array(
        [
            SPECS["run03_initial_state"]["X0"],
            SPECS["run03_initial_state"]["S0"],
            SPECS["run03_initial_state"]["P0"],
            SPECS["run03_initial_state"]["CL0"],
            SPECS["run03_initial_state"]["V0"],
        ],
        dtype=float,
    )

    def feed(t):
        idx = int(np.searchsorted(edges, t, side="right") - 1)
        idx = min(max(idx, 0), n - 1)
        return float(rates[idx])

    def rhs(t, y):
        X, S, Pc, CL, V = [max(float(y[i]), 0.0) for i in range(4)] + [float(y[4])]
        F = max(feed(t), 0.0)
        Ve = max(V, 1e-9)
        mu, qp = mu_qp(S, Pc, t)
        qo = SPECS["Y_xo"] * mu + SPECS["m_o"]
        return [
            mu * X - (F / Ve) * X,
            -(1 / SPECS["Y_xs"]) * mu * X - SPECS["m_s"] * X + (F / Ve) * (SPECS["S_in"] - S),
            qp * X - (F / Ve) * Pc,
            SPECS["kLa"] * (SPECS["C_star"] - CL) - qo * X - (F / Ve) * CL,
            F,
        ]

    t = np.linspace(0.0, tf, 401)
    sol = solve_ivp(rhs, (0.0, tf), y0, t_eval=t, method="Radau", rtol=1e-8, atol=1e-10)
    assert sol.success
    return sol.y.T


def main():
    ref = GT["reference_feed"]["feed_rates_L_per_h"]
    Y = simulate(ref)
    mass = float(Y[-1, 4] * Y[-1, 2])
    S_max = float(np.max(Y[:, 1]))
    CL_min = float(np.min(Y[:, 3]))
    print("reference_mass", mass)
    print("reference_S_max", S_max)
    print("reference_CL_min", CL_min)
    print("min_required_mass", GT["tolerances"]["min_product_mass_g"])
    print("reference_meets_mass", mass >= GT["tolerances"]["min_product_mass_g"])

    # Naive baseline: constant feed + parameter means from CSV scales
    naive_rates = [0.12] * 10
    Yn = simulate(naive_rates)
    naive_mass = float(Yn[-1, 4] * Yn[-1, 2])
    print("naive_constant_feed_mass", naive_mass)
    print(
        "naive_fails_mass_gate",
        naive_mass < GT["tolerances"]["min_product_mass_g"],
    )

    df1 = pd.read_csv(ROOT / "environment" / "data" / "sensor_run_01.csv")
    fake = {
        "mu_max": float(df1["X"].mean() / 20.0),
        "K_S": float(df1["S"].mean()),
        "K_I": 10.0,
        "P_crit": float(df1["P"].max() + 5.0),
        "t_switch": 10.0,
        "alpha": 0.7,
        "beta": 0.05,
    }
    print("fake_csv_mean_params", fake)
    for key in ("mu_max", "P_crit", "alpha", "beta"):
        rel = abs(fake[key] - P[key]) / abs(P[key])
        print(f"fake_{key}_rel", rel, "fails", rel > GT["tolerances"]["relative_error_primary"])


if __name__ == "__main__":
    main()
