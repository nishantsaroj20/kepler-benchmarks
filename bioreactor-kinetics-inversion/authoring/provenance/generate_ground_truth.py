#!/usr/bin/env python3
"""Seeded generator for bioreactor-kinetics-inversion ground truth and telemetry."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[2]
PROV = Path(__file__).resolve().parent
TESTS = ROOT / "tests"
DATA = ROOT / "environment" / "data"

SEED = 20241111
RNG = np.random.default_rng(SEED)

GT = {
    "mu_max": 0.50,
    "K_S": 0.80,
    "K_I": 22.0,
    "P_crit": 30.0,
    "t_switch": 12.0,
    "alpha": 0.50,
    "beta": 0.10,
}

SPECS = {
    "Y_xs": 0.47,
    "Y_px": 0.35,
    "m_s": 0.022,
    "Y_xo": 0.28,
    "m_o": 0.006,
    "kLa": 360.0,
    "C_star": 0.008,
    "S_in": 60.0,
    "V_max": 10.0,
    "max_pump_rate": 0.45,
    "t_final": 30.0,
    "n_feed_intervals": 10,
    "gamma": 1.0,
    "run03_initial_state": {
        "X0": 1.0,
        "S0": 7.0,
        "P0": 0.3,
        "CL0": 0.0072,
        "V0": 3.5,
    },
    "units": {
        "time": "h",
        "X": "g/L",
        "S": "g/L",
        "P": "g/L",
        "C_L": "g/L",
        "V": "L",
        "F": "L/h",
    },
}


def full_params(p: dict | None = None) -> dict:
    base = dict(GT if p is None else p)
    base["gamma"] = float(SPECS["gamma"])
    return base


def mu_qp(S: float, P: float, t: float, p: dict) -> tuple[float, float]:
    mu_max = p["mu_max"] * (p["alpha"] if t >= p["t_switch"] else 1.0)
    S = max(S, 0.0)
    P = max(P, 0.0)
    denom = p["K_S"] + S + (S * S) / p["K_I"]
    inhib = max(1.0 - P / p["P_crit"], 0.0) ** p["gamma"]
    mu = mu_max * S / denom * inhib
    qp = SPECS["Y_px"] * mu if t < p["t_switch"] else p["beta"]
    return mu, qp


def rhs(t, y, feed_fn, p):
    X = max(float(y[0]), 0.0)
    S = max(float(y[1]), 0.0)
    P = max(float(y[2]), 0.0)
    CL = max(float(y[3]), 0.0)
    V = float(y[4])
    F = max(float(feed_fn(t)), 0.0)
    V_eff = max(V, 1e-9)
    mu, qp = mu_qp(S, P, t, p)
    q_o2 = SPECS["Y_xo"] * mu + SPECS["m_o"]
    return [
        mu * X - (F / V_eff) * X,
        -(1.0 / SPECS["Y_xs"]) * mu * X - SPECS["m_s"] * X + (F / V_eff) * (SPECS["S_in"] - S),
        qp * X - (F / V_eff) * P,
        SPECS["kLa"] * (SPECS["C_star"] - CL) - q_o2 * X - (F / V_eff) * CL,
        F,
    ]


def piecewise_feed(times, rates):
    times = np.asarray(times, dtype=float)
    rates = np.asarray(rates, dtype=float)

    def feed_fn(t):
        idx = int(np.searchsorted(times, t, side="right") - 1)
        idx = min(max(idx, 0), len(rates) - 1)
        return float(rates[idx])

    return feed_fn


def simulate(y0, t_eval, feed_fn, p):
    t_eval = np.asarray(t_eval, dtype=float)
    sol = solve_ivp(
        lambda t, y: rhs(t, y, feed_fn, p),
        (float(t_eval[0]), float(t_eval[-1])),
        np.asarray(y0, dtype=float),
        t_eval=t_eval,
        method="Radau",
        rtol=1e-8,
        atol=1e-10,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y.T


def sample_times(t0, tf, n, rng):
    n_early = max(n // 3, 8)
    n_mid = max(n // 3, 8)
    n_late = max(n - n_early - n_mid, 8)
    t_early = np.linspace(t0, 8.0, n_early, endpoint=False)
    t_mid = np.sort(rng.uniform(8.0, 16.0, size=n_mid))
    t_late = np.linspace(16.0, tf, n_late)
    times = np.unique(np.concatenate([[t0], t_early, t_mid, t_late, [tf]]))
    kept = [times[0]]
    for t in times[1:]:
        if t - kept[-1] >= 0.18:
            kept.append(float(t))
    kept[-1] = float(tf)
    return np.asarray(kept, dtype=float)


def add_noise(df, rng):
    out = df.copy()
    specs = {
        "X": (0.018, 0.015),
        "S": (0.020, 0.02),
        "P": (0.018, 0.015),
        "C_L": (0.040, 0.0001),
        "V": (0.004, 0.006),
        "F": (0.010, 0.001),
    }
    for col, (rel, abs_) in specs.items():
        noise = rng.normal(0.0, 1.0, size=len(out)) * (
            rel * np.abs(out[col].to_numpy()) + abs_
        )
        out[col] = out[col].to_numpy() + noise
    for col in ["X", "S", "P", "C_L", "V", "F"]:
        out.loc[out.index[0], col] = df.loc[df.index[0], col]
        out[col] = np.clip(out[col], 0.0, None)
    return out


def build_run(run_id, y0, bp, rates, t_final, n_samples, rng):
    p = full_params()
    feed_fn = piecewise_feed(bp, rates)
    t_dense = np.linspace(0.0, t_final, 801)
    Y = simulate(y0, t_dense, feed_fn, p)
    t_samp = sample_times(0.0, t_final, n_samples, rng)
    Ys = np.vstack([np.interp(t_samp, t_dense, Y[:, i]) for i in range(5)]).T
    F_true = np.array([feed_fn(t) for t in t_samp])
    df = pd.DataFrame(
        {
            "time_h": t_samp,
            "X": Ys[:, 0],
            "S": Ys[:, 1],
            "P": Ys[:, 2],
            "C_L": Ys[:, 3],
            "V": Ys[:, 4],
            "F": F_true,
        }
    )
    df = add_noise(df, rng)
    df.insert(0, "run_id", run_id)
    return df


def optimize_reference_feed(p_in):
    p = full_params(p_in)
    n = SPECS["n_feed_intervals"]
    tf = SPECS["t_final"]
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
    Fmax = SPECS["max_pump_rate"]
    Vmax = SPECS["V_max"]
    CL_lim = 0.15 * SPECS["C_star"]
    S_lim = 15.0
    headroom = Vmax - y0[4]

    def eval_feed(rates, n_grid=121):
        rates = np.clip(np.asarray(rates, dtype=float), 0.0, Fmax)
        dt = tf / n
        if float(np.sum(rates) * dt) > headroom:
            rates = rates * (headroom / (float(np.sum(rates) * dt) + 1e-12))
        feed_fn = piecewise_feed(edges[:-1], rates)
        t_eval = np.linspace(0.0, tf, n_grid)
        try:
            Y = simulate(y0, t_eval, feed_fn, p)
        except RuntimeError:
            return -1e9, None
        V, S, CL, P = Y[:, 4], Y[:, 1], Y[:, 3], Y[:, 2]
        pen = (
            250.0 * np.sum(np.maximum(S - S_lim, 0.0) ** 2)
            + 800.0 * np.sum(np.maximum(CL_lim - CL, 0.0) ** 2)
            + 200.0 * np.sum(np.maximum(V - Vmax, 0.0) ** 2)
        )
        mass = float(V[-1] * P[-1])
        return mass - pen, {
            "product_mass": mass,
            "V_final": float(V[-1]),
            "P_final": float(P[-1]),
            "S_max": float(np.max(S)),
            "CL_min": float(np.min(CL)),
            "feasible_hard": bool(
                np.all(S <= S_lim + 0.25)
                and np.all(CL >= CL_lim - 1e-4)
                and np.all(V <= Vmax + 1e-6)
            ),
            "feed_rates_used": rates.tolist(),
        }

    candidates = [
        np.array([0.05, 0.09, 0.14, 0.18, 0.16, 0.12, 0.10, 0.08, 0.06, 0.04]),
        np.array([0.04, 0.07, 0.11, 0.15, 0.17, 0.14, 0.11, 0.08, 0.05, 0.03]),
        np.array([0.06, 0.10, 0.13, 0.15, 0.13, 0.11, 0.09, 0.07, 0.05, 0.03]),
        np.linspace(0.04, 0.14, n),
        np.full(n, 0.09),
    ]
    best_val = -1e99
    best_rates = candidates[0].copy()
    for x0 in candidates:
        val, meta = eval_feed(x0)
        if meta is not None and val > best_val:
            best_val = val
            best_rates = np.array(meta["feed_rates_used"], dtype=float)

    def objective(r):
        v, _ = eval_feed(r, n_grid=81)
        return -v

    res = minimize(
        objective,
        best_rates,
        method="SLSQP",
        bounds=[(0.0, Fmax)] * n,
        options={"maxiter": 45, "ftol": 1e-8, "disp": False},
    )
    val, meta = eval_feed(res.x, n_grid=201)
    if meta is not None and val > best_val:
        best_rates = np.array(meta["feed_rates_used"], dtype=float)

    _, final_meta = eval_feed(best_rates, n_grid=301)
    best_rates = np.array(final_meta["feed_rates_used"], dtype=float)
    intervals = [
        {
            "index": i,
            "t_start": float(edges[i]),
            "t_end": float(edges[i + 1]),
            "F": float(best_rates[i]),
        }
        for i in range(n)
    ]
    return {
        "interval_hours": float(tf / n),
        "feed_rates_L_per_h": [float(x) for x in best_rates],
        "intervals": intervals,
        "reference_metrics": {
            k: final_meta[k]
            for k in (
                "product_mass",
                "V_final",
                "P_final",
                "S_max",
                "CL_min",
                "feasible_hard",
            )
        },
    }


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    TESTS.mkdir(parents=True, exist_ok=True)

    # Run 01 starts with a short batch (F=0) so S declines through K_S.
    y0_1 = np.array([1.0, 9.0, 0.2, 0.007, 3.2], dtype=float)
    bp1 = np.array([0.0, 3.0, 7.0, 11.0, 15.0, 20.0, 24.0])
    rates1 = np.array([0.0, 0.12, 0.18, 0.09, 0.05, 0.04, 0.02])
    df1 = build_run("run01", y0_1, bp1, rates1, t_final=28.0, n_samples=44, rng=RNG)

    y0_2 = np.array([1.2, 4.5, 0.5, 0.007, 4.0], dtype=float)
    bp2 = np.array([0.0, 3.5, 7.0, 11.0, 15.0, 19.0, 24.0, 28.0])
    rates2 = np.array([0.04, 0.11, 0.17, 0.09, 0.06, 0.10, 0.04, 0.02])
    df2 = build_run("run02", y0_2, bp2, rates2, t_final=30.0, n_samples=46, rng=RNG)

    df1.to_csv(DATA / "sensor_run_01.csv", index=False, float_format="%.8g")
    df2.to_csv(DATA / "sensor_run_02.csv", index=False, float_format="%.8g")
    with open(DATA / "reactor_specs.json", "w", encoding="utf-8") as f:
        json.dump(SPECS, f, indent=2)
        f.write("\n")

    ref_feed = optimize_reference_feed(GT)
    ref_mass = float(ref_feed["reference_metrics"]["product_mass"])
    min_product_mass = 0.90 * ref_mass

    gt_payload = {
        "seed": SEED,
        "integrator": "scipy.integrate.solve_ivp method=Radau rtol=1e-8 atol=1e-10",
        "parameters": full_params(GT),
        "graded_keys": [
            "mu_max",
            "K_S",
            "K_I",
            "P_crit",
            "t_switch",
            "alpha",
            "beta",
        ],
        "known_constants_in_specs": ["gamma"],
        "tolerances": {
            "relative_error_primary": 0.10,
            "relative_error_haldane": 0.20,
            "primary_keys": ["mu_max", "P_crit", "alpha", "beta"],
            "haldane_keys": ["K_S", "K_I"],
            "t_switch_abs_h": 0.5,
            "S_constraint_g_per_L": 15.0,
            "S_constraint_slack_g_per_L": 0.5,
            "CL_fraction_of_Cstar": 0.15,
            "min_product_mass_g": min_product_mass,
            "n_feed_intervals": 10,
        },
        "reference_feed": ref_feed,
        "identification_runs": {
            "run01": {
                "y0": y0_1.tolist(),
                "feed_breakpoints_h": bp1.tolist(),
                "feed_rates_L_per_h": rates1.tolist(),
                "t_final_h": 28.0,
                "n_samples": int(len(df1)),
            },
            "run02": {
                "y0": y0_2.tolist(),
                "feed_breakpoints_h": bp2.tolist(),
                "feed_rates_L_per_h": rates2.tolist(),
                "t_final_h": 30.0,
                "n_samples": int(len(df2)),
            },
        },
    }
    for dest in (PROV / "ground_truth_params.json", TESTS / "ground_truth_params.json"):
        with open(dest, "w", encoding="utf-8") as f:
            json.dump(gt_payload, f, indent=2)
            f.write("\n")

    print("Pmax run01/run02:", float(df1["P"].max()), float(df2["P"].max()))
    print("Reference product mass:", ref_mass)
    print("Min product mass:", min_product_mass)
    print("Reference metrics:", json.dumps(ref_feed["reference_metrics"], indent=2))


if __name__ == "__main__":
    main()
