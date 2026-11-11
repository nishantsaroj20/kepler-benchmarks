#!/usr/bin/env python3
"""Reference solver: staged kinetic ID + constrained feed optimization."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares, minimize

DATA = Path("/workspace/data")
OUT = Path("/workspace/output")
if not DATA.exists():
    ROOT = Path(__file__).resolve().parents[1]
    DATA = ROOT / "environment" / "data"
    OUT = ROOT / "output"


def load_specs():
    with open(DATA / "reactor_specs.json", encoding="utf-8") as f:
        return json.load(f)


def mu_qp(S, P, t, p, specs):
    mu_max = p["mu_max"] * (p["alpha"] if t >= p["t_switch"] else 1.0)
    S = max(S, 0.0)
    P = max(P, 0.0)
    denom = p["K_S"] + S + (S * S) / p["K_I"]
    gamma = float(specs["gamma"])
    inhib = max(1.0 - P / p["P_crit"], 0.0) ** gamma
    mu = mu_max * S / denom * inhib
    qp = specs["Y_px"] * mu if t < p["t_switch"] else p["beta"]
    return mu, qp


def rhs(t, y, feed_fn, p, specs):
    X = max(float(y[0]), 0.0)
    S = max(float(y[1]), 0.0)
    P = max(float(y[2]), 0.0)
    CL = max(float(y[3]), 0.0)
    V = float(y[4])
    F = max(float(feed_fn(t)), 0.0)
    V_eff = max(V, 1e-9)
    mu, qp = mu_qp(S, P, t, p, specs)
    q_o2 = specs["Y_xo"] * mu + specs["m_o"]
    return [
        mu * X - (F / V_eff) * X,
        -(1.0 / specs["Y_xs"]) * mu * X - specs["m_s"] * X + (F / V_eff) * (specs["S_in"] - S),
        qp * X - (F / V_eff) * P,
        specs["kLa"] * (specs["C_star"] - CL) - q_o2 * X - (F / V_eff) * CL,
        F,
    ]


def hold_feed(times, feeds):
    times = np.asarray(times, dtype=float)
    feeds = np.asarray(feeds, dtype=float)

    def feed_fn(t):
        if t <= times[0]:
            return float(feeds[0])
        if t >= times[-1]:
            return float(feeds[-1])
        i = int(np.searchsorted(times, t, side="right") - 1)
        i = min(max(i, 0), len(feeds) - 2)
        return float(feeds[i])

    return feed_fn


def piecewise_constant(edges, rates):
    edges = np.asarray(edges, dtype=float)
    rates = np.asarray(rates, dtype=float)

    def feed_fn(t):
        idx = int(np.searchsorted(edges, t, side="right") - 1)
        idx = min(max(idx, 0), len(rates) - 1)
        return float(rates[idx])

    return feed_fn


def simulate(y0, t_eval, feed_fn, p, specs, rtol=1e-6, atol=1e-8):
    t_eval = np.asarray(t_eval, dtype=float)
    sol = solve_ivp(
        lambda t, y: rhs(t, y, feed_fn, p, specs),
        (float(t_eval[0]), float(t_eval[-1])),
        np.asarray(y0, dtype=float),
        t_eval=t_eval,
        method="BDF",
        rtol=rtol,
        atol=atol,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y.T


def pack(vec):
    return {
        "mu_max": float(vec[0]),
        "K_S": float(vec[1]),
        "K_I": float(vec[2]),
        "P_crit": float(vec[3]),
        "t_switch": float(vec[4]),
        "alpha": float(vec[5]),
        "beta": float(vec[6]),
    }


def y0_from(df):
    return np.array(
        [
            max(float(df["X"].iloc[0]), 1e-6),
            max(float(df["S"].iloc[0]), 0.0),
            max(float(df["P"].iloc[0]), 0.0),
            max(float(df["C_L"].iloc[0]), 1e-6),
            max(float(df["V"].iloc[0]), 1e-3),
        ],
        dtype=float,
    )


def residuals(vec, dfs, specs, weights=None):
    if weights is None:
        weights = np.array([1.3, 1.0, 1.8, 0.3])
    p = pack(vec)
    chunks = []
    for df in dfs:
        t = df["time_h"].to_numpy(dtype=float)
        feed_fn = hold_feed(t, df["F"].to_numpy(dtype=float))
        try:
            Y = simulate(y0_from(df), t, feed_fn, p, specs, rtol=1e-5, atol=1e-7)
        except Exception:
            chunks.append(np.full(len(df) * 4, 1e3))
            continue
        scale = np.array(
            [
                max(float(np.std(df["X"])), 0.15),
                max(float(np.std(df["S"])), 0.2),
                max(float(np.std(df["P"])), 0.15),
                max(float(np.std(df["C_L"])), 4e-4),
            ]
        )
        obs = df[["X", "S", "P", "C_L"]].to_numpy(dtype=float)
        chunks.append((weights * (Y[:, :4] - obs) / scale).ravel())
    return np.concatenate(chunks)


def estimate_parameters(df1, df2, specs):
    dfs = [df1, df2]
    lb = np.array([0.25, 0.30, 8.0, 20.0, 9.0, 0.35, 0.04], dtype=float)
    ub = np.array([0.75, 2.00, 45.0, 45.0, 16.0, 0.70, 0.18], dtype=float)

    def sse(vec):
        r = residuals(vec, dfs, specs)
        return float(np.dot(r, r))

    # Coarse grid on the weakly identifiable pair (t_switch, P_crit)
    seeds = []
    for tsw in np.linspace(10.5, 13.5, 7):
        for pc in np.linspace(26.0, 34.0, 5):
            seeds.append(np.array([0.48, 0.85, 22.0, pc, tsw, 0.52, 0.10]))
    seeds.extend(
        [
            np.array([0.50, 0.80, 22.0, 30.0, 12.0, 0.50, 0.10]),
            np.array([0.45, 1.00, 25.0, 28.0, 11.5, 0.55, 0.09]),
            np.array([0.55, 0.70, 18.0, 32.0, 12.5, 0.45, 0.12]),
            np.array([0.52, 0.75, 20.0, 31.0, 12.0, 0.48, 0.11]),
            np.array([0.47, 0.90, 24.0, 29.0, 11.8, 0.53, 0.095]),
        ]
    )

    best = None
    best_val = 1e99
    for i, x0 in enumerate(seeds):
        res = least_squares(
            lambda v: residuals(v, dfs, specs),
            np.clip(x0, lb, ub),
            bounds=(lb, ub),
            max_nfev=55,
            xtol=1e-7,
            ftol=1e-7,
            method="trf",
        )
        val = sse(res.x)
        if val < best_val:
            best_val = val
            best = res.x
        if (i + 1) % 20 == 0:
            print(f"  grid progress {i+1}/{len(seeds)} best_sse={best_val:.4f}", flush=True)

    # Local polish around best
    rng = np.random.default_rng(7)
    candidates = [best]
    for _ in range(12):
        candidates.append(np.clip(best * (1.0 + rng.normal(0.0, 0.05, size=7)), lb, ub))

    w_prod = np.array([1.4, 1.1, 2.0, 0.25])
    for x0 in candidates:
        res = least_squares(
            lambda v: residuals(v, dfs, specs, weights=w_prod),
            x0,
            bounds=(lb, ub),
            max_nfev=180,
            xtol=1e-10,
            ftol=1e-10,
            method="trf",
        )
        val = sse(res.x)
        if val < best_val:
            best_val = val
            best = res.x

    # Refit K_S / alpha with biomass-heavy weights while keeping other params near best
    def resid_ks_alpha(z):
        vec = best.copy()
        vec[1] = z[0]  # K_S
        vec[5] = z[1]  # alpha
        return residuals(vec, dfs, specs, weights=np.array([2.0, 1.6, 1.0, 0.2]))

    z0 = np.array([best[1], best[5]])
    z_lb = np.array([lb[1], lb[5]])
    z_ub = np.array([ub[1], ub[5]])
    res_z = least_squares(
        resid_ks_alpha,
        z0,
        bounds=(z_lb, z_ub),
        max_nfev=120,
        xtol=1e-10,
        ftol=1e-10,
        method="trf",
    )
    trial = best.copy()
    trial[1] = res_z.x[0]
    trial[5] = res_z.x[1]
    # Final joint polish from that point
    res_final = least_squares(
        lambda v: residuals(v, dfs, specs, weights=np.array([1.6, 1.3, 1.7, 0.25])),
        trial,
        bounds=(lb, ub),
        max_nfev=200,
        xtol=1e-10,
        ftol=1e-10,
        method="trf",
    )
    if sse(res_final.x) <= best_val * 1.03:
        best = res_final.x
        best_val = sse(best)
    elif sse(trial) <= best_val * 1.03:
        best = trial
        best_val = sse(best)

    print(f"final sse={best_val:.6f}", flush=True)
    return pack(best)


def optimize_feed(p, specs):
    n = int(specs["n_feed_intervals"])
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
    Fmax = float(specs["max_pump_rate"])
    Vmax = float(specs["V_max"])
    CL_lim = 0.15 * float(specs["C_star"])
    S_lim = 15.0
    headroom = Vmax - y0[4]
    dt = tf / n

    def eval_rates(rates, n_grid=121):
        rates = np.clip(np.asarray(rates, dtype=float), 0.0, Fmax)
        vol = float(np.sum(rates) * dt)
        if vol > headroom:
            rates = rates * (headroom / (vol + 1e-12))
        feed_fn = piecewise_constant(edges[:-1], rates)
        t_eval = np.linspace(0.0, tf, n_grid)
        try:
            Y = simulate(y0, t_eval, feed_fn, p, specs)
        except Exception:
            return -1e9, rates
        V, S, CL, P = Y[:, 4], Y[:, 1], Y[:, 3], Y[:, 2]
        pen = (
            300.0 * np.sum(np.maximum(S - S_lim, 0.0) ** 2)
            + 900.0 * np.sum(np.maximum(CL_lim - CL, 0.0) ** 2)
            + 250.0 * np.sum(np.maximum(V - Vmax, 0.0) ** 2)
        )
        return float(V[-1] * P[-1]) - pen, rates

    starts = [
        np.array([0.05, 0.09, 0.14, 0.18, 0.16, 0.12, 0.10, 0.08, 0.06, 0.04]),
        np.array([0.04, 0.07, 0.11, 0.15, 0.17, 0.14, 0.11, 0.08, 0.05, 0.03]),
        np.full(n, 0.09),
        np.linspace(0.04, 0.14, n),
        np.array([0.03, 0.06, 0.10, 0.14, 0.18, 0.15, 0.11, 0.08, 0.05, 0.03]),
    ]
    best_score = -1e99
    best_rates = starts[0].copy()
    for x0 in starts:
        score, rates = eval_rates(x0)
        if score > best_score:
            best_score, best_rates = score, rates

    def objective(r):
        s, _ = eval_rates(r, n_grid=81)
        return -s

    bounds = [(0.0, Fmax)] * n
    for x0 in (best_rates, starts[0], starts[1]):
        res = minimize(
            objective,
            x0,
            method="SLSQP",
            bounds=bounds,
            options={"maxiter": 60, "ftol": 1e-9, "disp": False},
        )
        score, rates = eval_rates(res.x, n_grid=201)
        if score > best_score:
            best_score, best_rates = score, rates

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
        "n_intervals": n,
        "t_final": tf,
        "feed_rates_L_per_h": [float(x) for x in best_rates],
        "intervals": intervals,
    }


def main():
    specs = load_specs()
    df1 = pd.read_csv(DATA / "sensor_run_01.csv")
    df2 = pd.read_csv(DATA / "sensor_run_02.csv")
    print("Estimating parameters (grid + least squares)...", flush=True)
    params = estimate_parameters(df1, df2, specs)
    print(json.dumps(params, indent=2), flush=True)
    print("Optimizing feed...", flush=True)
    feed = optimize_feed(params, specs)
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "estimated_parameters.json", "w", encoding="utf-8") as f:
        json.dump(params, f, indent=2)
        f.write("\n")
    with open(OUT / "optimal_feed.json", "w", encoding="utf-8") as f:
        json.dump(feed, f, indent=2)
        f.write("\n")
    print("Wrote", OUT)


if __name__ == "__main__":
    main()
