#!/usr/bin/env python3
"""Estimate Langmuir+LDF+D_ax from two breakthroughs and size a packed bed."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

DATA_ROOT = Path("/workspace/data")
OUT_ROOT = Path("/workspace/output")
if not DATA_ROOT.exists():
    DATA_ROOT = Path(__file__).resolve().parents[1] / "environment" / "data"
    OUT_ROOT = Path(__file__).resolve().parents[1] / "output"

N_FIT = 48
N_DESIGN = 48


def q_star(c, qm, K):
    return qm * K * c / (1.0 + K * c)


def simulate_ports(
    *,
    L,
    u,
    c_in,
    t_obs,
    qm,
    K,
    k_ldf,
    D_ax,
    eps,
    rho_p,
    n_nodes=N_FIT,
):
    dz = L / n_nodes
    factor = ((1.0 - eps) / eps) * rho_p
    y0 = np.zeros(2 * n_nodes)
    mid = n_nodes // 2
    t_obs = np.asarray(t_obs, dtype=float)
    t_end = float(np.max(t_obs)) if t_obs.size else 0.0
    if t_end <= 0.0:
        z = np.zeros_like(t_obs)
        return z, z

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

    t_eval = np.unique(np.concatenate(([0.0], t_obs)))
    sol = solve_ivp(
        rhs,
        (0.0, t_end),
        y0,
        method="BDF",
        t_eval=t_eval,
        rtol=1e-5,
        atol=1e-8,
    )
    if not sol.success:
        bad = np.full(t_obs.shape, 1e3)
        return bad, bad
    c_out = np.interp(t_obs, sol.t, sol.y[n_nodes - 1, :])
    c_mid = np.interp(t_obs, sol.t, sol.y[mid, :])
    return c_out, c_mid


def load_inputs():
    specs = json.loads((DATA_ROOT / "bed_specs.json").read_text())
    exp_a = pd.read_csv(DATA_ROOT / "breakthrough_exp_A.csv")
    exp_b = pd.read_csv(DATA_ROOT / "breakthrough_exp_B.csv")
    return specs, exp_a, exp_b


def pack4(qm, K, k_ldf, D_ax):
    return np.log(np.array([qm, K, k_ldf, D_ax], dtype=float))


def unpack4(theta):
    qm, K, k_ldf, D_ax = np.exp(theta)
    return float(qm), float(K), float(k_ldf), float(D_ax)


def pack2(k_ldf, D_ax):
    return np.log(np.array([k_ldf, D_ax], dtype=float))


def unpack2(theta):
    k_ldf, D_ax = np.exp(theta)
    return float(k_ldf), float(D_ax)


def isotherm_seed(specs, exp_a, exp_b):
    """Crude capacity / curvature seeds from stoichiometric uptake estimates."""
    eps = specs["bed_voidage"]
    rho = specs["particle_density_kg_m3"]
    L = specs["lab_bed_length_m"]
    qs = []
    cs = []
    for df, key in ((exp_a, "experiment_A"), (exp_b, "experiment_B")):
        meta = specs[key]
        c_in = float(meta["c_in_mol_m3"])
        u = float(meta["superficial_velocity_m_s"])
        t = df["time_s"].to_numpy(dtype=float)
        c = np.clip(df["c_outlet_mol_m3"].to_numpy(dtype=float), 0.0, c_in)
        order = np.argsort(t)
        t = t[order]
        c = c[order]
        # Approximate unused fraction integral for solid loading.
        unused = np.trapezoid(np.clip(c_in - c, 0.0, None), t)
        q_est = (eps / ((1.0 - eps) * rho)) * (u / L) * unused
        qs.append(max(q_est, 1e-4))
        cs.append(c_in)
    c1, c2 = cs
    q1, q2 = qs
    # Two-point Langmuir inversion: q = qm K c / (1 + K c)
    # q/c = qm K - K q  =>  linearize in (q, q/c)
    # From q1,q2,c1,c2:
    # 1/q = 1/qm + 1/(qm K c)
    a = np.array([[1.0, 1.0 / c1], [1.0, 1.0 / c2]], dtype=float)
    b = np.array([1.0 / q1, 1.0 / q2], dtype=float)
    try:
        x = np.linalg.solve(a, b)
        qm = 1.0 / x[0]
        K = 1.0 / (qm * x[1])
    except np.linalg.LinAlgError:
        qm, K = 0.4, 2.5
    qm = float(np.clip(qm, 0.1, 1.2))
    K = float(np.clip(K, 0.4, 10.0))
    return qm, K


def residual_all(theta, specs, exp_a, exp_b):
    qm, K, k_ldf, D_ax = unpack4(theta)
    return _port_residuals(specs, exp_a, exp_b, qm, K, k_ldf, D_ax)


def residual_transport(theta, specs, exp_a, exp_b, qm, K):
    k_ldf, D_ax = unpack2(theta)
    return _port_residuals(specs, exp_a, exp_b, qm, K, k_ldf, D_ax)


def _port_residuals(specs, exp_a, exp_b, qm, K, k_ldf, D_ax):
    eps = specs["bed_voidage"]
    rho_p = specs["particle_density_kg_m3"]
    L = specs["lab_bed_length_m"]
    chunks = []
    for df, meta_key in ((exp_a, "experiment_A"), (exp_b, "experiment_B")):
        meta = specs[meta_key]
        t = df["time_s"].to_numpy(dtype=float)
        y_out = df["c_outlet_mol_m3"].to_numpy(dtype=float)
        y_mid = df["c_midbed_mol_m3"].to_numpy(dtype=float)
        pred_out, pred_mid = simulate_ports(
            L=L,
            u=meta["superficial_velocity_m_s"],
            c_in=meta["c_in_mol_m3"],
            t_obs=t,
            qm=qm,
            K=K,
            k_ldf=k_ldf,
            D_ax=D_ax,
            eps=eps,
            rho_p=rho_p,
            n_nodes=N_FIT,
        )
        scale = max(float(meta["c_in_mol_m3"]), 1e-6)
        # Emphasize front region where transport parameters matter most.
        front = np.exp(-(((y_out / scale) - 0.5) ** 2) / (2 * 0.35**2))
        w = 0.55 + 0.90 * front
        chunks.append(w * (pred_out - y_out) / scale)
        chunks.append(0.85 * w * (pred_mid - y_mid) / scale)
    return np.concatenate(chunks)


def estimate_parameters(specs, exp_a, exp_b):
    qm0, K0 = isotherm_seed(specs, exp_a, exp_b)

    # Stage 1: transport only on a grid of isotherm seeds near the capacity guess.
    transport_starts = [
        pack2(0.02, 1.5e-4),
        pack2(0.03, 2.5e-4),
        pack2(0.04, 3.5e-4),
        pack2(0.025, 2.0e-4),
    ]
    isotherm_seeds = [
        (qm0, K0),
        (qm0 * 0.92, K0 * 1.1),
        (qm0 * 1.08, K0 * 0.9),
    ]
    best = None
    for qm_s, K_s in isotherm_seeds:
        for th0 in transport_starts:
            res = least_squares(
                residual_transport,
                th0,
                args=(specs, exp_a, exp_b, qm_s, K_s),
                bounds=(pack2(0.005, 2.0e-5), pack2(0.2, 2.0e-3)),
                xtol=1e-8,
                ftol=1e-8,
                gtol=1e-8,
                max_nfev=35,
            )
            cost = float(np.sum(res.fun**2))
            k_ldf, D_ax = unpack2(res.x)
            cand = (cost, qm_s, K_s, k_ldf, D_ax)
            if best is None or cost < best[0]:
                best = cand

    # Stage 2: joint polish from the best staged candidate.
    theta0 = pack4(best[1], best[2], best[3], best[4])
    polish_starts = [
        theta0,
        pack4(best[1], best[2], best[3] * 1.15, best[4] * 0.85),
        pack4(best[1] * 1.05, best[2] * 0.95, best[3] * 0.9, best[4] * 1.1),
    ]
    best4 = None
    for th0 in polish_starts:
        res = least_squares(
            residual_all,
            th0,
            args=(specs, exp_a, exp_b),
            bounds=(
                pack4(0.08, 0.4, 0.005, 2.0e-5),
                pack4(1.2, 12.0, 0.25, 2.5e-3),
            ),
            xtol=1e-8,
            ftol=1e-8,
            gtol=1e-8,
            max_nfev=50,
        )
        cost = float(np.sum(res.fun**2))
        if best4 is None or cost < best4[0]:
            best4 = (cost, res.x)

    qm, K, k_ldf, D_ax = unpack4(best4[1])
    return {"qm": qm, "K": K, "k_LDF": k_ldf, "D_ax": D_ax, "fit_cost": best4[0]}


def design_length(params, specs):
    duty = specs["design_duty"]
    eps = specs["bed_voidage"]
    rho_p = specs["particle_density_kg_m3"]
    c_in = float(duty["c_in_mol_m3"])
    u = float(duty["superficial_velocity_m_s"])
    t_service = float(duty["service_time_s"])
    c_limit = float(duty["outlet_mole_fraction_limit"]) * c_in
    t_obs = np.linspace(0.0, t_service, 121)

    def max_outlet(L):
        c_out, _ = simulate_ports(
            L=L,
            u=u,
            c_in=c_in,
            t_obs=t_obs,
            qm=params["qm"],
            K=params["K"],
            k_ldf=params["k_LDF"],
            D_ax=params["D_ax"],
            eps=eps,
            rho_p=rho_p,
            n_nodes=N_DESIGN,
        )
        return float(np.max(c_out))

    lo, hi = 0.05, 0.12
    if max_outlet(lo) <= c_limit:
        while lo > 1e-3 and max_outlet(lo) <= c_limit:
            hi = lo
            lo *= 0.5
    else:
        while max_outlet(hi) > c_limit:
            hi *= 1.55
            if hi > 25.0:
                raise RuntimeError("Unable to find feasible packed length")
        lo = hi / 1.55

    for _ in range(26):
        mid = 0.5 * (lo + hi)
        if max_outlet(mid) > c_limit:
            lo = mid
        else:
            hi = mid
    # Small safety margin against grid/noise bias in the inverse design step.
    return float(hi * 1.04)


def main():
    specs, exp_a, exp_b = load_inputs()
    params = estimate_parameters(specs, exp_a, exp_b)
    L_min = design_length(params, specs)
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    payload = {
        "qm": params["qm"],
        "K": params["K"],
        "k_LDF": params["k_LDF"],
        "D_ax": params["D_ax"],
        "L_min_m": L_min,
    }
    (OUT_ROOT / "adsorption_identification.json").write_text(
        json.dumps(payload, indent=2) + "\n"
    )
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
