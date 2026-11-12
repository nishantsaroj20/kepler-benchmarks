#!/usr/bin/env python3
"""Reference column solver: identify Langmuir+LDF+D_ax and size L_min.

CLI:
  python3 column_solver.py --case /path/to/case.json --out /path/to/result.json
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

N_FIT = 32
N_DESIGN = 36

# NumPy 1.x uses trapz; 2.x uses trapezoid.
_trapz = getattr(np, "trapezoid", None) or np.trapz


def q_star(c, qm, K):
    return qm * K * c / (1.0 + K * c)


def simulate_outlet(
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
    also_mid=False,
):
    dz = L / n_nodes
    factor = ((1.0 - eps) / eps) * rho_p
    y0 = np.zeros(2 * n_nodes)
    mid = n_nodes // 2
    t_obs = np.asarray(t_obs, dtype=float)
    t_end = float(np.max(t_obs)) if t_obs.size else 0.0
    if t_end <= 0.0:
        z = np.zeros_like(t_obs)
        return (z, z) if also_mid else z

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
        rtol=3e-5,
        atol=1e-8,
    )
    if not sol.success:
        bad = np.full(t_obs.shape, 1e3)
        return (bad, bad) if also_mid else bad
    c_out = np.interp(t_obs, sol.t, sol.y[n_nodes - 1, :])
    if not also_mid:
        return c_out
    c_mid = np.interp(t_obs, sol.t, sol.y[mid, :])
    return c_out, c_mid


def pack4(qm, K, k_ldf, D_ax):
    return np.log(np.array([qm, K, k_ldf, D_ax], dtype=float))


def unpack4(theta):
    qm, K, k_ldf, D_ax = np.exp(theta)
    return float(qm), float(K), float(k_ldf), float(D_ax)


def load_case(case_path: Path):
    case = json.loads(case_path.read_text())
    case_dir = case_path.parent
    frames = {}
    for key in ("experiment_A", "experiment_B"):
        meta = case[key]
        frames[key] = pd.read_csv(case_dir / meta["csv"])
    return case, frames


def uptake_qstar(case, frames):
    """PDE-consistent completed-breakthrough solid loadings."""
    eps = case["bed_voidage"]
    rho = case["particle_density_kg_m3"]
    L = case["lab_bed_length_m"]
    qs, cs = [], []
    for key in ("experiment_A", "experiment_B"):
        meta = case[key]
        df = frames[key]
        c_in = float(meta["c_in_mol_m3"])
        u = float(meta["interstitial_velocity_m_s"])
        t = df["time_s"].to_numpy(dtype=float)
        c = df["c_outlet_mol_m3"].to_numpy(dtype=float)
        order = np.argsort(t)
        t, c = t[order], c[order]
        unused = _trapz(np.clip(c_in - c, 0.0, None), t)
        rhs = u * unused / L - c_in
        q_est = rhs * eps / ((1.0 - eps) * rho)
        qs.append(max(float(q_est), 1e-4))
        cs.append(c_in)
    return np.array(cs), np.array(qs)


def isotherm_seed(case, frames):
    cs, qs = uptake_qstar(case, frames)
    c1, c2 = float(cs[0]), float(cs[1])
    q1, q2 = float(qs[0]), float(qs[1])
    denom = q2 - q1
    if abs(denom) < 1e-12:
        qm, K = max(q1, q2) * 1.25, 1.0
    else:
        K = (q1 / c1 - q2 / c2) / denom
        if K <= 1e-6:
            a = np.array([[1.0, 1.0 / c1], [1.0, 1.0 / c2]], dtype=float)
            b = np.array([1.0 / q1, 1.0 / q2], dtype=float)
            try:
                x = np.linalg.solve(a, b)
                qm = 1.0 / x[0]
                K = 1.0 / (qm * x[1])
            except np.linalg.LinAlgError:
                qm, K = max(q1, q2) * 1.25, 1.0
        else:
            qm = q1 * (1.0 + K * c1) / (K * c1)
    return float(np.clip(qm, 0.08, 2.5)), float(np.clip(K, 0.08, 20.0))


def clean_series(df):
    t = df["time_s"].to_numpy(dtype=float)
    y_out = df["c_outlet_mol_m3"].to_numpy(dtype=float)
    y_mid = df["c_midbed_mol_m3"].to_numpy(dtype=float)
    finite = np.isfinite(y_mid)
    c_ref = max(float(np.nanmax(y_out)), 1e-6)
    ok = finite & (y_mid >= -0.02 * c_ref) & (y_mid <= 1.25 * c_ref)
    mid = (t[ok], y_mid[ok]) if ok.sum() >= max(8, len(y_mid) // 3) else None
    return t, y_out, mid


def outlet_residuals(case, frames, qm, K, k_ldf, D_ax, use_mid=True):
    eps = case["bed_voidage"]
    rho_p = case["particle_density_kg_m3"]
    L = case["lab_bed_length_m"]
    chunks = []
    for key in ("experiment_A", "experiment_B"):
        meta = case[key]
        t, y_out, mid = clean_series(frames[key])
        pred = simulate_outlet(
            L=L,
            u=meta["interstitial_velocity_m_s"],
            c_in=meta["c_in_mol_m3"],
            t_obs=t,
            qm=qm,
            K=K,
            k_ldf=k_ldf,
            D_ax=D_ax,
            eps=eps,
            rho_p=rho_p,
            n_nodes=N_FIT,
            also_mid=use_mid and mid is not None,
        )
        if use_mid and mid is not None:
            pred_out, pred_mid_full = pred
        else:
            pred_out = pred
            pred_mid_full = None
        scale = max(float(meta["c_in_mol_m3"]), 1e-6)
        # Emphasize the front; downweight saturated late points.
        frac = np.clip(y_out / scale, 0.0, 1.2)
        w = 0.35 + 1.1 * np.exp(-((frac - 0.45) ** 2) / (2 * 0.28**2))
        chunks.append(w * (pred_out - y_out) / scale)
        if pred_mid_full is not None:
            t_m, y_m = mid
            pred_m = np.interp(t_m, t, pred_mid_full)
            frac_m = np.clip(y_m / scale, 0.0, 1.2)
            w_m = 0.35 + 1.1 * np.exp(-((frac_m - 0.45) ** 2) / (2 * 0.28**2))
            chunks.append(0.55 * w_m * (pred_m - y_m) / scale)
    return np.concatenate(chunks)


def estimate_parameters(case, frames):
    qm0, K0 = isotherm_seed(case, frames)
    # Multi-start joint fit in log space around the uptake seed.
    starts = []
    for qm in (qm0, qm0 * 0.88, qm0 * 1.12):
        for K in (K0, K0 * 0.7, K0 * 1.35, max(0.12, K0 * 0.4)):
            for k_ldf, D_ax in (
                (0.016, 1.2e-4),
                (0.022, 2.5e-4),
                (0.032, 4.5e-4),
                (0.045, 8.0e-4),
                (0.060, 1.2e-3),
            ):
                starts.append(pack4(qm, max(K, 0.1), k_ldf, D_ax))
    # Cap starts for runtime (~1 CPU, per-run limit) while keeping diversity.
    starts = starts[::5][:12]

    best = None
    for th0 in starts:
        res = least_squares(
            lambda th: outlet_residuals(case, frames, *unpack4(th), use_mid=True),
            th0,
            bounds=(
                pack4(0.08, 0.08, 0.005, 3.0e-5),
                pack4(2.5, 20.0, 0.25, 2.5e-3),
            ),
            xtol=3e-7,
            ftol=3e-7,
            gtol=3e-7,
            max_nfev=30,
        )
        cost = float(np.sum(res.fun**2))
        if best is None or cost < best[0]:
            best = (cost, res.x)

    # Polish best with a bit more budget (outlet-heavy).
    res = least_squares(
        lambda th: outlet_residuals(case, frames, *unpack4(th), use_mid=True),
        best[1],
        bounds=(
            pack4(0.08, 0.08, 0.005, 3.0e-5),
            pack4(2.5, 20.0, 0.25, 2.5e-3),
        ),
        xtol=1e-8,
        ftol=1e-8,
        gtol=1e-8,
        max_nfev=50,
    )
    qm, K, k_ldf, D_ax = unpack4(res.x)
    return {"qm": qm, "K": K, "k_LDF": k_ldf, "D_ax": D_ax}


def design_length(params, case):
    duty = case["design_duty"]
    eps = case["bed_voidage"]
    rho_p = case["particle_density_kg_m3"]
    c_in = float(duty["c_in_mol_m3"])
    u = float(duty["interstitial_velocity_m_s"])
    t_service = float(duty["service_time_s"])
    c_limit = float(duty["outlet_mole_fraction_limit"]) * c_in
    t_obs = np.linspace(0.0, t_service, 97)

    def max_outlet(L):
        c_out = simulate_outlet(
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
            also_mid=False,
        )
        return float(np.max(c_out))

    lo, hi = 0.05, 0.12
    if max_outlet(lo) <= c_limit:
        while lo > 1e-3 and max_outlet(lo) <= c_limit:
            hi = lo
            lo *= 0.5
    else:
        while max_outlet(hi) > c_limit:
            lo = hi
            hi *= 1.55
            if hi > 30.0:
                raise RuntimeError("Unable to find feasible packed length")
    for _ in range(22):
        mid = 0.5 * (lo + hi)
        if max_outlet(mid) > c_limit:
            lo = mid
        else:
            hi = mid
    # Safety against MTZ under-prediction from noisy ID.
    return float(hi * 1.05)


def solve_case(case_path: Path) -> dict:
    case, frames = load_case(case_path)
    params = estimate_parameters(case, frames)
    L_min = design_length(params, case)
    return {
        "case_id": case.get("case_id", case_path.parent.name),
        "qm": params["qm"],
        "K": params["K"],
        "k_LDF": params["k_LDF"],
        "D_ax": params["D_ax"],
        "L_min_m": L_min,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    result = solve_case(Path(args.case))
    for key in ("qm", "K", "k_LDF", "D_ax", "L_min_m"):
        if not math.isfinite(result[key]) or result[key] <= 0:
            raise SystemExit(f"non-finite or non-positive {key}")
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
