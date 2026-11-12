#!/usr/bin/env python3
"""Independent estimator structurally different from solution/column_solver.py.

Differences vs the oracle reference:
- Integrator: LSODA (oracle: BDF)
- Spatial mesh: N=48 (oracle: 32 fit / 36 design)
- Parameterization: log10 space (oracle: natural log)
- Residual: equal-weight outlet + mid-bed (oracle: front-Gaussian outlet weights
  with 0.55 mid scale); mid points filtered by MAD spike rule (oracle: hard
  amplitude clip)
- Search: coarse log-space grid search then short polish (oracle: multi-start
  ln-space least_squares only)
- Design: LSODA N=40 binary search (oracle: BDF N=36)

Produces authoring/evidence/independent_estimator_report.json with measured errors.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[2]
HIDDEN = ROOT / "tests" / "hidden"
N = 48
_trapz = getattr(np, "trapezoid", None) or np.trapz


def q_star(c, qm, K):
    return qm * K * c / (1.0 + K * c)


def simulate(L, u, c_in, t_obs, qm, K, k_ldf, D_ax, eps, rho, n=N, also_mid=False):
    dz = L / n
    fac = ((1.0 - eps) / eps) * rho
    y0 = np.zeros(2 * n)
    mid = n // 2
    t_obs = np.asarray(t_obs, dtype=float)
    t_end = float(np.max(t_obs)) if t_obs.size else 0.0
    if t_end <= 0.0:
        z = np.zeros_like(t_obs)
        return (z, z) if also_mid else z

    def rhs(_t, y):
        c = y[:n]
        q = y[n:]
        dq = k_ldf * (q_star(c, qm, K) - q)
        dc_conv = np.empty(n)
        dc_conv[0] = u * (c[0] - c_in) / dz
        dc_conv[1:] = u * (c[1:] - c[:-1]) / dz
        c_im1 = np.empty(n)
        c_ip1 = np.empty(n)
        c_im1[0] = c_in
        c_im1[1:] = c[:-1]
        c_ip1[:-1] = c[1:]
        c_ip1[-1] = c[-1]
        d2 = (c_ip1 - 2.0 * c + c_im1) / (dz * dz)
        return np.concatenate([-dc_conv + D_ax * d2 - fac * dq, dq])

    te = np.unique(np.concatenate(([0.0], t_obs)))
    sol = solve_ivp(
        rhs, (0.0, t_end), y0, method="LSODA", t_eval=te, rtol=1e-5, atol=1e-8
    )
    if not sol.success:
        bad = np.full_like(t_obs, 1e3, dtype=float)
        return (bad, bad) if also_mid else bad
    c_out = np.interp(t_obs, sol.t, sol.y[n - 1])
    if not also_mid:
        return c_out
    c_mid = np.interp(t_obs, sol.t, sol.y[mid])
    return c_out, c_mid


def load_case(case_dir: Path):
    case = json.loads((case_dir / "case.json").read_text())
    frames = {
        "A": pd.read_csv(case_dir / case["experiment_A"]["csv"]),
        "B": pd.read_csv(case_dir / case["experiment_B"]["csv"]),
    }
    return case, frames


def seed_isotherm(case, frames):
    eps = case["bed_voidage"]
    rho = case["particle_density_kg_m3"]
    L = case["lab_bed_length_m"]
    qs, cs = [], []
    for tag, key in (("A", "experiment_A"), ("B", "experiment_B")):
        meta = case[key]
        df = frames[tag]
        c_in = float(meta["c_in_mol_m3"])
        u = float(meta["interstitial_velocity_m_s"])
        t = df["time_s"].to_numpy(float)
        c = df["c_outlet_mol_m3"].to_numpy(float)
        order = np.argsort(t)
        t, c = t[order], c[order]
        unused = _trapz(np.clip(c_in - c, 0.0, None), t)
        rhs = u * unused / L - c_in
        qs.append(max(rhs * eps / ((1.0 - eps) * rho), 1e-4))
        cs.append(c_in)
    c1, c2 = cs
    q1, q2 = qs
    denom = q2 - q1
    if abs(denom) < 1e-12:
        qm, K = max(q1, q2) * 1.2, 1.0
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
                qm, K = max(q1, q2) * 1.2, 1.0
        else:
            qm = q1 * (1.0 + K * c1) / (K * c1)
    return float(np.clip(qm, 0.08, 2.5)), float(np.clip(K, 0.08, 20.0))


def pack4(qm, K, k, D):
    return np.log10(np.array([qm, K, k, D], dtype=float))


def unpack4(theta):
    return tuple(float(x) for x in 10.0 ** np.asarray(theta, dtype=float))


def mid_mask(y_mid, y_out):
    """MAD-based spike/dropout filter (≠ oracle hard amplitude clip)."""
    finite = np.isfinite(y_mid)
    if finite.sum() < 8:
        return finite
    # Smooth-ish baseline via rolling median on finite points.
    ym = y_mid.copy()
    med = np.nanmedian(ym[finite])
    mad = np.nanmedian(np.abs(ym[finite] - med)) + 1e-12
    # Keep points within ~6 MAD of local 5-point median when possible.
    ok = finite.copy()
    idx = np.where(finite)[0]
    for i in idx:
        lo = max(0, i - 2)
        hi = min(len(ym), i + 3)
        local = ym[lo:hi]
        local = local[np.isfinite(local)]
        if local.size == 0:
            ok[i] = False
            continue
        loc_med = np.median(local)
        if abs(ym[i] - loc_med) > 6.0 * mad and abs(ym[i] - loc_med) > 0.15 * max(
            float(np.nanmax(y_out)), 1e-6
        ):
            ok[i] = False
    if ok.sum() < max(6, len(y_mid) // 4):
        return finite
    return ok


def residual(theta, case, frames):
    qm, K, k_ldf, D_ax = unpack4(theta)
    eps = case["bed_voidage"]
    rho = case["particle_density_kg_m3"]
    L = case["lab_bed_length_m"]
    chunks = []
    for tag, key in (("A", "experiment_A"), ("B", "experiment_B")):
        meta = case[key]
        df = frames[tag]
        t = df["time_s"].to_numpy(float)
        y_out = df["c_outlet_mol_m3"].to_numpy(float)
        y_mid = df["c_midbed_mol_m3"].to_numpy(float)
        scale = max(float(meta["c_in_mol_m3"]), 1e-6)
        mask = mid_mask(y_mid, y_out)
        pred_out, pred_mid = simulate(
            L,
            meta["interstitial_velocity_m_s"],
            meta["c_in_mol_m3"],
            t,
            qm,
            K,
            k_ldf,
            D_ax,
            eps,
            rho,
            also_mid=True,
        )
        # Equal weights (≠ oracle front-Gaussian). Mid scaled slightly up for D_ax.
        chunks.append((pred_out - y_out) / scale)
        if mask.any():
            chunks.append(1.15 * (pred_mid[mask] - y_mid[mask]) / scale)
    return np.concatenate(chunks)


def design_L(params, case):
    duty = case["design_duty"]
    eps = case["bed_voidage"]
    rho = case["particle_density_kg_m3"]
    c_in = float(duty["c_in_mol_m3"])
    u = float(duty["interstitial_velocity_m_s"])
    t_service = float(duty["service_time_s"])
    c_limit = float(duty["outlet_mole_fraction_limit"]) * c_in
    t_obs = np.linspace(0.0, t_service, 81)

    def max_out(L):
        return float(
            np.max(
                simulate(
                    L,
                    u,
                    c_in,
                    t_obs,
                    params["qm"],
                    params["K"],
                    params["k_LDF"],
                    params["D_ax"],
                    eps,
                    rho,
                    n=40,
                )
            )
        )

    lo, hi = 0.04, 0.15
    if max_out(lo) <= c_limit:
        while lo > 1e-3 and max_out(lo) <= c_limit:
            hi = lo
            lo *= 0.55
    else:
        while max_out(hi) > c_limit:
            lo = hi
            hi *= 1.6
            if hi > 40:
                raise RuntimeError("no feasible L")
    for _ in range(22):
        mid = 0.5 * (lo + hi)
        if max_out(mid) > c_limit:
            lo = mid
        else:
            hi = mid
    return float(hi * 1.05)


def estimate_case(case_dir: Path) -> dict:
    case, frames = load_case(case_dir)
    qm0, K0 = seed_isotherm(case, frames)
    # Structure ≠ oracle multi-start least_squares: two-stage log-space grid
    # search (isotherm scout → kinetics refine) spanning both sides of the
    # capacity seed, then a short local polish of the winners.
    iso_scouts = []
    for qm in (
        qm0 * 0.97,
        qm0,
        qm0 * 1.08,
        qm0 * 1.10,
        qm0 * 1.12,
        qm0 * 1.18,
        qm0 * 1.25,
    ):
        for K in (
            max(0.12, K0 * 0.55),
            max(0.12, K0 * 0.65),
            max(0.12, K0 * 0.75),
            K0,
            min(20.0, K0 * 1.20),
            min(20.0, K0 * 1.35),
            min(20.0, K0 * 1.55),
        ):
            for k_ldf, D_ax in (
                (0.022, 1.4e-4),
                (0.031, 2.0e-4),
                (0.030, 3.0e-4),
                (0.045, 7.0e-4),
            ):
                th = pack4(qm, max(K, 0.1), k_ldf, D_ax)
                iso_scouts.append(
                    (float(np.sum(residual(th, case, frames) ** 2)), qm, max(K, 0.1))
                )
    iso_scouts.sort(key=lambda item: item[0])
    # Unique (qm,K) survivors.
    isotherms = []
    seen = set()
    for _c, qm, K in iso_scouts:
        key = (round(qm, 4), round(K, 4))
        if key in seen:
            continue
        seen.add(key)
        isotherms.append((qm, K))
        if len(isotherms) >= 6:
            break

    grid = []
    for qm, K in isotherms:
        for k_ldf in (0.018, 0.022, 0.028, 0.031, 0.038, 0.05):
            for D_ax in (1.0e-4, 1.4e-4, 2.0e-4, 2.5e-4, 3.5e-4, 6.0e-4, 9.0e-4):
                grid.append(pack4(qm, K, k_ldf, D_ax))

    scored = []
    for th in grid:
        scored.append((float(np.sum(residual(th, case, frames) ** 2)), th))
    scored.sort(key=lambda item: item[0])

    lo = pack4(0.08, 0.08, 0.005, 3e-5)
    hi = pack4(2.5, 20.0, 0.25, 2.5e-3)
    best = None
    for _cost0, th0 in scored[:5]:
        res = least_squares(
            residual,
            th0,
            args=(case, frames),
            bounds=(lo, hi),
            xtol=1e-7,
            ftol=1e-7,
            gtol=1e-7,
            max_nfev=30,
        )
        cost = float(np.sum(res.fun**2))
        if best is None or cost < best[0]:
            best = (cost, res.x)

    qm, K, k_ldf, D_ax = unpack4(best[1])
    params = {"qm": qm, "K": K, "k_LDF": k_ldf, "D_ax": D_ax}
    return {"case_id": case["case_id"], **params, "L_min_m": design_L(params, case)}


def main():
    manifest = json.loads((HIDDEN / "manifest.json").read_text())
    rows = []
    all_ok = True
    for entry in manifest["hidden_cases"]:
        case_id = entry["case_id"]
        print("estimating", case_id, flush=True)
        est = estimate_case(HIDDEN / case_id)
        tol = entry["tolerances"]
        errs = {}
        ok = True
        for key, rel in (
            ("qm", tol["qm_rel"]),
            ("K", tol["K_rel"]),
            ("k_LDF", tol["k_LDF_rel"]),
            ("D_ax", tol["D_ax_rel"]),
        ):
            err = abs(est[key] - entry["parameters"][key]) / abs(entry["parameters"][key])
            errs[key] = err
            ok = ok and err <= rel
        case = json.loads((HIDDEN / case_id / "case.json").read_text())
        c_in = case["design_duty"]["c_in_mol_m3"]
        c_limit = case["design_duty"]["outlet_mole_fraction_limit"] * c_in
        max_out = float(
            np.max(
                simulate(
                    est["L_min_m"],
                    case["design_duty"]["interstitial_velocity_m_s"],
                    c_in,
                    np.linspace(0.0, case["design_duty"]["service_time_s"], 121),
                    entry["parameters"]["qm"],
                    entry["parameters"]["K"],
                    entry["parameters"]["k_LDF"],
                    entry["parameters"]["D_ax"],
                    case["bed_voidage"],
                    case["particle_density_kg_m3"],
                    n=56,
                )
            )
        )
        feas = max_out <= c_limit * tol["outlet_feasibility_slack"]
        not_huge = est["L_min_m"] <= tol["L_min_oversize_ceil"] * entry["L_min_m"]
        ok = ok and feas and not_huge
        all_ok = all_ok and ok
        rows.append(
            {
                "case_id": case_id,
                "estimate": est,
                "relative_errors": {k: round(v, 6) for k, v in errs.items()},
                "L_min_m_est": est["L_min_m"],
                "L_min_m_ref": entry["L_min_m"],
                "true_model_max_outlet": max_out,
                "c_limit": c_limit,
                "pass": ok,
            }
        )
        print(
            case_id,
            "PASS" if ok else "FAIL",
            {k: round(v, 4) for k, v in errs.items()},
            flush=True,
        )
    report = {
        "estimator": (
            "independent_estimator.py (LSODA N=48, equal-weight outlet+mid, "
            "log10 grid search + polish, MAD mid filter; ≠ oracle BDF N=32 "
            "front-weighted ln-space multi-start least_squares)"
        ),
        "all_hidden_pass": all_ok,
        "rows": rows,
    }
    Path(__file__).with_name("independent_estimator_report.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    print(json.dumps({"all_hidden_pass": all_ok}, indent=2))
    if not all_ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
