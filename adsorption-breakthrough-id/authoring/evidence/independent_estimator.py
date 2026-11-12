#!/usr/bin/env python3
"""Independent column solver — structurally different from solution/column_solver.py.

Differences vs the oracle reference:
- Fit grid: N=40 (oracle: 32); design grid N=36 (oracle: 36) but Radau-free BDF
- Parameterization: log10 (oracle: natural log)
- Residual: equal-weight outlet + mid (oracle: front-Gaussian + 0.55 mid)
- Mid filter: MAD spike rule (oracle: hard amplitude clip)
- Starts: capacity seed × different kinetics lattice / stride

Feasibility / grading validation (report) uses the SAME gates as the verifier:
BDF, N=64, 120 s wall clock, relative tolerances, 1.30 L ceil, 1.03 slack.

CLI (verifier-compatible):
  python3 independent_estimator.py --case CASE.json --out OUT.json
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[2]
HIDDEN = ROOT / "tests" / "hidden"
N_FIT = 40
_trapz = getattr(np, "trapezoid", None) or np.trapz


def q_star(c, qm, K):
    return qm * K * c / (1.0 + K * c)


def simulate(
    L, u, c_in, t_obs, qm, K, k_ldf, D_ax, eps, rho, n=N_FIT, also_mid=False, method="BDF"
):
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
    sol = solve_ivp(rhs, (0.0, t_end), y0, method=method, t_eval=te, rtol=3e-5, atol=1e-8)
    if not sol.success:
        bad = np.full_like(t_obs, 1e3, dtype=float)
        return (bad, bad) if also_mid else bad
    c_out = np.interp(t_obs, sol.t, sol.y[n - 1])
    if not also_mid:
        return c_out
    return c_out, np.interp(t_obs, sol.t, sol.y[mid])


def pack4(qm, K, k, D):
    return np.log10(np.array([qm, K, k, D], dtype=float))


def unpack4(theta):
    return tuple(float(x) for x in 10.0 ** np.asarray(theta, dtype=float))


def load_case(case_path: Path):
    case = json.loads(case_path.read_text())
    case_dir = case_path.parent
    frames = {}
    for key in ("experiment_A", "experiment_B"):
        frames[key] = pd.read_csv(case_dir / case[key]["csv"])
    return case, frames


def seed_isotherm(case, frames):
    eps = case["bed_voidage"]
    rho = case["particle_density_kg_m3"]
    L = case["lab_bed_length_m"]
    qs, cs = [], []
    for key in ("experiment_A", "experiment_B"):
        meta = case[key]
        df = frames[key]
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


def mid_mask(y_mid, y_out):
    finite = np.isfinite(y_mid)
    if finite.sum() < 8:
        return finite
    ym = y_mid.copy()
    med = np.nanmedian(ym[finite])
    mad = np.nanmedian(np.abs(ym[finite] - med)) + 1e-12
    ok = finite.copy()
    for i in np.where(finite)[0]:
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
    for key in ("experiment_A", "experiment_B"):
        meta = case[key]
        df = frames[key]
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
        chunks.append((pred_out - y_out) / scale)
        if mask.any():
            chunks.append(1.1 * (pred_mid[mask] - y_mid[mask]) / scale)
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
                    n=36,
                    method="BDF",
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


def estimate_case(case_path: Path) -> dict:
    case, frames = load_case(case_path)
    qm0, K0 = seed_isotherm(case, frames)
    starts = []
    for qm in (qm0, qm0 * 0.9, qm0 * 1.1):
        for K in (K0, K0 * 0.7, K0 * 1.3, max(0.12, K0 * 0.45)):
            for k_ldf, D_ax in (
                (0.02, 1.8e-4),
                (0.03, 3.5e-4),
                (0.045, 7.0e-4),
                (0.055, 1.1e-3),
            ):
                starts.append(pack4(qm, max(K, 0.1), k_ldf, D_ax))
    starts = starts[::5][:8]
    lo = pack4(0.08, 0.08, 0.005, 3e-5)
    hi = pack4(2.5, 20.0, 0.25, 2.5e-3)
    best = None
    for th0 in starts:
        res = least_squares(
            residual,
            th0,
            args=(case, frames),
            bounds=(lo, hi),
            xtol=3e-7,
            ftol=3e-7,
            gtol=3e-7,
            max_nfev=22,
        )
        cost = float(np.sum(res.fun**2))
        if best is None or cost < best[0]:
            best = (cost, res.x)
    res = least_squares(
        residual,
        best[1],
        args=(case, frames),
        bounds=(lo, hi),
        xtol=1e-8,
        ftol=1e-8,
        gtol=1e-8,
        max_nfev=35,
    )
    qm, K, k_ldf, D_ax = unpack4(res.x)
    params = {"qm": qm, "K": K, "k_LDF": k_ldf, "D_ax": D_ax}
    out = {
        "case_id": case.get("case_id"),
        **params,
        "L_min_m": design_L(params, case),
    }
    return out


def main_cli():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    payload = estimate_case(Path(args.case))
    Path(args.out).write_text(json.dumps(payload, indent=2) + "\n")


def main_report():
    """Run independent solver through the real verifier isolation + BDF/N=64 gates."""
    # Import verifier helpers (same gates / isolation as pytest).
    import sys

    sys.path.insert(0, str(ROOT / "tests"))
    import test_evaluation as te  # noqa: E402

    manifest = json.loads((HIDDEN / "manifest.json").read_text())
    timeout_s = float(manifest["per_run_timeout_s"])
    program = Path(__file__).resolve()
    rows = []
    all_ok = True
    for entry in manifest["hidden_cases"]:
        case_id = entry["case_id"]
        case_src = HIDDEN / case_id
        case = json.loads((case_src / "case.json").read_text())
        print("verifier-path", case_id, flush=True)
        t0 = time.monotonic()
        try:
            run = te.run_program_isolated(program, case_src, timeout_s)
            elapsed = run["elapsed_s"]
            payload = run["payload"]
            ok, detail = te.answer_passes_gates(payload, entry, case)
            under_budget = elapsed <= timeout_s
            ok = ok and under_budget
        except BaseException as exc:
            elapsed = time.monotonic() - t0
            payload = {}
            detail = {"error": str(exc)}
            ok = False
            under_budget = elapsed <= timeout_s
        all_ok = all_ok and ok
        rows.append(
            {
                "case_id": case_id,
                "elapsed_s": round(elapsed, 3),
                "under_120s": under_budget,
                "estimate": payload,
                "relative_errors": {
                    k: round(v, 6) for k, v in detail.get("relative_errors", {}).items()
                },
                "feasible_bdf_n64": detail.get("feasible"),
                "not_oversized": detail.get("not_oversized"),
                "max_outlet_bdf_n64": detail.get("max_outlet"),
                "c_limit": detail.get("c_limit"),
                "pass": ok,
            }
        )
        print(
            case_id,
            "PASS" if ok else "FAIL",
            "t=",
            round(elapsed, 2),
            detail.get("relative_errors"),
            flush=True,
        )
    report = {
        "estimator": (
            "independent_estimator.py (BDF N=40 equal-weight outlet+mid, log10 "
            "least_squares, MAD mid filter; ≠ oracle BDF N=32 front-weighted ln-space). "
            "Validation uses the real verifier path: isolated subprocess + BDF/N=64 "
            "feasibility gates + 120s wall clock."
        ),
        "gates": {
            "integrator": "BDF",
            "n_nodes": 64,
            "per_run_timeout_s": timeout_s,
            "L_min_oversize_ceil": manifest["tolerances"]["L_min_oversize_ceil"],
            "outlet_feasibility_slack": manifest["tolerances"]["outlet_feasibility_slack"],
        },
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
    import sys

    if len(sys.argv) > 1 and sys.argv[1] == "--report":
        sys.argv = [sys.argv[0]] + sys.argv[2:]
        main_report()
    elif "--case" in sys.argv:
        main_cli()
    else:
        main_report()
