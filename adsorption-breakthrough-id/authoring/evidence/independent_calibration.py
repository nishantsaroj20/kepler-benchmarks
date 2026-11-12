#!/usr/bin/env python3
"""Independent tolerance check: alternate MOL grid feasibility for oracle L_min.

Runs solution/column_solver.py on each hidden case, records measured parameter
errors and L_min / outlet margins, and confirms feasibility under spatial grids
that differ from the generator — so the 1.30 oversize ceiling and 1.03 outlet
slack are tied to measured evidence rather than an arbitrary mesh.
"""

from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "tests" / "hidden" / "manifest.json"
SOLVER = ROOT / "solution" / "column_solver.py"
ESTIMATOR_REPORT = Path(__file__).with_name("independent_estimator_report.json")


def q_star(c, qm, K):
    return qm * K * c / (1.0 + K * c)


def max_outlet(case, params, L, n_nodes):
    duty = case["design_duty"]
    eps = case["bed_voidage"]
    rho = case["particle_density_kg_m3"]
    c_in = duty["c_in_mol_m3"]
    u = duty["interstitial_velocity_m_s"]
    t_service = duty["service_time_s"]
    dz = L / n_nodes
    factor = ((1.0 - eps) / eps) * rho
    y0 = np.zeros(2 * n_nodes)

    def rhs(_t, y):
        c = y[:n_nodes]
        q = y[n_nodes:]
        dq = params["k_LDF"] * (q_star(c, params["qm"], params["K"]) - q)
        dc = np.empty(n_nodes)
        dc[0] = u * (c[0] - c_in) / dz
        dc[1:] = u * (c[1:] - c[:-1]) / dz
        cim = np.empty(n_nodes)
        cip = np.empty(n_nodes)
        cim[0] = c_in
        cim[1:] = c[:-1]
        cip[:-1] = c[1:]
        cip[-1] = c[-1]
        d2 = (cip - 2.0 * c + cim) / (dz * dz)
        return np.concatenate([-dc + params["D_ax"] * d2 - factor * dq, dq])

    sol = solve_ivp(
        rhs,
        (0.0, t_service),
        y0,
        method="BDF",
        t_eval=np.linspace(0.0, t_service, 181),
        rtol=1e-6,
        atol=1e-9,
    )
    assert sol.success, sol.message
    return float(np.max(sol.y[n_nodes - 1, :]))


def main():
    manifest = json.loads(MANIFEST.read_text())
    global_tol = manifest["tolerances"]
    oversize_ceil = float(global_tol["L_min_oversize_ceil"])
    slack = float(global_tol["outlet_feasibility_slack"])
    rows = []
    for entry in manifest["hidden_cases"]:
        case_id = entry["case_id"]
        case_path = ROOT / "tests" / "hidden" / case_id / "case.json"
        out_path = Path("/tmp") / f"ind_{case_id}.json"
        print("calibrating", case_id, flush=True)
        t0 = time.monotonic()
        proc = subprocess.run(
            [sys.executable, str(SOLVER), "--case", str(case_path), "--out", str(out_path)],
            capture_output=True,
            text=True,
            env={
                "OMP_NUM_THREADS": "1",
                "OPENBLAS_NUM_THREADS": "1",
                "MKL_NUM_THREADS": "1",
                "NUMEXPR_NUM_THREADS": "1",
                "PATH": str(Path(sys.executable).parent) + ":/usr/bin:/bin",
            },
        )
        elapsed = time.monotonic() - t0
        assert proc.returncode == 0, proc.stderr or proc.stdout
        payload = json.loads(out_path.read_text())
        case = json.loads(case_path.read_text())
        tol = entry["tolerances"]
        rel_errors = {}
        param_ok = True
        for key, rel in (
            ("qm", tol["qm_rel"]),
            ("K", tol["K_rel"]),
            ("k_LDF", tol["k_LDF_rel"]),
            ("D_ax", tol["D_ax_rel"]),
        ):
            err = abs(payload[key] - entry["parameters"][key]) / abs(entry["parameters"][key])
            rel_errors[key] = round(float(err), 6)
            param_ok = param_ok and err <= rel
        c_limit = (
            case["design_duty"]["outlet_mole_fraction_limit"]
            * case["design_duty"]["c_in_mol_m3"]
        )
        # Independent grids (48 and 80) — neither needs to match generator n_nodes.
        m48 = max_outlet(case, entry["parameters"], payload["L_min_m"], 48)
        m80 = max_outlet(case, entry["parameters"], payload["L_min_m"], 80)
        case_slack = float(tol["outlet_feasibility_slack"])
        case_ceil = float(tol["L_min_oversize_ceil"])
        feasible = m48 <= c_limit * case_slack and m80 <= c_limit * case_slack
        L_ref = float(entry["L_min_m"])
        L_est = float(payload["L_min_m"])
        L_ratio = L_est / L_ref
        not_huge = L_est <= case_ceil * L_ref
        rows.append(
            {
                "case_id": case_id,
                "elapsed_s": round(elapsed, 2),
                "oracle_estimate": {
                    "qm": payload["qm"],
                    "K": payload["K"],
                    "k_LDF": payload["k_LDF"],
                    "D_ax": payload["D_ax"],
                    "L_min_m": L_est,
                },
                "relative_errors": rel_errors,
                "L_min_m_ref": L_ref,
                "L_min_oversize_ratio": round(L_ratio, 6),
                "param_ok": param_ok,
                "feasible_alt_grids": feasible,
                "not_trivially_oversized": not_huge,
                "max_outlet_n48": m48,
                "max_outlet_n80": m80,
                "outlet_over_limit_n48": round(m48 / c_limit, 6),
                "outlet_over_limit_n80": round(m80 / c_limit, 6),
                "c_limit": c_limit,
                "pass": param_ok and feasible and not_huge,
            }
        )

    max_L_ratio = max(r["L_min_oversize_ratio"] for r in rows)
    max_outlet_frac = max(
        max(r["outlet_over_limit_n48"], r["outlet_over_limit_n80"]) for r in rows
    )

    estimator_cite = None
    if ESTIMATOR_REPORT.is_file():
        est = json.loads(ESTIMATOR_REPORT.read_text())
        estimator_cite = {
            "source": "authoring/evidence/independent_estimator_report.json",
            "all_hidden_pass": est.get("all_hidden_pass"),
            "relative_errors": {
                row["case_id"]: row["relative_errors"] for row in est.get("rows", [])
            },
        }

    report = {
        "schema": "independent_calibration_report.v2",
        "oracle_independent_pass": all(r["pass"] for r in rows),
        "tolerances": {
            "L_min_oversize_ceil": oversize_ceil,
            "outlet_feasibility_slack": slack,
            "qm_rel": global_tol["qm_rel"],
            "K_rel": global_tol["K_rel"],
            "k_LDF_rel": global_tol["k_LDF_rel"],
            "D_ax_rel": global_tol["D_ax_rel"],
        },
        "calibration": {
            "L_min_oversize_ceil": {
                "value": oversize_ceil,
                "meaning": (
                    "Reject L_min larger than ceil × reference design length so "
                    "trivially padded beds fail while still allowing modest safety "
                    "factors and coarser MOL grids."
                ),
                "why": (
                    "Oracle and independent estimators routinely land a few percent "
                    "above the generator binary-search L_min because of 1.05 design "
                    "margins and grid differences; ceil=1.30 leaves headroom (~30%) "
                    "without accepting giant padded answers."
                ),
                "measured_max_oracle_L_ratio": round(max_L_ratio, 6),
                "headroom_vs_ceil": round(oversize_ceil - max_L_ratio, 6),
            },
            "outlet_feasibility_slack": {
                "value": slack,
                "meaning": (
                    "On the true parameters, max outlet at the reported L_min must "
                    "be ≤ slack × purity limit under an independent MOL check."
                ),
                "why": (
                    "Alternate spatial grids (n=48, n=80) differ slightly from the "
                    "generator mesh; slack=1.03 (3%) absorbs that numerical scatter "
                    "without relaxing the purity duty in a material way."
                ),
                "measured_max_outlet_over_limit": round(max_outlet_frac, 6),
                "headroom_vs_slack": round(slack - max_outlet_frac, 6),
                "alt_grids": [48, 80],
            },
        },
        "independent_estimator": estimator_cite,
        "rows": rows,
        "n_rows": len(rows),
        "visible_condition_coverage_note": (
            "Every hidden condition class appears in environment/data/sample_cases/"
        ),
    }
    Path(__file__).with_name("independent_calibration_report.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    print(
        json.dumps(
            {
                "oracle_independent_pass": report["oracle_independent_pass"],
                "n_rows": report["n_rows"],
                "max_L_ratio": max_L_ratio,
                "max_outlet_frac": max_outlet_frac,
            },
            indent=2,
        )
    )
    if not report["oracle_independent_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
