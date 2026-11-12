#!/usr/bin/env python3
"""Independent tolerance check: alternate MOL grid feasibility for oracle L_min.

Confirms that reference parameter recoveries meet stated relative tolerances and
that reported L_min values are feasible under a different spatial grid than the
generator used for truth — so grading is not tied to one proprietary mesh.
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
    rows = []
    for entry in manifest["hidden_cases"]:
        case_id = entry["case_id"]
        case_path = ROOT / "tests" / "hidden" / case_id / "case.json"
        out_path = Path("/tmp") / f"ind_{case_id}.json"
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
        assert proc.returncode == 0, proc.stderr
        payload = json.loads(out_path.read_text())
        case = json.loads(case_path.read_text())
        tol = entry["tolerances"]
        param_ok = True
        for key, rel in (
            ("qm", tol["qm_rel"]),
            ("K", tol["K_rel"]),
            ("k_LDF", tol["k_LDF_rel"]),
            ("D_ax", tol["D_ax_rel"]),
        ):
            err = abs(payload[key] - entry["parameters"][key]) / abs(entry["parameters"][key])
            param_ok = param_ok and err <= rel
        c_limit = (
            case["design_duty"]["outlet_mole_fraction_limit"]
            * case["design_duty"]["c_in_mol_m3"]
        )
        # Independent grids (48 and 80) — neither needs to match generator n_nodes.
        m48 = max_outlet(case, entry["parameters"], payload["L_min_m"], 48)
        m80 = max_outlet(case, entry["parameters"], payload["L_min_m"], 80)
        slack = tol["outlet_feasibility_slack"]
        feasible = m48 <= c_limit * slack and m80 <= c_limit * slack
        not_huge = payload["L_min_m"] <= tol["L_min_oversize_ceil"] * entry["L_min_m"]
        rows.append(
            {
                "case_id": case_id,
                "elapsed_s": round(elapsed, 2),
                "param_ok": param_ok,
                "feasible_alt_grids": feasible,
                "not_trivially_oversized": not_huge,
                "max_outlet_n48": m48,
                "max_outlet_n80": m80,
                "c_limit": c_limit,
                "pass": param_ok and feasible and not_huge,
            }
        )
    report = {
        "oracle_independent_pass": all(r["pass"] for r in rows),
        "rows": rows,
        "visible_condition_coverage_note": (
            "Every hidden condition class appears in environment/data/sample_cases/"
        ),
    }
    Path(__file__).with_name("independent_calibration_report.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
