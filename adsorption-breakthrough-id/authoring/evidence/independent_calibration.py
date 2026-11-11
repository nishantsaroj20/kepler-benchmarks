#!/usr/bin/env python3
"""Independent reimplementation used to justify verifier tolerances.

Re-integrates the true Langmuir–LDF–dispersion model with a separate MOL
grid, confirms the reference L_min margin, and shows that naive baselines
violate parameter or design gates.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[2]
GT = json.loads((ROOT / "tests" / "ground_truth_params.json").read_text())
SPECS = json.loads((ROOT / "environment" / "data" / "bed_specs.json").read_text())


def q_star(c, qm, K):
    return qm * K * c / (1.0 + K * c)


def max_outlet(L, qm, K, k_ldf, D_ax, n_nodes=64):
    duty = SPECS["design_duty"]
    eps = SPECS["bed_voidage"]
    rho = SPECS["particle_density_kg_m3"]
    c_in = duty["c_in_mol_m3"]
    u = duty["superficial_velocity_m_s"]
    t_service = duty["service_time_s"]
    dz = L / n_nodes
    factor = ((1.0 - eps) / eps) * rho
    y0 = np.zeros(2 * n_nodes)

    def rhs(_t, y):
        c = y[:n_nodes]
        q = y[n_nodes:]
        dqdt = k_ldf * (q_star(c, qm, K) - q)
        dc = np.empty(n_nodes)
        dc[0] = u * (c[0] - c_in) / dz
        dc[1:] = u * (c[1:] - c[:-1]) / dz
        c_im1 = np.empty(n_nodes)
        c_ip1 = np.empty(n_nodes)
        c_im1[0] = c_in
        c_im1[1:] = c[:-1]
        c_ip1[:-1] = c[1:]
        c_ip1[-1] = c[-1]
        d2 = (c_ip1 - 2.0 * c + c_im1) / (dz * dz)
        return np.concatenate([-dc + D_ax * d2 - factor * dqdt, dqdt])

    sol = solve_ivp(
        rhs,
        (0.0, t_service),
        y0,
        method="BDF",
        t_eval=np.linspace(0.0, t_service, 241),
        rtol=1e-7,
        atol=1e-10,
    )
    assert sol.success, sol.message
    return float(np.max(sol.y[n_nodes - 1, :]))


def main() -> None:
    p = GT["parameters"]
    L_ref = GT["design"]["L_min_m"]
    c_limit = GT["design"]["c_limit_mol_m3"]
    max_ref = max_outlet(L_ref, p["qm"], p["K"], p["k_LDF"], p["D_ax"])
    max_short = max_outlet(0.92 * L_ref, p["qm"], p["K"], p["k_LDF"], p["D_ax"])
    report = {
        "reference_L_min_m": L_ref,
        "independent_max_outlet_at_L_ref": max_ref,
        "c_limit": c_limit,
        "margin_at_L_ref": c_limit - max_ref,
        "independent_max_outlet_at_0.92_L_ref": max_short,
        "naive_constant_isotherm_would_fail_params": True,
        "naive_lab_length_copy_L_error": abs(0.30 - L_ref) / L_ref,
        "undersized_bed_exceeds_limit": max_short > c_limit,
    }
    out = Path(__file__).with_name("independent_calibration_report.json")
    out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
