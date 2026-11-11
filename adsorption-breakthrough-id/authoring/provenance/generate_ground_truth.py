#!/usr/bin/env python3
"""Seeded MOL generator for 1D Langmuir+LDF+axial-dispersion breakthrough data.

Writes ground-truth JSON (provenance + tests), noisy CSVs, and bed specs.
Does not place hidden parameters in environment specs.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "environment" / "data"
PROV_DIR = ROOT / "authoring" / "provenance"
TESTS_DIR = ROOT / "tests"

SEED = 20261111

# Hidden isotherm / transport parameters (never shipped in environment specs).
GT = {
    "qm": 0.42,  # mol/kg
    "K": 2.8,  # m^3/mol
    "k_LDF": 0.028,  # 1/s — slow enough that MTZ width is kinetic
    "D_ax": 2.5e-4,  # m^2/s — Pe ~ O(10-50), dispersion visible
}

# Known bed packing constants.
KNOWN = {
    "bed_voidage": 0.40,
    "particle_density_kg_m3": 800.0,
    "lab_bed_length_m": 0.30,
}

EXPERIMENT_A = {
    "experiment_id": "A",
    "c_in_mol_m3": 0.45,
    "superficial_velocity_m_s": 0.035,
    "t_end_s": 3200.0,
}

EXPERIMENT_B = {
    "experiment_id": "B",
    "c_in_mol_m3": 1.40,
    "superficial_velocity_m_s": 0.070,
    "t_end_s": 1400.0,
}

DESIGN = {
    "c_in_mol_m3": 0.85,
    "superficial_velocity_m_s": 0.050,
    "service_time_s": 1000.0,
    "outlet_mole_fraction_limit": 0.05,
}

N_NODES_TRUTH = 72
NOISE_REL = 0.028
NOISE_ABS = 0.003


def q_star(c: np.ndarray | float, qm: float, K: float) -> np.ndarray | float:
    return qm * K * c / (1.0 + K * c)


def simulate_bed(
    *,
    L: float,
    u: float,
    c_in: float,
    t_end: float,
    qm: float,
    K: float,
    k_ldf: float,
    D_ax: float,
    eps: float,
    rho_p: float,
    n_nodes: int = N_NODES_TRUTH,
    t_eval: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Integrate mobile + solid phases with method of lines.

    Inlet: Dirichlet c = c_in. Outlet: zero Neumann. Convection: first-order
    upwind. Dispersion: second-order central with ghost cells.
    """
    dz = L / n_nodes
    factor = ((1.0 - eps) / eps) * rho_p
    y0 = np.zeros(2 * n_nodes, dtype=float)

    def rhs(_t: float, y: np.ndarray) -> np.ndarray:
        c = y[:n_nodes]
        q = y[n_nodes:]
        dqdt = k_ldf * (q_star(c, qm, K) - q)

        dc_conv = np.empty(n_nodes)
        # Upwind with Dirichlet inlet ghost.
        c_left = c_in
        dc_conv[0] = u * (c[0] - c_left) / dz
        dc_conv[1:] = u * (c[1:] - c[:-1]) / dz

        # Axial dispersion with Dirichlet inlet ghost and Neumann outlet ghost.
        c_im1 = np.empty(n_nodes)
        c_ip1 = np.empty(n_nodes)
        c_im1[0] = c_in
        c_im1[1:] = c[:-1]
        c_ip1[:-1] = c[1:]
        c_ip1[-1] = c[-1]
        d2c = (c_ip1 - 2.0 * c + c_im1) / (dz * dz)

        dcdt = -dc_conv + D_ax * d2c - factor * dqdt
        return np.concatenate([dcdt, dqdt])

    if t_eval is None:
        t_eval = np.linspace(0.0, t_end, 241)

    sol = solve_ivp(
        rhs,
        (0.0, float(t_end)),
        y0,
        method="BDF",
        t_eval=t_eval,
        rtol=1e-6,
        atol=1e-9,
    )
    if not sol.success:
        raise RuntimeError(f"MOL integration failed: {sol.message}")

    c_hist = sol.y[:n_nodes, :]
    mid_idx = n_nodes // 2
    c_outlet = c_hist[-1, :]
    c_mid = c_hist[mid_idx, :]
    return sol.t, c_outlet, c_mid


def sample_noisy_series(
    rng: np.random.Generator,
    t: np.ndarray,
    c_outlet: np.ndarray,
    c_mid: np.ndarray,
    c_in: float,
    n_samples: int,
) -> pd.DataFrame:
    # Irregular sparse times, denser near expected front.
    weights = np.linspace(0.35, 1.65, t.size)
    weights /= weights.sum()
    idx = np.sort(rng.choice(np.arange(t.size), size=n_samples, replace=False, p=weights))
    ts = t[idx]
    out = c_outlet[idx]
    mid = c_mid[idx]
    out_noisy = out + rng.normal(0.0, NOISE_REL * c_in + NOISE_ABS, size=out.size)
    mid_noisy = mid + rng.normal(0.0, NOISE_REL * c_in + NOISE_ABS, size=mid.size)
    out_noisy = np.clip(out_noisy, 0.0, None)
    mid_noisy = np.clip(mid_noisy, 0.0, None)
    return pd.DataFrame(
        {
            "time_s": ts,
            "c_outlet_mol_m3": out_noisy,
            "c_midbed_mol_m3": mid_noisy,
        }
    )


def find_L_min(
    *,
    qm: float,
    K: float,
    k_ldf: float,
    D_ax: float,
    eps: float,
    rho_p: float,
    c_in: float,
    u: float,
    t_service: float,
    frac_limit: float,
    n_nodes: int = 48,
) -> float:
    c_limit = frac_limit * c_in

    def violates(L: float) -> bool:
        t_eval = np.linspace(0.0, t_service, 181)
        _t, c_out, _mid = simulate_bed(
            L=L,
            u=u,
            c_in=c_in,
            t_end=t_service,
            qm=qm,
            K=K,
            k_ldf=k_ldf,
            D_ax=D_ax,
            eps=eps,
            rho_p=rho_p,
            n_nodes=n_nodes,
            t_eval=t_eval,
        )
        return float(np.max(c_out)) > c_limit + 1e-12

    # Bracket: start from a small guess and grow until feasible.
    lo = 0.05
    hi = 0.10
    if not violates(lo):
        # Even tiny bed works — shrink downward.
        while lo > 1e-3 and not violates(lo):
            hi = lo
            lo *= 0.5
    else:
        while violates(hi):
            hi *= 1.6
            if hi > 20.0:
                raise RuntimeError("Could not bracket feasible design length")
        lo = hi / 1.6

    for _ in range(28):
        mid = 0.5 * (lo + hi)
        if violates(mid):
            lo = mid
        else:
            hi = mid
    return float(hi)


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    PROV_DIR.mkdir(parents=True, exist_ok=True)
    TESTS_DIR.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(SEED)
    eps = KNOWN["bed_voidage"]
    rho_p = KNOWN["particle_density_kg_m3"]
    L_lab = KNOWN["lab_bed_length_m"]

    frames = {}
    for exp, n_samp in ((EXPERIMENT_A, 55), (EXPERIMENT_B, 48)):
        t_dense = np.linspace(0.0, exp["t_end_s"], 401)
        t, c_out, c_mid = simulate_bed(
            L=L_lab,
            u=exp["superficial_velocity_m_s"],
            c_in=exp["c_in_mol_m3"],
            t_end=exp["t_end_s"],
            qm=GT["qm"],
            K=GT["K"],
            k_ldf=GT["k_LDF"],
            D_ax=GT["D_ax"],
            eps=eps,
            rho_p=rho_p,
            n_nodes=N_NODES_TRUTH,
            t_eval=t_dense,
        )
        frames[exp["experiment_id"]] = sample_noisy_series(
            rng, t, c_out, c_mid, exp["c_in_mol_m3"], n_samp
        )

    L_min = find_L_min(
        qm=GT["qm"],
        K=GT["K"],
        k_ldf=GT["k_LDF"],
        D_ax=GT["D_ax"],
        eps=eps,
        rho_p=rho_p,
        c_in=DESIGN["c_in_mol_m3"],
        u=DESIGN["superficial_velocity_m_s"],
        t_service=DESIGN["service_time_s"],
        frac_limit=DESIGN["outlet_mole_fraction_limit"],
        n_nodes=56,
    )

    # Confirm feasible and near-infeasible margins for evidence.
    t_chk = np.linspace(0.0, DESIGN["service_time_s"], 201)
    _t, c_ok, _ = simulate_bed(
        L=L_min,
        u=DESIGN["superficial_velocity_m_s"],
        c_in=DESIGN["c_in_mol_m3"],
        t_end=DESIGN["service_time_s"],
        qm=GT["qm"],
        K=GT["K"],
        k_ldf=GT["k_LDF"],
        D_ax=GT["D_ax"],
        eps=eps,
        rho_p=rho_p,
        n_nodes=56,
        t_eval=t_chk,
    )
    c_limit = DESIGN["outlet_mole_fraction_limit"] * DESIGN["c_in_mol_m3"]
    max_ok = float(np.max(c_ok))

    gt_payload = {
        "seed": SEED,
        "parameters": {
            "qm": GT["qm"],
            "K": GT["K"],
            "k_LDF": GT["k_LDF"],
            "D_ax": GT["D_ax"],
        },
        "design": {
            "L_min_m": L_min,
            "c_limit_mol_m3": c_limit,
            "max_outlet_at_L_min": max_ok,
        },
        "tolerances": {
            "qm_rel": 0.10,
            "K_rel": 0.12,
            "k_LDF_rel": 0.20,
            "D_ax_rel": 0.30,
            "L_min_rel": 0.10,
        },
        "numerics": {
            "n_nodes_truth": N_NODES_TRUTH,
            "integrator": "BDF",
            "noise_rel": NOISE_REL,
            "noise_abs": NOISE_ABS,
        },
    }

    for path in (PROV_DIR / "ground_truth_params.json", TESTS_DIR / "ground_truth_params.json"):
        path.write_text(json.dumps(gt_payload, indent=2) + "\n")

    frames["A"].to_csv(DATA_DIR / "breakthrough_exp_A.csv", index=False)
    frames["B"].to_csv(DATA_DIR / "breakthrough_exp_B.csv", index=False)

    specs = {
        "bed_voidage": KNOWN["bed_voidage"],
        "particle_density_kg_m3": KNOWN["particle_density_kg_m3"],
        "lab_bed_length_m": KNOWN["lab_bed_length_m"],
        "axial_ports": ["outlet", "midbed"],
        "experiment_A": {
            "csv": "/workspace/data/breakthrough_exp_A.csv",
            "c_in_mol_m3": EXPERIMENT_A["c_in_mol_m3"],
            "superficial_velocity_m_s": EXPERIMENT_A["superficial_velocity_m_s"],
            "t_end_s": EXPERIMENT_A["t_end_s"],
        },
        "experiment_B": {
            "csv": "/workspace/data/breakthrough_exp_B.csv",
            "c_in_mol_m3": EXPERIMENT_B["c_in_mol_m3"],
            "superficial_velocity_m_s": EXPERIMENT_B["superficial_velocity_m_s"],
            "t_end_s": EXPERIMENT_B["t_end_s"],
        },
        "design_duty": {
            "c_in_mol_m3": DESIGN["c_in_mol_m3"],
            "superficial_velocity_m_s": DESIGN["superficial_velocity_m_s"],
            "service_time_s": DESIGN["service_time_s"],
            "outlet_mole_fraction_limit": DESIGN["outlet_mole_fraction_limit"],
        },
        "model_notes": {
            "isotherm": "Langmuir",
            "kinetics": "linear_driving_force",
            "dispersion": "axial_Fickian",
            "inlet_bc": "Dirichlet",
            "outlet_bc": "Neumann_zero_gradient",
        },
    }
    (DATA_DIR / "bed_specs.json").write_text(json.dumps(specs, indent=2) + "\n")

    print(json.dumps({"L_min_m": L_min, "max_outlet_at_L_min": max_ok, "c_limit": c_limit}, indent=2))
    print("Wrote data, specs, and ground_truth_params.json")


if __name__ == "__main__":
    main()
