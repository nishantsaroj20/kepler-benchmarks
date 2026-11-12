#!/usr/bin/env python3
"""Generate visible sample cases and deliberately shifted hidden cases.

Hidden shifts are NOT reseeds of the visible set. Each hidden case changes
regime, material, sensor integrity, or operating range (see SHIFT_NOTES).
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[2]
VIS_ROOT = ROOT / "environment" / "data" / "sample_cases"
HID_ROOT = ROOT / "tests" / "hidden"
PROV_DIR = ROOT / "authoring" / "provenance"
EVID_DIR = ROOT / "authoring" / "evidence"

N_NODES = 64
NOISE_REL = 0.025
NOISE_ABS = 0.0025

# Tolerances applied by the verifier to program outputs.
TOL = {
    "qm_rel": 0.12,
    "K_rel": 0.25,
    "k_LDF_rel": 0.50,
    "D_ax_rel": 0.45,
    "L_min_rel": 0.20,
    "L_min_undersize_floor": 0.95,
}

SHIFT_NOTES = {
    "H1_new_material": (
        "New packing material: lower voidage and denser particles than any "
        "visible sample, with a different Langmuir pair."
    ),
    "H2_dispersion_regime": (
        "High superficial velocity / large D_ax dispersion-dominated MTZ; "
        "outside the Pe range of the visible samples."
    ),
    "H3_sensor_faults": (
        "Same broad class as visible kinetics, but mid-bed channel has "
        "injected dropouts and spike faults the solver must ignore/robustify."
    ),
    "H4_high_cin": (
        "Feed concentrations far above the visible range (near-saturation "
        "Langmuir curvature not seen in samples)."
    ),
    "H5_near_linear": (
        "Near-linear isotherm (small K) plus longer service window — edge "
        "case absent from visible duties."
    ),
}


def q_star(c, qm, K):
    return qm * K * c / (1.0 + K * c)


def simulate_bed(
    *,
    L,
    u,
    c_in,
    t_end,
    qm,
    K,
    k_ldf,
    D_ax,
    eps,
    rho_p,
    n_nodes=N_NODES,
    t_eval=None,
):
    dz = L / n_nodes
    factor = ((1.0 - eps) / eps) * rho_p
    y0 = np.zeros(2 * n_nodes)

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

    if t_eval is None:
        t_eval = np.linspace(0.0, t_end, 301)
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
        raise RuntimeError(sol.message)
    mid = n_nodes // 2
    return sol.t, sol.y[n_nodes - 1, :], sol.y[mid, :]


def find_L_min(params, bed, duty, n_nodes=56):
    c_in = duty["c_in_mol_m3"]
    u = duty["superficial_velocity_m_s"]
    t_service = duty["service_time_s"]
    c_limit = duty["outlet_mole_fraction_limit"] * c_in
    t_eval = np.linspace(0.0, t_service, 181)

    def violates(L):
        _t, c_out, _m = simulate_bed(
            L=L,
            u=u,
            c_in=c_in,
            t_end=t_service,
            qm=params["qm"],
            K=params["K"],
            k_ldf=params["k_LDF"],
            D_ax=params["D_ax"],
            eps=bed["bed_voidage"],
            rho_p=bed["particle_density_kg_m3"],
            n_nodes=n_nodes,
            t_eval=t_eval,
        )
        return float(np.max(c_out)) > c_limit + 1e-12

    lo, hi = 0.04, 0.10
    if not violates(lo):
        while lo > 1e-3 and not violates(lo):
            hi = lo
            lo *= 0.5
    else:
        while violates(hi):
            lo = hi
            hi *= 1.55
            if hi > 30:
                raise RuntimeError("cannot bracket L_min")
        # lo is violating (or the last violating expand step), hi is feasible
    for _ in range(28):
        mid = 0.5 * (lo + hi)
        if violates(mid):
            lo = mid
        else:
            hi = mid
    return float(hi)


def sample_csv(rng, t, c_out, c_mid, c_in, n_samples, fault_mode=None):
    weights = np.linspace(0.4, 1.6, t.size)
    weights /= weights.sum()
    idx = np.sort(rng.choice(np.arange(t.size), size=n_samples, replace=False, p=weights))
    ts = t[idx]
    out = c_out[idx] + rng.normal(0.0, NOISE_REL * c_in + NOISE_ABS, size=n_samples)
    mid = c_mid[idx] + rng.normal(0.0, NOISE_REL * c_in + NOISE_ABS, size=n_samples)
    out = np.clip(out, 0.0, None)
    mid = np.clip(mid, 0.0, None)
    if fault_mode == "midbed_spikes_dropouts":
        # Inject ~15% dropouts (NaN) and a few large spikes on mid-bed only.
        n_drop = max(3, n_samples // 7)
        drop_idx = rng.choice(n_samples, size=n_drop, replace=False)
        mid = mid.astype(float)
        mid[drop_idx] = np.nan
        spike_idx = rng.choice([i for i in range(n_samples) if i not in set(drop_idx)], size=3, replace=False)
        mid[spike_idx] = mid[spike_idx] + 0.55 * c_in + rng.uniform(0.1, 0.3, size=3)
    return pd.DataFrame(
        {
            "time_s": ts,
            "c_outlet_mol_m3": out,
            "c_midbed_mol_m3": mid,
        }
    )


def stoichiometric_time(bed, params, exp):
    """Equilibrium-front time consistent with the interstitial PDE form used here."""
    eps = bed["bed_voidage"]
    rho = bed["particle_density_kg_m3"]
    L = bed["lab_bed_length_m"]
    c_in = exp["c_in_mol_m3"]
    u = exp["superficial_velocity_m_s"]  # advection speed u in the PDE
    qs = q_star(c_in, params["qm"], params["K"])
    return (L / u) * (1.0 + ((1.0 - eps) / eps) * rho * qs / c_in)


def simulate_until_breakthrough(bed, params, exp, target_frac=0.92):
    """Grow the horizon until outlet reaches target_frac of c_in."""
    exp = dict(exp)
    t_st = stoichiometric_time(bed, params, exp)
    t_end = float(max(exp.get("t_end_s", 0.0), 1.5 * t_st))
    c_in = exp["c_in_mol_m3"]
    last = None
    for _ in range(8):
        t_dense = np.linspace(0.0, t_end, 361)
        t, c_out, c_mid = simulate_bed(
            L=bed["lab_bed_length_m"],
            u=exp["superficial_velocity_m_s"],
            c_in=c_in,
            t_end=t_end,
            qm=params["qm"],
            K=params["K"],
            k_ldf=params["k_LDF"],
            D_ax=params["D_ax"],
            eps=bed["bed_voidage"],
            rho_p=bed["particle_density_kg_m3"],
            t_eval=t_dense,
        )
        last = (t, c_out, c_mid, t_end)
        if float(np.max(c_out)) >= target_frac * c_in:
            exp["t_end_s"] = float(t_end)
            return exp, t, c_out, c_mid
        t_end *= 1.7
    raise RuntimeError(
        f"incomplete breakthrough after growth: max={float(np.max(last[1]))} "
        f"cin={c_in} t_end={last[3]}"
    )


def write_case(case_dir: Path, case_id: str, bed, params, exp_a, exp_b, duty, rng, fault_mode=None):
    case_dir.mkdir(parents=True, exist_ok=True)
    frames = {}
    resolved = {}
    for tag, exp, n_samp in (("A", exp_a, 52), ("B", exp_b, 48)):
        exp_res, t, c_out, c_mid = simulate_until_breakthrough(bed, params, exp)
        resolved[tag] = exp_res
        frames[tag] = sample_csv(
            rng,
            t,
            c_out,
            c_mid,
            exp_res["c_in_mol_m3"],
            n_samp,
            fault_mode=fault_mode,
        )
    exp_a, exp_b = resolved["A"], resolved["B"]

    frames["A"].to_csv(case_dir / "breakthrough_exp_A.csv", index=False)
    frames["B"].to_csv(case_dir / "breakthrough_exp_B.csv", index=False)

    L_min = find_L_min(params, bed, duty)
    case = {
        "case_id": case_id,
        "bed_voidage": bed["bed_voidage"],
        "particle_density_kg_m3": bed["particle_density_kg_m3"],
        "lab_bed_length_m": bed["lab_bed_length_m"],
        "experiment_A": {
            "csv": "breakthrough_exp_A.csv",
            "c_in_mol_m3": exp_a["c_in_mol_m3"],
            "superficial_velocity_m_s": exp_a["superficial_velocity_m_s"],
            "t_end_s": exp_a["t_end_s"],
        },
        "experiment_B": {
            "csv": "breakthrough_exp_B.csv",
            "c_in_mol_m3": exp_b["c_in_mol_m3"],
            "superficial_velocity_m_s": exp_b["superficial_velocity_m_s"],
            "t_end_s": exp_b["t_end_s"],
        },
        "design_duty": duty,
        "model_notes": {
            "isotherm": "Langmuir",
            "kinetics": "linear_driving_force",
            "dispersion": "axial_Fickian",
            "inlet_bc": "Dirichlet",
            "outlet_bc": "Neumann_zero_gradient",
        },
    }
    (case_dir / "case.json").write_text(json.dumps(case, indent=2) + "\n")
    truth = {
        "case_id": case_id,
        "parameters": params,
        "L_min_m": L_min,
        "tolerances": TOL,
    }
    return truth


def main():
    if VIS_ROOT.exists():
        shutil.rmtree(VIS_ROOT)
    if HID_ROOT.exists():
        shutil.rmtree(HID_ROOT)
    VIS_ROOT.mkdir(parents=True)
    HID_ROOT.mkdir(parents=True)
    PROV_DIR.mkdir(parents=True, exist_ok=True)
    EVID_DIR.mkdir(parents=True, exist_ok=True)

    # ----- Visible samples (agent-facing) -----
    visible_specs = []
    visible_truth_for_provenance = []

    vis_cases = [
        {
            "case_id": "V1_base",
            "seed": 20261112,
            "bed": {"bed_voidage": 0.40, "particle_density_kg_m3": 800.0, "lab_bed_length_m": 0.30},
            "params": {"qm": 0.42, "K": 2.8, "k_LDF": 0.028, "D_ax": 2.5e-4},
            "exp_a": {"c_in_mol_m3": 0.45, "superficial_velocity_m_s": 0.035, "t_end_s": 3200.0},
            "exp_b": {"c_in_mol_m3": 1.40, "superficial_velocity_m_s": 0.070, "t_end_s": 1400.0},
            "duty": {
                "c_in_mol_m3": 0.85,
                "superficial_velocity_m_s": 0.050,
                "service_time_s": 1000.0,
                "outlet_mole_fraction_limit": 0.05,
            },
        },
        {
            "case_id": "V2_moderate_u",
            "seed": 20261113,
            "bed": {"bed_voidage": 0.42, "particle_density_kg_m3": 760.0, "lab_bed_length_m": 0.28},
            "params": {"qm": 0.38, "K": 2.2, "k_LDF": 0.035, "D_ax": 1.8e-4},
            "exp_a": {"c_in_mol_m3": 0.55, "superficial_velocity_m_s": 0.040, "t_end_s": 2800.0},
            "exp_b": {"c_in_mol_m3": 1.20, "superficial_velocity_m_s": 0.065, "t_end_s": 1200.0},
            "duty": {
                "c_in_mol_m3": 0.75,
                "superficial_velocity_m_s": 0.045,
                "service_time_s": 900.0,
                "outlet_mole_fraction_limit": 0.05,
            },
        },
    ]

    for spec in vis_cases:
        rng = np.random.default_rng(spec["seed"])
        d = VIS_ROOT / spec["case_id"]
        truth = write_case(
            d,
            spec["case_id"],
            spec["bed"],
            spec["params"],
            spec["exp_a"],
            spec["exp_b"],
            spec["duty"],
            rng,
        )
        visible_specs.append({"case_id": spec["case_id"], "path": f"sample_cases/{spec['case_id']}/case.json"})
        visible_truth_for_provenance.append(truth)

    # Agent-facing index (no hidden answers).
    index = {
        "program_path": "/workspace/output/column_solver.py",
        "per_run_timeout_s": 120,
        "sample_cases": visible_specs,
        "cli": {
            "case": "--case",
            "out": "--out",
            "example": (
                "python3 /workspace/output/column_solver.py "
                "--case /workspace/data/sample_cases/V1_base/case.json "
                "--out /workspace/output/sample_V1_result.json"
            ),
        },
    }
    (ROOT / "environment" / "data" / "run_index.json").write_text(json.dumps(index, indent=2) + "\n")

    # ----- Hidden cases (verifier-only) -----
    hidden = [
        {
            "case_id": "H1_new_material",
            "seed": 91001,
            "bed": {"bed_voidage": 0.32, "particle_density_kg_m3": 1250.0, "lab_bed_length_m": 0.26},
            "params": {"qm": 0.55, "K": 3.6, "k_LDF": 0.022, "D_ax": 1.4e-4},
            "exp_a": {"c_in_mol_m3": 0.50, "superficial_velocity_m_s": 0.030, "t_end_s": 3600.0},
            "exp_b": {"c_in_mol_m3": 1.10, "superficial_velocity_m_s": 0.055, "t_end_s": 1600.0},
            "duty": {
                "c_in_mol_m3": 0.70,
                "superficial_velocity_m_s": 0.040,
                "service_time_s": 1100.0,
                "outlet_mole_fraction_limit": 0.04,
            },
        },
        {
            "case_id": "H2_dispersion_regime",
            "seed": 91002,
            "bed": {"bed_voidage": 0.40, "particle_density_kg_m3": 820.0, "lab_bed_length_m": 0.22},
            "params": {"qm": 0.40, "K": 2.5, "k_LDF": 0.045, "D_ax": 9.0e-4},
            "exp_a": {"c_in_mol_m3": 0.60, "superficial_velocity_m_s": 0.090, "t_end_s": 900.0},
            "exp_b": {"c_in_mol_m3": 1.30, "superficial_velocity_m_s": 0.120, "t_end_s": 700.0},
            "duty": {
                "c_in_mol_m3": 0.95,
                "superficial_velocity_m_s": 0.100,
                "service_time_s": 500.0,
                "outlet_mole_fraction_limit": 0.06,
            },
        },
        {
            "case_id": "H3_sensor_faults",
            "seed": 91003,
            "fault_mode": "midbed_spikes_dropouts",
            "bed": {"bed_voidage": 0.39, "particle_density_kg_m3": 790.0, "lab_bed_length_m": 0.30},
            "params": {"qm": 0.44, "K": 2.9, "k_LDF": 0.026, "D_ax": 2.2e-4},
            "exp_a": {"c_in_mol_m3": 0.48, "superficial_velocity_m_s": 0.036, "t_end_s": 3000.0},
            "exp_b": {"c_in_mol_m3": 1.35, "superficial_velocity_m_s": 0.068, "t_end_s": 1300.0},
            "duty": {
                "c_in_mol_m3": 0.80,
                "superficial_velocity_m_s": 0.048,
                "service_time_s": 950.0,
                "outlet_mole_fraction_limit": 0.05,
            },
        },
        {
            "case_id": "H4_high_cin",
            "seed": 91004,
            "bed": {"bed_voidage": 0.41, "particle_density_kg_m3": 780.0, "lab_bed_length_m": 0.27},
            "params": {"qm": 0.48, "K": 1.9, "k_LDF": 0.031, "D_ax": 2.0e-4},
            # Still outside the visible c_in range (<=1.4), but spaced for Langmuir ID.
            "exp_a": {"c_in_mol_m3": 1.90, "superficial_velocity_m_s": 0.038, "t_end_s": 2200.0},
            "exp_b": {"c_in_mol_m3": 3.40, "superficial_velocity_m_s": 0.072, "t_end_s": 1100.0},
            "duty": {
                "c_in_mol_m3": 2.60,
                "superficial_velocity_m_s": 0.055,
                "service_time_s": 800.0,
                "outlet_mole_fraction_limit": 0.05,
            },
        },
        {
            "case_id": "H5_near_linear",
            "seed": 91005,
            "bed": {"bed_voidage": 0.43, "particle_density_kg_m3": 740.0, "lab_bed_length_m": 0.32},
            "params": {"qm": 0.60, "K": 0.35, "k_LDF": 0.020, "D_ax": 1.6e-4},
            "exp_a": {"c_in_mol_m3": 0.40, "superficial_velocity_m_s": 0.028, "t_end_s": 4000.0},
            "exp_b": {"c_in_mol_m3": 0.90, "superficial_velocity_m_s": 0.050, "t_end_s": 2200.0},
            "duty": {
                "c_in_mol_m3": 0.65,
                "superficial_velocity_m_s": 0.035,
                "service_time_s": 1600.0,
                "outlet_mole_fraction_limit": 0.05,
            },
        },
    ]

    hidden_truths = []
    for spec in hidden:
        rng = np.random.default_rng(spec["seed"])
        d = HID_ROOT / spec["case_id"]
        truth = write_case(
            d,
            spec["case_id"],
            spec["bed"],
            spec["params"],
            spec["exp_a"],
            spec["exp_b"],
            spec["duty"],
            rng,
            fault_mode=spec.get("fault_mode"),
        )
        truth["shift"] = SHIFT_NOTES[spec["case_id"]]
        hidden_truths.append(truth)

    manifest = {
        "per_run_timeout_s": 120,
        "program_artifact": "/workspace/output/column_solver.py",
        "tolerances": TOL,
        "hidden_cases": hidden_truths,
        "shift_summary": SHIFT_NOTES,
        "visible_case_ids": [c["case_id"] for c in vis_cases],
        "notes": (
            "Hidden cases deliberately differ from visible samples via new "
            "material properties, dispersion-dominated regime, injected mid-bed "
            "sensor faults, out-of-range feed concentrations, and a near-linear "
            "isotherm edge case. No visible case appears in the hidden set; "
            "shifts are not reseeds."
        ),
    }
    (HID_ROOT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (PROV_DIR / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (PROV_DIR / "visible_truth.json").write_text(
        json.dumps(visible_truth_for_provenance, indent=2) + "\n"
    )

    shift_doc = {
        "built_by": "authoring/provenance/generate_all_cases.py",
        "hidden_root": "tests/hidden/",
        "manifest": "tests/hidden/manifest.json",
        "shifts": SHIFT_NOTES,
    }
    (EVID_DIR / "hidden_shift_notes.json").write_text(json.dumps(shift_doc, indent=2) + "\n")

    print(json.dumps({
        "visible": [t["case_id"] for t in visible_truth_for_provenance],
        "hidden": [(t["case_id"], round(t["L_min_m"], 4)) for t in hidden_truths],
    }, indent=2))


if __name__ == "__main__":
    main()
