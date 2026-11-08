"""Build the public generalized plants from the averaged circuit."""

from __future__ import annotations

import hashlib
import json

import numpy as np

from .equilibrium import solve_equilibrium
from .indexing import IOF, VF, base
from .jacobian_analytic import linearize_analytic
from .params import (
    A_NOISE_PU,
    A_REF_PU,
    A_V_PU,
    I_BASE,
    N_X_PHYS,
    S_RATED,
    V_BASE,
    W_E_KINF,
    W_E_POLE,
    W_E_ZERO,
    W_NAMES,
    W_U_KINF,
    W_U_POLE,
    W_U_ZERO,
    WEIGHT_STATE_NAMES,
    Y_NAMES,
    Z_NAMES,
    U_NAMES,
    Case,
    scaling_diagonal,
    state_names_si,
)
from .rhs_complex import rhs_complex
from .rhs_real import powers_from_filters, rhs_real


def output_matrices(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Cy such that y = Cy x_deviation, with y = (schedule - measured)/S_rated.

    The map does not depend on u. Dyu is identically zero by this construction.
    """
    Cy = np.zeros((4, N_X_PHYS))
    for conv in (0, 1):
        b = base(conv)
        vd = float(x[b + VF])
        vq = float(x[b + VF + 1])
        iod = float(x[b + IOF])
        ioq = float(x[b + IOF + 1])
        # dP = 1.5 (vd dio_d + iod dv_d + vq dio_q + ioq dv_q)
        Cy[2 * conv, b + VF] = -1.5 * iod / S_RATED
        Cy[2 * conv, b + VF + 1] = -1.5 * ioq / S_RATED
        Cy[2 * conv, b + IOF] = -1.5 * vd / S_RATED
        Cy[2 * conv, b + IOF + 1] = -1.5 * vq / S_RATED
        # dQ = 1.5 (vd dio_q + ioq dv_d - vq dio_d - iod dv_q)
        Cy[2 * conv + 1, b + VF] = -1.5 * ioq / S_RATED
        Cy[2 * conv + 1, b + VF + 1] = 1.5 * iod / S_RATED
        Cy[2 * conv + 1, b + IOF] = 1.5 * vq / S_RATED
        Cy[2 * conv + 1, b + IOF + 1] = -1.5 * vd / S_RATED
    Dyu = np.zeros((4, 4))
    return Cy, Dyu


def _dyw() -> np.ndarray:
    Dyw = np.zeros((4, 10))
    for i in range(4):
        Dyw[i, 2 + i] = A_REF_PU
        Dyw[i, 6 + i] = -A_NOISE_PU
    return Dyw


def build_generalized(case: Case) -> dict:
    sol = solve_equilibrium(case)
    if not sol["success"]:
        raise RuntimeError(f"equilibrium failed for {case.case_id}: {sol['message']}")
    x = sol["x"]
    u_amp = sol["u_amp"]
    vg = sol["vg"]
    A_si, Bu_amp, Bv = linearize_analytic(x, u_amp, vg, case)
    Cy_si, Dyu = output_matrices(x)
    T = np.diag(scaling_diagonal())
    Tinv = np.diag(1.0 / np.array(scaling_diagonal()))
    A = Tinv @ A_si @ T
    Bu = Tinv @ Bu_amp * I_BASE
    Bw_v = Tinv @ Bv * (A_V_PU * V_BASE)
    Cy = Cy_si @ T
    Dyw = _dyw()
    n_w = 10
    Bw = np.zeros((N_X_PHYS, n_w))
    Bw[:, 0:2] = Bw_v

    # Append weighting filters. Error weights see y; effort weights see u.
    n_e = 4
    n_u = 4
    n = N_X_PHYS + n_e + n_u
    Ag = np.zeros((n, n))
    Ag[:N_X_PHYS, :N_X_PHYS] = A
    Ag[N_X_PHYS : N_X_PHYS + n_e, :N_X_PHYS] = Cy
    Ag[N_X_PHYS : N_X_PHYS + n_e, N_X_PHYS : N_X_PHYS + n_e] = -W_E_POLE * np.eye(n_e)
    Ag[N_X_PHYS + n_e :, N_X_PHYS + n_e :] = -W_U_POLE * np.eye(n_u)

    Bwg = np.zeros((n, n_w))
    Bwg[:N_X_PHYS, :] = Bw
    Bwg[N_X_PHYS : N_X_PHYS + n_e, :] = Dyw
    Bug = np.zeros((n, n_u))
    Bug[:N_X_PHYS, :] = Bu
    Bug[N_X_PHYS + n_e :, :] = np.eye(n_u)

    Cyg = np.zeros((4, n))
    Cyg[:, :N_X_PHYS] = Cy
    Dyug = np.zeros((4, 4))

    k_e_state = W_E_KINF * (W_E_ZERO - W_E_POLE)
    k_u_state = W_U_KINF * (W_U_ZERO - W_U_POLE)
    Cz = np.zeros((8, n))
    Cz[0:4, :N_X_PHYS] = W_E_KINF * Cy
    Cz[0:4, N_X_PHYS : N_X_PHYS + n_e] = k_e_state * np.eye(n_e)
    Cz[4:8, N_X_PHYS + n_e :] = k_u_state * np.eye(n_u)
    Dzw = np.zeros((8, n_w))
    Dzw[0:4, :] = W_E_KINF * Dyw
    Dzu = np.zeros((8, n_u))
    Dzu[4:8, :] = W_U_KINF * np.eye(n_u)

    p_meas, q_meas = powers_from_filters(x)
    plant = {
        "case_id": case.case_id,
        "role": case.role,
        "A": Ag,
        "Bw": Bwg,
        "Bu": Bug,
        "Cy": Cyg,
        "Dyw": Dyw,
        "Dyu": Dyug,
        "Cz": Cz,
        "Dzw": Dzw,
        "Dzu": Dzu,
        "physical": {
            "A_si": A_si,
            "Bu_amp": Bu_amp,
            "Bv": Bv,
            "Cy_si": Cy_si,
            "A_scaled_physical": A,
            "Bu_pu": Bu,
            "Cy_scaled": Cy,
        },
        "equilibrium": sol,
        "scheduled_P_W": p_meas,
        "scheduled_Q_var": q_meas,
        "u_bias_ampere": u_amp,
        "u_bias_pu": u_amp / I_BASE,
        "state_names": state_names_si() + WEIGHT_STATE_NAMES,
        "w_names": list(W_NAMES),
        "u_names": list(U_NAMES),
        "y_names": list(Y_NAMES),
        "z_names": list(Z_NAMES),
        "n_physical": N_X_PHYS,
        "n_weight": n_e + n_u,
    }
    return plant


def plant_public_dict(plant: dict) -> dict:
    """JSON-ready authoritative generalized plant, scaled contract."""
    eq = plant["equilibrium"]
    bal = eq["balance"]
    pub = {
        "case_id": plant["case_id"],
        "role": plant["role"],
        "plant_set_id": "crc-two-gfl-lcl-v1",
        "coordinate": "scaled_similarity_of_SI_plus_weight_states",
        "state_names": plant["state_names"],
        "w_names": plant["w_names"],
        "u_names": plant["u_names"],
        "y_names": plant["y_names"],
        "z_names": plant["z_names"],
        "n_physical_states": plant["n_physical"],
        "n_weight_states": plant["n_weight"],
        "feedback_sign": (
            "y = (power schedule - filtered measured power)/S_rated. "
            "Positive y requests more injection. u is a pu current-reference perturbation "
            "about the published bias, added before the fixed inner current loop."
        ),
        "Dyu_structurally_zero": True,
        "matrices": {
            key: plant[key].tolist()
            for key in ("A", "Bw", "Bu", "Cy", "Dyw", "Dyu", "Cz", "Dzw", "Dzu")
        },
        "equilibrium": {
            "rhs_inf_si": eq["rhs_inf_si"],
            "p_W": eq["p_W"].tolist(),
            "q_var": eq["q_var"].tolist(),
            "u_bias_pu": plant["u_bias_pu"].tolist(),
            "u_bias_ampere": plant["u_bias_ampere"].tolist(),
            "delta_rad": bal["delta_rad"],
            "v_t_peak_V": bal["v_t_peak_V"],
            "power_imbalance_W": bal["imbalance_W"],
            "p_loss_W": bal["p_loss_W"],
            "p_grid_absorbed_W": bal["p_grid_absorbed_W"],
            "p_converter_W": bal["p_converter_W"],
            "x_si": eq["x"].tolist(),
            "vg_V": eq["vg"].tolist(),
        },
    }
    blob = json.dumps(pub["matrices"], separators=(",", ":")).encode()
    pub["matrix_sha256"] = hashlib.sha256(blob).hexdigest()
    return pub


def compare_paths(case: Case, plant: dict) -> dict:
    """Finite-difference Jacobian of both residuals versus the analytic matrix."""
    eq = plant["equilibrium"]
    x = eq["x"]
    u = eq["u_amp"]
    vg = eq["vg"]
    A_si = plant["physical"]["A_si"]
    fr = rhs_real(x, u, vg, case)
    fc = rhs_complex(x, u, vg, case)
    rhs_err = float(np.max(np.abs(fr - fc)))

    def fd(fun, eps: float) -> np.ndarray:
        M = np.zeros_like(A_si)
        for i in range(x.size):
            step = eps * max(1.0, abs(float(x[i])))
            xp = x.copy()
            xm = x.copy()
            xp[i] += step
            xm[i] -= step
            M[:, i] = (fun(xp) - fun(xm)) / (2.0 * step)
        return M

    rels = {}
    for eps in (1e-4, 1e-5, 1e-6, 1e-7):
        fd_real = fd(lambda z: rhs_real(z, u, vg, case), eps)
        rels[f"real_{eps}"] = float(
            np.linalg.norm(A_si - fd_real) / np.linalg.norm(A_si)
        )
    fd_c = fd(lambda z: rhs_complex(z, u, vg, case), 1e-5)
    rel_c = float(np.linalg.norm(A_si - fd_c) / np.linalg.norm(A_si))
    # Realization-invariant check: singular values of y/u at shared frequencies.
    Cy = plant["physical"]["Cy_si"]
    Bu = plant["physical"]["Bu_amp"]
    ws = np.logspace(0, 4, 25)
    sig_err = 0.0
    eye = np.eye(A_si.shape[0])
    # Build FD state matrix from the complex residual and compare FRF.
    for w in ws:
        Gr = Cy @ np.linalg.solve(1j * w * eye - A_si, Bu)
        Gc = Cy @ np.linalg.solve(1j * w * eye - fd_c, Bu)
        sr = np.linalg.svd(Gr, compute_uv=False)
        sc = np.linalg.svd(Gc, compute_uv=False)
        sig_err = max(sig_err, float(np.max(np.abs(sr - sc))))
    return {
        "rhs_inf_difference": rhs_err,
        "fd_relative_by_step": rels,
        "fd_complex_rel_1e-5": rel_c,
        "sigma_max_abs_diff": sig_err,
    }
