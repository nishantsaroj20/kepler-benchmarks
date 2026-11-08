"""Centralized comparators that are not admissible architectures.

sb10ad: SLICOT full-order H-infinity output feedback, one controller per case.
State-feedback LMI: optimistic lower bound on the achievable gamma for one case.
"""

from __future__ import annotations

import numpy as np
import cvxpy as cp
from slycot.synthesis import sb10ad


def full_order_hinf(plant: dict, gamma0: float = 5.0) -> dict:
    A = np.array(plant["A"], float, order="F")
    Bw = np.array(plant["Bw"], float)
    Bu = np.array(plant["Bu"], float)
    Cz = np.array(plant["Cz"], float)
    Cy = np.array(plant["Cy"], float)
    Dzw = np.array(plant["Dzw"], float)
    Dzu = np.array(plant["Dzu"], float)
    Dyw = np.array(plant["Dyw"], float)
    Dyu = np.array(plant["Dyu"], float)
    B = np.asfortranarray(np.hstack([Bw, Bu]))
    C = np.asfortranarray(np.vstack([Cz, Cy]))
    D = np.asfortranarray(np.vstack([np.hstack([Dzw, Dzu]), np.hstack([Dyw, Dyu])]))
    n = A.shape[0]
    m = B.shape[1]
    np_ = C.shape[0]
    ncon = Bu.shape[1]
    nmeas = Cy.shape[0]
    try:
        gamma, Ak, Bk, Ck, Dk, *_rest = sb10ad(
            n, m, np_, ncon, nmeas, gamma0, A, B, C, D, job=3, gtol=1e-3, actol=1e-6
        )
    except Exception as exc:
        return {"status": "failed", "error": str(exc), "gamma_reported": None}
    return {
        "status": "returned",
        "gamma_reported": float(gamma),
        "Ak": np.array(Ak, float),
        "Bk": np.array(Bk, float),
        "Ck": np.array(Ck, float),
        "Dk": np.array(Dk, float),
        "n_controller": int(np.array(Ak).shape[0]),
    }


def state_feedback_feasible(plant: dict, gamma: float, eps: float = 1e-6) -> dict:
    A = np.asarray(plant["A"], float)
    B1 = np.asarray(plant["Bw"], float)
    B2 = np.asarray(plant["Bu"], float)
    C1 = np.asarray(plant["Cz"], float)
    D11 = np.asarray(plant["Dzw"], float)
    D12 = np.asarray(plant["Dzu"], float)
    n = A.shape[0]
    m1 = B1.shape[1]
    p1 = C1.shape[0]
    P = cp.Variable((n, n), symmetric=True)
    Y = cp.Variable((B2.shape[1], n))
    M11 = A @ P + P @ A.T + B2 @ Y + Y.T @ B2.T
    M12 = B1
    M13 = P @ C1.T + Y.T @ D12.T
    top = cp.hstack([M11, M12, M13])
    mid = cp.hstack([B1.T, -gamma * np.eye(m1), D11.T])
    bot = cp.hstack([C1 @ P + D12 @ Y, D11, -gamma * np.eye(p1)])
    M = cp.vstack([top, mid, bot])
    prob = cp.Problem(cp.Minimize(0.0), [M << -eps * np.eye(n + m1 + p1), P >> eps * np.eye(n)])
    try:
        prob.solve(solver=cp.CLARABEL, verbose=False)
    except Exception as exc:
        return {"feasible": False, "status": f"error:{exc}"}
    return {"feasible": prob.status in {"optimal", "optimal_inaccurate"}, "status": str(prob.status)}


def state_feedback_gamma(plant: dict) -> dict:
    """Coarse bisection. A finite value is an optimistic lower bound, not a controller we ship."""
    lo, hi = 0.05, 0.2
    feas_hi = None
    for _ in range(8):
        trial = state_feedback_feasible(plant, hi)
        if trial["feasible"]:
            feas_hi = hi
            break
        hi *= 2.0
        if hi > 200:
            return {"status": "no_upper", "last": trial}
    assert feas_hi is not None
    for _ in range(8):
        mid = 0.5 * (lo + hi)
        trial = state_feedback_feasible(plant, mid)
        if trial["feasible"]:
            hi = mid
        else:
            lo = mid
    return {"status": "bracket", "lower": lo, "upper": hi, "note": "state_feedback_optimistic_lower_bound"}
