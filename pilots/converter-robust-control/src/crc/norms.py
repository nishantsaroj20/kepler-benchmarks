"""H-infinity / L-infinity checks.

Three numerical routes:

1. SLICOT AB13DD via slycot (Bruinsma-Steinbuch, equilibrated).
   This is an L-infinity norm. It is the H-infinity norm only if the
   system is stable. python-control's linfnorm documentation states the
   same distinction; that package is quarantined on PyPI and is not imported.

2. An independent Hamiltonian bisection (Boyd / Zhou-Doyle-Glover test:
   no imaginary-axis eigenvalue of the Hamiltonian, and sigma_max(D) < gamma).

3. A bounded-real lemma LMI in CVXPY/Clarabel. Feasibility is accepted
   only after the returned P reproduces an LMI whose largest eigenvalue
   is negative. This is an upper certificate, not a norm estimate.

A frequency grid is a lower bound and a search diagnostic. It is never
reported as the H-infinity norm.
"""

from __future__ import annotations

import numpy as np
import cvxpy as cp
from slycot.analysis import ab13dd

from .params import EIG_PERT_REL


def spectral_abscissa(A: np.ndarray) -> np.ndarray:
    return np.linalg.eigvals(np.asarray(A, dtype=float))


def stability_report(A: np.ndarray, rel_pert: float = EIG_PERT_REL) -> dict:
    """Bauer-Fike style guard on the rightmost eigenvalues.

    PASS if abscissa + uncertainty <= pass limit is decided by the caller.
    This function only returns the abscissa and an uncertainty estimate.
    """
    A = np.asarray(A, dtype=float)
    n = A.shape[0]
    eigvals, evecs = np.linalg.eig(A)
    # Eigenvalue condition numbers: 1 / |w* v| for unit-norm left/right vectors.
    try:
        left = np.linalg.inv(evecs)
    except np.linalg.LinAlgError:
        return {
            "abscissa": float(np.max(eigvals.real)),
            "uncertainty": float("inf"),
            "eig_cond_max": float("inf"),
            "n": n,
            "status": "unresolved",
            "reason": "eigenvector matrix is singular",
        }
    conds = []
    for i in range(n):
        dot = abs(left[i, :] @ evecs[:, i])
        conds.append(np.inf if dot < 1e-16 else 1.0 / dot)
    conds_a = np.asarray(conds)
    scale = float(np.linalg.norm(A, ord="fro"))
    # Absolute eigenvalue perturbation bound used as a budget, not a proof.
    uncertainty = float(np.max(conds_a) * rel_pert * max(scale, 1.0))
    abscissa = float(np.max(eigvals.real))
    return {
        "abscissa": abscissa,
        "uncertainty": uncertainty,
        "eig_cond_max": float(np.max(conds_a)),
        "frobenius": scale,
        "n": n,
        "rightmost_hz": float(abs(eigvals[np.argmax(eigvals.real)].imag) / (2.0 * np.pi)),
    }


def classify_stability(report: dict, limit: float) -> str:
    a = report["abscissa"]
    d = report["uncertainty"]
    if not np.isfinite(a) or not np.isfinite(d):
        return "unresolved"
    if a + d <= limit:
        return "pass"
    if a - d > limit:
        return "fail"
    return "unresolved"


def frequency_lower_bound(
    A: np.ndarray,
    B: np.ndarray,
    C: np.ndarray,
    D: np.ndarray,
    n_freq: int = 80,
) -> dict:
    """Max singular value on a fixed grid. A lower bound up to FRF roundoff."""
    A = np.asarray(A, float)
    B = np.asarray(B, float)
    C = np.asarray(C, float)
    D = np.asarray(D, float)
    n = A.shape[0]
    sig_d = float(np.linalg.svd(D, compute_uv=False)[0]) if D.size else 0.0
    if n == 0:
        return {"lower": sig_d, "fpeak_rad_s": float("inf"), "kind": "feedthrough"}
    ws = np.concatenate(
        [
            np.array([0.0, 1e-3, 1e-2]),
            np.logspace(-1, 5, n_freq),
            np.array([1e6]),
        ]
    )
    peak = sig_d
    fpeak = float("inf")
    eye = np.eye(n)
    for w in ws:
        try:
            G = C @ np.linalg.solve(1j * w * eye - A, B) + D
        except np.linalg.LinAlgError:
            continue
        s = float(np.linalg.svd(G, compute_uv=False)[0])
        if s > peak:
            peak = s
            fpeak = float(w)
    return {"lower": peak, "fpeak_rad_s": fpeak, "kind": "frequency_grid_lower_bound"}


def _sigma_at(A, B, C, D, w: float) -> float:
    n = A.shape[0]
    if n == 0:
        return float(np.linalg.svd(D, compute_uv=False)[0])
    G = C @ np.linalg.solve(1j * w * np.eye(n) - A, B) + D
    return float(np.linalg.svd(G, compute_uv=False)[0])


def hamiltonian_imag_hit(A, B, C, D, gamma: float, real_tol: float) -> bool:
    """True when the gamma level is not a strict upper bound.

    sigma_max(D) >= gamma, a failed factorization, or a Hamiltonian
    eigenvalue with |real| below the tolerance counts as a hit.
    """
    D = np.asarray(D, float)
    sig_d = float(np.linalg.svd(D, compute_uv=False)[0]) if D.size else 0.0
    if sig_d >= gamma * (1.0 - 1e-12):
        return True
    n = A.shape[0]
    if n == 0:
        return False
    R = (gamma**2) * np.eye(B.shape[1]) - D.T @ D
    # Keep a little distance from a singular R.
    if np.linalg.eigvalsh((R + R.T) * 0.5)[0] <= 1e-14 * gamma**2:
        return True
    n = A.shape[0]
    Atil = np.block([[A, np.zeros((n, n))], [-C.T @ C, -A.T]])
    Bhat = np.vstack([B, C.T @ D])
    Chat = np.hstack([D.T @ C, B.T])
    H = Atil + Bhat @ np.linalg.solve(R, Chat)
    ev = np.linalg.eigvals(H)
    scale = max(1.0, float(np.linalg.norm(H, ord=np.inf)))
    return bool(np.min(np.abs(ev.real)) <= real_tol * scale)


def refined_frequency_lower_bound(A, B, C, D) -> dict:
    """Lower bound from a grid plus bounded scalar searches on log-frequency.

    This is still not an upper bound. It is the independent peak estimate
    used beside SLICOT. A narrow missed resonance would leave it low.
    """
    from scipy.optimize import minimize_scalar

    grid = frequency_lower_bound(A, B, C, D, n_freq=60)
    best = float(grid["lower"])
    fbest = float(grid["fpeak_rad_s"]) if np.isfinite(grid["fpeak_rad_s"]) else 0.0

    def objective(log_w: float) -> float:
        return -_sigma_at(A, B, C, D, 10.0 ** log_w)

    for lo, hi in ((-4.0, -1.0), (-1.0, 1.0), (1.0, 2.5), (2.5, 4.0), (4.0, 5.5)):
        res = minimize_scalar(objective, bounds=(lo, hi), method="bounded", options={"xatol": 1e-3})
        val = float(-res.fun)
        if val > best:
            best = val
            fbest = float(10.0 ** res.x)
    dc = _sigma_at(A, B, C, D, 0.0)
    if dc > best:
        best = float(dc)
        fbest = 0.0
    return {"lower": best, "fpeak_rad_s": fbest, "grid_lower": float(grid["lower"]), "kind": "refined_frequency_lower_bound"}


def hamiltonian_norm(
    A: np.ndarray,
    B: np.ndarray,
    C: np.ndarray,
    D: np.ndarray,
    tol: float = 1e-6,
    real_tol: float = 1e-8,
    max_iter: int = 48,
) -> dict:
    """Bisection. `upper` is a candidate strict upper bound; `lower` is not above the norm
    if the imaginary-axis test and the frequency samples are accurate.
    """
    A = np.atleast_2d(np.asarray(A, float))
    B = np.atleast_2d(np.asarray(B, float))
    C = np.atleast_2d(np.asarray(C, float))
    D = np.atleast_2d(np.asarray(D, float))
    if B.shape[0] == 1 and A.shape[0] != 1 and B.shape[1] == A.shape[0]:
        # np.atleast_2d on a row vector is already 2d; shapes are checked below.
        pass
    grid = frequency_lower_bound(A, B, C, D)
    lower = float(grid["lower"])
    if lower == 0.0:
        lower = 0.0
    # The imaginary-axis test is trustworthy only if it fires below the
    # suspected peak and clears above it. Otherwise the bracket is inactive.
    probe_lo = max(lower * 0.5, 1e-8)
    probe_hi = max(lower * 1.5, lower + 1e-3)
    if not (
        hamiltonian_imag_hit(A, B, C, D, probe_lo, real_tol)
        and not hamiltonian_imag_hit(A, B, C, D, probe_hi, real_tol)
    ):
        return {
            "lower": lower,
            "upper": None,
            "fpeak_grid": grid["fpeak_rad_s"],
            "grid_lower": lower,
            "status": "imag_axis_test_inactive",
            "iterations": 0,
        }
    # Expand an upper bound.
    upper = probe_hi
    hit = True
    for _ in range(30):
        hit = hamiltonian_imag_hit(A, B, C, D, upper, real_tol)
        if not hit:
            break
        upper = max(2.0 * upper, lower + 1.0)
        if upper > 1e8:
            return {
                "lower": lower,
                "upper": float("inf"),
                "fpeak_grid": grid["fpeak_rad_s"],
                "status": "unbounded_or_failed",
                "iterations": max_iter,
            }
    for it in range(max_iter):
        if upper == float("inf"):
            break
        if (upper - lower) <= tol * max(1.0, upper):
            break
        mid = 0.5 * (lower + upper)
        if hamiltonian_imag_hit(A, B, C, D, mid, real_tol):
            lower = mid
        else:
            upper = mid
    return {
        "lower": float(lower),
        "upper": float(upper),
        "fpeak_grid": grid["fpeak_rad_s"],
        "grid_lower": float(grid["lower"]),
        "status": "bracket",
        "iterations": it + 1,
    }


def slicot_linf(A, B, C, D, tol: float = 1e-8) -> dict:
    A = np.array(A, dtype=float, order="F")
    B = np.array(B, dtype=float, order="F")
    C = np.array(C, dtype=float, order="F")
    D = np.array(D, dtype=float, order="F")
    n = A.shape[0]
    m = B.shape[1]
    p = C.shape[0]
    if n == 0:
        sig = float(np.linalg.svd(D, compute_uv=False)[0]) if D.size else 0.0
        return {"value": sig, "fpeak": float("inf"), "kind": "feedthrough"}
    gpeak, fpeak = ab13dd(
        "C",
        "I",
        "S",
        "D",
        n,
        m,
        p,
        A,
        np.eye(n),
        B,
        C,
        D,
        tol,
    )
    return {"value": float(gpeak), "fpeak": float(fpeak), "kind": "slicot_ab13dd_linf"}


def bounded_real_certificate(A, B, C, D, gamma: float, eps: float = 1e-7) -> dict:
    """Upper certificate: feasible and residual-checked implies ||G||_inf < gamma."""
    A = np.asarray(A, float)
    B = np.asarray(B, float)
    C = np.asarray(C, float)
    D = np.asarray(D, float)
    n = A.shape[0]
    m = B.shape[1]
    p = C.shape[0]
    if n == 0:
        sig = float(np.linalg.svd(D, compute_uv=False)[0]) if D.size else 0.0
        ok = sig < gamma - eps
        return {"feasible": ok, "certified": ok, "status": "static", "lmi_max_eig": sig - gamma}
    P = cp.Variable((n, n), symmetric=True)
    M11 = A.T @ P + P @ A
    M12 = P @ B
    M13 = C.T
    top = cp.hstack([M11, M12, M13])
    mid = cp.hstack([M12.T, -gamma * np.eye(m), D.T])
    bot = cp.hstack([M13.T, D, -gamma * np.eye(p)])
    M = cp.vstack([top, mid, bot])
    prob = cp.Problem(
        cp.Minimize(0.0),
        [M << -eps * np.eye(n + m + p), P >> eps * np.eye(n)],
    )
    try:
        prob.solve(solver=cp.CLARABEL, verbose=False, tol_gap_abs=1e-9, tol_gap_rel=1e-9)
    except Exception as exc:  # solver crash is not a scientific pass
        return {"feasible": False, "certified": False, "status": f"error:{exc}", "lmi_max_eig": None}
    if prob.status not in {"optimal", "optimal_inaccurate"} or P.value is None:
        return {
            "feasible": False,
            "certified": False,
            "status": str(prob.status),
            "lmi_max_eig": None,
        }
    Pv = np.array(P.value, dtype=float)
    Pv = 0.5 * (Pv + Pv.T)
    M11n = A.T @ Pv + Pv @ A
    M12n = Pv @ B
    Mn = np.block(
        [
            [M11n, M12n, C.T],
            [M12n.T, -gamma * np.eye(m), D.T],
            [C, D, -gamma * np.eye(p)],
        ]
    )
    lmi_max = float(np.linalg.eigvalsh(0.5 * (Mn + Mn.T))[-1])
    p_min = float(np.linalg.eigvalsh(Pv)[0])
    certified = lmi_max < -0.1 * eps and p_min > 0.1 * eps
    return {
        "feasible": True,
        "certified": bool(certified),
        "status": str(prob.status),
        "lmi_max_eig": lmi_max,
        "P_min_eig": p_min,
        "gamma": float(gamma),
    }


def evaluate_transfer(A, B, C, D, stab_limit: float) -> dict:
    """Stability first. Norms are H-infinity only when stability passes or fails clearly.

    An unresolved stability result does not get a finite H-infinity pass.
    """
    stab = stability_report(A)
    stab_class = classify_stability(stab, stab_limit)
    grid = frequency_lower_bound(A, B, C, D)
    out = {
        "stability": stab,
        "stability_class": stab_class,
        "grid_lower": grid,
    }
    if stab_class == "fail":
        # Still compute L-infinity so the unstable-vs-Hinf distinction is recorded.
        try:
            linf = slicot_linf(A, B, C, D)
        except Exception as exc:
            linf = {"value": None, "error": str(exc)}
        out["slicot_linf"] = linf
        out["hinf"] = {"status": "infinite", "reason": "unstable"}
        out["pass_hinf"] = False
        return out
    if stab_class == "unresolved":
        out["hinf"] = {"status": "unresolved", "reason": "stability margin unresolved"}
        out["pass_hinf"] = False
        return out
    ham = hamiltonian_norm(A, B, C, D)
    try:
        linf = slicot_linf(A, B, C, D)
        linf_err = None
    except Exception as exc:
        linf = {"value": None}
        linf_err = str(exc)
    out["hamiltonian"] = ham
    out["slicot_linf"] = linf
    out["slicot_error"] = linf_err
    upper = ham.get("upper")
    if linf.get("value") is None or upper is None or not np.isfinite(upper):
        out["hinf"] = {
            "status": "slicot_only" if linf.get("value") is not None else "unverified",
            "reason": linf_err or ham.get("status"),
        }
        out["pass_hinf"] = False
        return out
    mid = 0.5 * (ham["lower"] + upper)
    spread = abs(mid - linf["value"])
    out["method_spread"] = float(spread)
    out["hinf_upper"] = float(max(ham["upper"], linf["value"]))
    out["hinf_lower"] = float(min(ham["lower"], linf["value"], grid["lower"]))
    # SLICOT can sit slightly outside a loose bracket; require both uppers.
    out["pass_hinf_number"] = float(max(ham["upper"], linf["value"]))
    return out
