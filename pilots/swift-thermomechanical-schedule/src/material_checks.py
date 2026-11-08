"""Material-point and conduction checks that do not use the coupled FE residual.

The return map implemented here is a third, scalar implementation used only for
these checks. It is not imported by either finite-element solver.
"""
from __future__ import annotations

import math

import numpy as np

from problem_spec import ALPHA_TH, CHI, E_MPA, EDOT0, L_HALF, NU, RHO_C, T_GRIP, T_REF, TM, TRUTH


def lame():
    mu = E_MPA / (2.0 * (1.0 + NU))
    lam = E_MPA * NU / ((1.0 + NU) * (1.0 - 2.0 * NU))
    kappa = lam + 2.0 * mu / 3.0
    return mu, kappa


def sig_h(ep, card=TRUTH):
    sy0 = card["K"] * card["eps0"] ** card["n"]
    swift = card["K"] * (card["eps0"] + ep) ** card["n"]
    voce = sy0 + card["Q"] * (1.0 - math.exp(-card["beta"] * ep))
    return card["alpha"] * swift + (1.0 - card["alpha"]) * voce


def flow(ep, rate, T, card=TRUTH):
    rate_fac = 1.0 + card["C"] * math.log(max(rate / EDOT0, 1.0))
    theta = min(max((T - T_REF) / (TM - T_REF), 0.0), 1.0)
    return sig_h(ep, card) * rate_fac * (1.0 - theta ** card["m"])


def return_dg(q_trial, ep, T, dt, adiabatic, card=TRUTH):
    """Scalar bisection for q - 3 mu dg - Y = 0. Independent of both FE codes."""
    mu, _ = lame()
    lo, hi = 0.0, q_trial / (3.0 * mu)
    if flow(ep, 0.0, T, card) >= q_trial:
        return 0.0, T, flow(ep, 0.0, T, card)

    def residual(dg):
        rate = dg / dt if dt > 0.0 else 0.0
        if not adiabatic:
            yv = flow(ep + dg, rate, T, card)
            return q_trial - 3.0 * mu * dg - yv, T
        # implicit local heating, fixed-point inside the residual
        temp = T
        yv = flow(ep + dg, rate, temp, card)
        for _ in range(12):
            yv = flow(ep + dg, rate, temp, card)
            temp = T + CHI * yv * dg / RHO_C
        return q_trial - 3.0 * mu * dg - yv, temp

    t_out = T
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        f, t_out = residual(mid)
        if f > 0.0:
            lo = mid
        else:
            hi = mid
    dg = 0.5 * (lo + hi)
    f, t_out = residual(dg)
    return dg, t_out, flow(ep + dg, dg / dt if dt > 0.0 else 0.0, t_out, card)


def check_elastic_uniaxial():
    """Small-strain uniaxial Kirchhoff stress equals E * log-stretch, lateral stress ~ 0."""
    mu, kappa = lame()
    stretch = 1.002
    # Solve lateral stretch so tau_rr = tau_tt = 0.
    lo, hi = 0.99, 1.0
    for _ in range(40):
        lam_r = 0.5 * (lo + hi)
        e_z = math.log(stretch)
        e_r = math.log(lam_r)
        tr = e_z + 2.0 * e_r
        tau_r = kappa * tr + 2.0 * mu * (e_r - tr / 3.0)
        if tau_r > 0.0:
            hi = lam_r
        else:
            lo = lam_r
    lam_r = 0.5 * (lo + hi)
    e_z = math.log(stretch)
    e_r = math.log(lam_r)
    tr = e_z + 2 * e_r
    tau_z = kappa * tr + 2 * mu * (e_z - tr / 3.0)
    tau_r = kappa * tr + 2 * mu * (e_r - tr / 3.0)
    e_eng = stretch - 1.0
    return {
        "tau_z_MPa": tau_z,
        "tau_r_MPa": tau_r,
        "E_times_eng_MPa": E_MPA * e_eng,
        "rel_vs_E_eng": (tau_z - E_MPA * e_eng) / (E_MPA * e_eng),
    }


def check_plastic_consistency():
    mu, _ = lame()
    # Large elastic trial, then return. End von Mises must match the flow stress.
    q_trial = 900.0
    dt = 0.05
    dg, T, yv = return_dg(q_trial, 0.05, 320.0, dt, adiabatic=True)
    q_end = q_trial - 3.0 * mu * dg
    dT = CHI * yv * dg / RHO_C
    return {
        "dg": dg,
        "T_K": T,
        "Y_MPa": yv,
        "q_end_MPa": q_end,
        "abs_consistency_MPa": abs(q_end - yv),
        "adiabatic_dT_from_Y_K": dT,
        "temperature_matches_Y": abs((T - 320.0) - dT),
    }


def check_incompressibility_and_objectivity():
    """Isochoric return preserves det(Fp)=1. A rigid rotation rotates Kirchhoff stress."""
    mu, kappa = lame()
    # Prescribed F with det != 1, plastic correction on the deviator only.
    F = np.array([[1.25, 0.05, 0.0], [0.02, 0.85, 0.0], [0.0, 0.0, 0.94]], dtype=float)
    # Polar-ish: use be = F F^T as a purely elastic trial (Cpinv = I).
    be = F @ F.T
    w, v = np.linalg.eigh(0.5 * (be + be.T))
    e = 0.5 * np.log(np.maximum(w, 1e-16))
    tr = float(e.sum())
    dev = e - tr / 3.0
    q = math.sqrt(1.5) * 2.0 * mu * math.sqrt(float(dev @ dev))
    dg, _, _ = return_dg(q, 0.0, T_REF, 0.1, adiabatic=False)
    dev2 = dev * (1.0 - 3.0 * mu * dg / q)
    e2 = dev2 + tr / 3.0
    # det Fe = exp(tr e) unchanged by the deviatoric corrector, so det Fp = det F / det Fe
    det_fe = math.exp(tr)
    det_fp = abs(np.linalg.det(F)) / det_fe
    # rebuilt be determinant
    be2 = (v * np.exp(2.0 * e2)) @ v.T
    det_be = float(np.linalg.det(be2))
    tau_p = kappa * tr + 2.0 * mu * dev2
    tau = (v * tau_p) @ v.T
    # Objectivity: rotate the elastic Finger tensor.
    ang = 0.4
    Q = np.array([[math.cos(ang), -math.sin(ang), 0.0], [math.sin(ang), math.cos(ang), 0.0], [0.0, 0.0, 1.0]])
    be_r = Q @ be2 @ Q.T
    wr, vr = np.linalg.eigh(0.5 * (be_r + be_r.T))
    er = 0.5 * np.log(np.maximum(wr, 1e-16))
    trr = float(er.sum())
    devr = er - trr / 3.0
    tau_pr = kappa * trr + 2.0 * mu * devr
    tau_rot_computed = (vr * tau_pr) @ vr.T
    tau_rot_expected = Q @ tau @ Q.T
    return {
        "det_Fp": det_fp,
        "det_be_over_detFe2": det_be / det_fe ** 2,
        "objectivity_abs_MPa": float(np.max(np.abs(tau_rot_computed - tau_rot_expected))),
        "dg": dg,
    }


def series_temperature(z, t, t_init, t_grip=T_GRIP, L=L_HALF, n_terms=40):
    """1D axial conduction, Neumann at z=0, Dirichlet at z=L, uniform initial temperature."""
    theta0 = t_init - t_grip
    acc = np.zeros_like(np.asarray(z, dtype=float))
    z = np.asarray(z, dtype=float)
    for n in range(n_terms):
        lam = (n + 0.5) * math.pi / L
        coef = theta0 * (2.0 / L) * ((-1.0) ** n) / lam
        acc = acc + coef * np.cos(lam * z) * math.exp(-ALPHA_TH * lam * lam * t)
    return acc + t_grip


def neck_average_series(t, t_init, z_max=1.5, n_z=40):
    z = np.linspace(0.0, z_max, n_z)
    # reference-volume weight for a nearly uniform bar ~ R(z)^2, approximated as uniform
    # because the imperfection is 0.5%. The check compares this idealisation to the FE.
    val = series_temperature(z, t, t_init)
    return float(np.trapezoid(val, z) / z_max)
